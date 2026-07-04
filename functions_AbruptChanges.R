fit_deterministic <- function(data){
  
  # Convert raw data into dynr format.
  # This step defines the structure of the panel data:
  # subject ID, time, and observed variable.
  datadynr <- dynr.data(
    data = data,
    id = "id",
    time = "time",
    observed = "X"
  )
  
  # Define how the latent states map to the observed variable.
  # Here, X is measured as a combination of two latent states:
  # Xlevel and Xa.
  meas <- prep.measurement(
    values.load = matrix(c(1, 0), nrow = 1),
    state.names = c("Xlevel", "Xa"),
    obs.names = "X"
  )
  
  # Set initial values for the latent states and their variability.
  # These are starting points for the estimation algorithm.
  initial <- prep.initial(
    values.inistate = c(0, 0),
    params.inistate = c("meanX0", "meanXa"),
    values.inicov = matrix(c(1, 0.1, 0.1, 0.5), 2),
    params.inicov = matrix(c("varX0", "covX0Xa", "covX0Xa", "varXa"), 2)
  )
  
  # Define measurement error and latent process noise.
  # In this deterministic version, latent noise is fixed to zero,
  # meaning that all variability is explained by the dynamics.
  noise <- prep.noise(
    values.latent = matrix(0, 2, 2),
    params.latent = matrix(NA, 2, 2),
    values.observed = matrix(0.1),
    params.observed = "meserr"
  )
  
  # Specify the dynamic model (continuous-time formulation).
  # Xlevel changes over time as a function of itself (beta)
  # and the secondary state Xa.
  # Xa is assumed to be constant.
  dynamics <- prep.formulaDynamics(
    formula = list(
      Xlevel ~ beta * Xlevel + Xa,
      Xa ~ 0
    ),
    startval = c(beta = -0.2),
    isContinuousTime = TRUE
  )
  
  # Combine all model components into a dynr model object.
  model <- dynr.model(dynamics, meas, noise, initial, datadynr)
  
  # Estimate the model parameters using dynr's Kalman filter-based algorithm.
  fit <- dynr.cook(model, verbose = FALSE)
  
  # Return both the fitted model and the model specification.
  list(fit = fit, model = model)
}

detect_shocks <- function(taste, data,
                          method = c("combined", "chi", "t"),
                          p = 0.05,
                          theta_threshold = 1.65,
                          level = "Xlevel"){
  
  method <- match.arg(method)
  
  # Extract diagnostics
  
  delta <- taste$delta.inn[level, ]
  theta <- scale(delta)[, 1]
  
  chi_p <- taste$chi.inn.pval
  t_p   <- taste$t.inn.pval[level, ]
  
  # Selection rule
  
  if(method == "chi"){
    pos <- which(chi_p < p)
  }
  
  if(method == "t"){
    pos <- which(t_p < p)
  }
  
  if(method == "combined"){
    pos <- which(
      chi_p < p &
        abs(theta) > theta_threshold
    )
  }
  
  # Output
  
  shocks <- data[pos, c("id", "time")]
  
  # Boundary filtering: the innovation at the first time point of each
  # subject lacks backward support (purely forward-looking), so it is
  # confounded with initialization effects and excluded.
  keep <- rep(TRUE, nrow(shocks))
  
  for(i in seq_len(nrow(shocks))){
    
    id_i <- shocks$id[i]
    t_i  <- shocks$time[i]
    
    sub_times <- data$time[data$id == id_i]
    
    first_t <- min(sub_times)
    
    if(t_i == first_t){
      keep[i] <- FALSE
    }
  }
  
  shocks <- shocks[keep, , drop = FALSE]
  
  # Metadata (does not change structure)
  attr(shocks, "statistic") <- method
  
  shocks
}

fit_stochastic_with_shocks <- function(data){
  
  # Convert the dataset into dynr format.
  # In this version, we include an additional time-varying covariate (delta_L),
  # which represents the magnitude of detected shocks.
  datadynr <- dynr.data(
    data = data,
    id = "id",
    time = "time",
    observed = "X",
    covariates = "delta_L"
  )
  
  # Define the measurement model.
  # The observed variable X is expressed as a function of two latent states:
  # Xlevel (main developmental process) and Xa (auxiliary component).
  meas <- prep.measurement(
    values.load = matrix(c(1, 0), 1),
    state.names = c("Xlevel", "Xa"),
    obs.names = "X"
  )
  
  # Set initial values for the latent states and their variability.
  # These provide starting points for the estimation algorithm.
  initial <- prep.initial(
    values.inistate = c(0, 0),
    params.inistate = c("meanX0", "meanXa"),
    values.inicov = matrix(c(1, 0.1, 0.1, 0.5), 2),
    params.inicov = matrix(c("varX0", "covX0Xa", "covX0Xa", "varXa"), 2)
  )
  
  # Define measurement error and latent process noise.
  # Here, latent noise is allowed (dynerr), meaning that the system
  # can vary stochastically around its deterministic dynamics.
  noise <- prep.noise(
    values.latent = matrix(c(0.1, 0, 0, 0), 2),
    params.latent = matrix(c("dynerr", NA, NA, NA), 2),
    values.observed = matrix(0.1),
    params.observed = "meserr"
  )
  
  # Specify the continuous-time dynamic model.
  # The main latent process (Xlevel) evolves depending on:
  # - its own previous value (beta)
  # - the secondary state Xa
  # - external shocks (delta_L)
  dynamics <- prep.formulaDynamics(
    formula = list(
      Xlevel ~ beta * Xlevel + Xa + delta_L,
      Xa ~ 0
    ),
    startval = c(beta = -0.2),
    isContinuousTime = TRUE
  )
  
  # Combine all model components into a dynr model object.
  model <- dynr.model(dynamics, meas, noise, initial, datadynr)
  
  # Estimate the model parameters using Kalman filter-based estimation.
  fit <- dynr.cook(model, verbose = FALSE)
  
  # Return both the fitted model and the model specification.
  list(fit = fit, model = model)
}

construct_shock_covariate <- function(data,
                                      shocks_locations,
                                      inn,
                                      level = "Xlevel",
                                      id_col = "id",
                                      time_col = "time") {
  
  # Create a working copy to avoid modifying the original dataset.
  data_out <- data
  
  # Initialize the shock covariate as zero (no shocks by default).
  data_out$delta_L <- 0
  
  # Match shock locations to dataset rows
  key_data   <- paste(data_out[[id_col]], data_out[[time_col]])
  key_shocks <- paste(shocks_locations[[id_col]], shocks_locations[[time_col]])
  
  idx <- match(key_shocks, key_data)
  
  valid <- !is.na(idx)
  
  if (!any(valid)) {
    warning("No matches were found between shocks_locations and data.")
    return(data_out)
  }
  
  # Validate latent level in innovation output
  delta_mat <- inn$delta.inn
  
  if (is.null(dimnames(delta_mat)) || !level %in% rownames(delta_mat)) {
    stop("The specified level was not found in inn$delta.inn.")
  }
  
  # Construct shock covariate from estimated innovations
  data_out$delta_L[idx[valid]] <- delta_mat[level, idx[valid]]
  
  # Return dataset with constructed covariate
  data_out
}

plot_subject <- function(data,
                         subject,
                         fit,
                         shocks_locations,
                         mode = c("latent", "both", "observed"),
                         ci = TRUE,
                         show_shocks = TRUE,
                         conf_level = 0.95,
                         latent_color = "#1D70B7",
                         latent_alpha = 0.15,
                         observed_color = "black",
                         shock_color = "red",
                         linetype = "dashed",
                         y_lims = NULL,
                         x_lims = NULL){
  
  # Validate and resolve the mode argument
  mode <- match.arg(mode)
  
  # Subset and sort the data for the selected subject
  df <- data[data$id == subject, ]
  df <- df[order(df$time), ]
  
  use_latent <- mode %in% c("latent", "both")
  use_observed <- mode %in% c("observed", "both")
  
  if(!use_latent && !use_observed){
    stop("At least one of 'latent' or 'observed' must be selected.")
  }
  
  p <- ggplot() +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # OBSERVED
  # Line + points share the same `color` mapping so ggplot merges them
  # into a single legend entry instead of creating a duplicate.
  if(use_observed){
    p <- p +
      geom_line(
        data = df,
        aes(x = time, y = X, color = "Observed"),
        linewidth = 0.7
      ) +
      geom_point(
        data = df,
        aes(x = time, y = X, color = "Observed"),
        size = 1.6
      )
  }
  
  # LATENT
  had_ci <- FALSE  # tracks whether a CI ribbon was actually drawn
  
  if(use_latent){
    
    # Locate this subject's rows in the original data to index into
    # the model's latent estimates (which are stored for all subjects)
    idx <- which(data$id == subject)
    
    eta <- fit$eta_smooth_final[1, idx]
    eta_var <- fit$error_cov_smooth_final[1, 1, idx]
    
    if(length(eta) != nrow(df)){
      stop("Mismatch between latent estimates and subject data.")
    }
    
    # Critical value for the requested confidence level
    z <- qnorm((1 + conf_level) / 2)
    
    df$eta <- eta
    
    if(ci){
      df$low <- eta - z * sqrt(eta_var)
      df$up  <- eta + z * sqrt(eta_var)
      
      p <- p +
        geom_ribbon(
          data = df,
          aes(x = time, ymin = low, ymax = up, fill = "Latent CI"),
          alpha = latent_alpha
        )
      
      had_ci <- TRUE
    }
    
    p <- p +
      geom_line(
        data = df,
        aes(x = time, y = eta, color = "Latent"),
        linewidth = 0.8
      )
  }
  
  # SHOCKS
  had_shocks <- FALSE  # tracks whether this subject actually has shocks
  
  if(show_shocks){
    
    s <- shocks_locations[shocks_locations$id == subject, ]
    
    if(nrow(s) > 0){
      p <- p +
        geom_vline(
          data = s,
          aes(xintercept = time, color = "Shocks"),
          linetype = linetype,
          linewidth = 0.6
        )
      had_shocks <- TRUE
    }
  }
  
  # AXIS LIMITS
  lims <- list()
  if(!is.null(x_lims)) lims$x <- x_lims
  if(!is.null(y_lims)) lims$y <- y_lims
  if(length(lims) > 0){
    p <- p + do.call(coord_cartesian, lims)
  }
  
  # DYNAMIC LEGEND
  
  # Which color levels actually exist in this plot, in the desired order.
  # Built dynamically so override.aes always matches the number of
  # legend entries actually present (avoids length-mismatch errors).
  present <- c("Latent", "Observed", "Shocks")[c(use_latent, use_observed, had_shocks)]
  
  # Full lookup tables for legend key appearance per series
  lty_lookup   <- c(Latent = "solid", Observed = "solid", Shocks = linetype)
  shape_lookup <- c(Latent = NA,      Observed = 16,       Shocks = NA)
  lwd_lookup   <- c(Latent = 0.8,     Observed = 0.7,       Shocks = 0.6)
  
  p <- p +
    scale_color_manual(
      breaks = present,
      values = c(
        "Observed" = observed_color,
        "Latent" = latent_color,
        "Shocks" = shock_color
      )
    ) +
    guides(
      color = guide_legend(
        order = 1,
        # Force each legend key to display the correct combination of
        # line type, point shape, and line width for its series
        override.aes = list(
          linetype = unname(lty_lookup[present]),
          shape = unname(shape_lookup[present]),
          linewidth = unname(lwd_lookup[present])
        )
      )
    ) +
    labs(x = "Time", y = "Level", color = NULL)
  
  # Only add the fill scale/legend if a CI ribbon was actually drawn;
  # otherwise this would trigger a "no shared levels" warning
  if(had_ci){
    p <- p +
      scale_fill_manual(values = c("Latent CI" = latent_color)) +
      guides(fill = guide_legend(order = 2)) +
      labs(fill = NULL)
  }
  
  p
}