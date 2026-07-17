#### Load required packages, data and external functions ####

library(dynr)
library(ggplot2)

source("functions_AbruptChanges.R")
data <- read.csv("data.csv")

names(data) <- c("id", "time", "X" )
data$time <- data$time - min(data$time)
data$X <- as.numeric(scale(data$X))

# Remove rows without info
data <- data[!is.na(data$X),]
data <- data[!is.na(data$time),]

# KEY SAFEGUARD: Dynr requires data sorted by case, then time
data <- data[order(data$id, data$time), ]

#### STEP 1. Fit deterministic continuous-time model ####

# The latent system follows a deterministic dual latent change score structure,
# combining state dependence (autoregressive component) and a latent change component.
#
# In this specification, process noise is excluded, so the evolution of the
# latent trajectory is fully determined by the structural dynamics of the model.
#
# As a consequence, any deviations from the expected latent trajectory are
# interpreted as candidate innovation shocks rather than stochastic process fluctuations.

res1 <- fit_deterministic(data)

#### STEP 2. Extract shock diagnostics and detect events ####

# The model computes innovation-based diagnostics to detect innovative outliers.
# Excludes only the first time point per subject, as it lacks backward support and gives an unreliable statistic.

shock_diagnostics <- dynr.taste(
  res1$model,
  res1$fit,
  which.state = "Xlevel"
)

shocks_locations <- detect_shocks(shock_diagnostics, data, method = "combined")

#### STEP 3. Construct shock covariate ####

# A time-varying covariate (delta_L) is created to represent
# the estimated magnitude of detected shocks.
#
# This covariate is initially set to zero and is updated only
# at time points where shocks have been identified.

data2 <- construct_shock_covariate(
  data = data,
  shocks_locations = shocks_locations,
  inn = shock_diagnostics,
  level = "Xlevel"
)

#### STEP 4. Fit stochastic model with shocks ####

# The detected shocks are incorporated into the dynamic model
# as an exogenous time-varying covariate.
#
# This allows the latent process to account for abrupt changes
# that are not explained by deterministic dynamics alone.

res2 <- fit_stochastic_with_shocks(data2)

summary(res2$fit)

#### STEP 5. Inspect detected shocks ####

# We examine:
#   - which individuals contain shocks
#   - the time points at which shocks occur
#   - the estimated magnitude of each shock

shock_ids <- unique(shocks_locations$id)
shock_ids

data2[data2$delta_L != 0,
      c("id", "time", "delta_L")]

#### STEP 6. Plot individual trajectory ####

# Visualize the estimated latent trajectory (and, optionally, the
# observed data) for individual subjects, marking shocks when present.
# We compare a subject known to have a shock against one that doesn't,
# each shown in two modes: latent-only and latent+observed.

# ---- Subject WITH a shock ----

# Latent estimate + CI band, shock marked as a vertical line
c_id <- 2L
plot_subject(
  data = data2,
  subject = shock_ids[c_id], # id of the subject
  fit = res2$fit,
  shocks_locations = shocks_locations,
  mode = "latent",
  ci = TRUE,
  show_shocks = TRUE,
  plot_group = TRUE
)

# Same subject, now overlaying observed data on top of the latent estimate
plot_subject(
  data = data2,
  subject = shock_ids[c_id], # same subject
  fit = res2$fit,
  shocks_locations = shocks_locations,
  mode = "both",
  ci = TRUE,
  show_shocks = TRUE,
  plot_group = TRUE
)

# ---- Subject WITHOUT a shock (control/baseline case) ----

# Latent estimate + CI band; show_shocks = TRUE has no visual effect here
# since this subject has no entries in shocks_locations
plot_subject(
  data = data2,
  subject = data2$id[10], # id of the subject
  fit = res2$fit,
  shocks_locations = shocks_locations,
  mode = "latent",
  ci = TRUE,
  show_shocks = TRUE, # no visual effect: this subject has no entries in shocks_locations
  plot_group = TRUE
)

# Same subject, latent + observed overlay for direct comparison
plot_subject(
  data = data2,
  subject = data2$id[10],
  fit = res2$fit,
  shocks_locations = shocks_locations,
  mode = "both",
  ci = TRUE,
  show_shocks = TRUE,
  plot_group = TRUE
)
