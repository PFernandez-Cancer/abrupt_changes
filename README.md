# Detecting Abrupt Changes in Developmental Trajectories

This repository contains a tutorial and R code for detecting **abrupt changes (innovative outliers, or "shocks")** in individual developmental trajectories using continuous-time state-space models fitted with the [`dynr`](https://cran.r-project.org/package=dynr) package.

The method fits a **deterministic** continuous-time dual latent change score model (no process noise), uses Kalman smoother diagnostics to flag time points where the observed trajectory deviates beyond what the model's dynamics predict, and then refits a **stochastic** model that explicitly incorporates the detected shocks as an exogenous covariate — separating "normal" process noise from discrete, abrupt jumps.

## 📄 Based on

This code accompanies:

> Romero-Suárez, M., Chow, S.-M., Cáncer, P. F., & Estrada, E. (2026). *Detecting Abrupt Changes in Developmental Trajectories with Continuous-Time State-Space Models* [manuscript in preparation].

See the paper for the full simulation study, the motivation behind the combined detection criterion (ρ² & θ), and the rationale for the θ = 1.65 threshold used by default.

## 🔍 View the tutorial

A rendered, step-by-step walkthrough is available here:

**[https://PFernandez-Cancer.github.io/abrupt_changes/AbruptChanges.html](https://tu-usuario.github.io/tu-repo/AbruptChanges.html)**

No installation needed — just open the link.

## 📁 Repository contents

| File | Description |
|---|---|
| `AbruptChanges.Rmd` | Annotated tutorial (R Markdown source) explaining each step of the detection pipeline |
| `AbruptChanges.html` | Rendered version of the tutorial (viewable via the link above) |
| `functions_AbruptChanges.R` | Helper functions used throughout the pipeline: `fit_deterministic()`, `detect_shocks()`, `construct_shock_covariate()`, `fit_stochastic_with_shocks()`, `plot_subject()` |
| `code_AbruptChanges.R` | Standalone script running the full pipeline end to end (equivalent to the Rmd, without the explanatory text) |
| `data.csv` | Example panel dataset (`id`, `time`, `X`) used throughout the tutorial |

## 🚀 Running it yourself

1. Clone this repository:
```bash
   git clone https://github.com/tu-usuario/tu-repo.git
```
2. Open the project in RStudio.
3. Make sure the required packages are installed:
```r
   install.packages(c("dynr", "ggplot2"))
```
4. Run `code_AbruptChanges.R` directly, or knit `AbruptChanges.Rmd` for the full walkthrough with explanations.

## 🧩 Overall logic

1. Fit a **deterministic** model (no process noise) so that any unexplained variability is "exposed" as a candidate shock.
2. Use `dynr.taste()` to compute innovation diagnostics and detect shocks in specific subjects and time points.
3. Build a covariate (`delta_L`) capturing the magnitude of each detected shock.
4. Refit a **stochastic** model that incorporates `delta_L`, so the model distinguishes ordinary process noise from explicit shocks.
5. Inspect which subjects were flagged and at which time points.
6. Visualize individual latent trajectories, with detected shocks marked.

## 📬 Contact

Questions about the method or the code: [marcos.romero@uam.es]
