# TB Vaccine Efficacy Bias Model

This repository contains the R code and data for analyzing bias in Tuberculosis Vaccine Efficacy (VE) estimates in high transmission settings.

The early/late split of prevalent infection at enrollment uses the stationary open-population solution (`init = "stationary"`, the default in `VaccineModel`). This is the specification used for all published results.

## Project Structure

*   **`R/`**
    *   `model.R`: Core stochastic model functions (`VaccineModel`). Unified to support Base, Fast Progression, and Variable POD variants.
    *   `utils.R`: Helper functions for sampling and VE calculation.
*   **`data/`**
    *   `epi_parameters.csv`: Epidemiological parameters (ARI, LTBI).
    *   `calibration_base.csv` / `calibration_fastprog.csv`: Calibrated progression ratios for the base and fast-progression models.
    *   `calibration_ariprev.csv`: Calibrated ratios under the superseded ARI/prevalence early/late split; not used for any published result (see `model.R`, `init`).
    *   `power_data.rds`: Generated data for power analysis.
*   **`analysis/`**
    *   `01_calibrate.R`: Script to calibrate progression ratios.
    *   `02_results_base.R`: Runs the **Base Model Analysis** (VE plots for all cohorts).
    *   `03_results_fastprog.R`: Runs the **Fast Progression Model Analysis**.
    *   `04_results_variablePOD.R`: Runs the **Variable POD Model Analysis**.
    *   `05_power_generation.R`: Runs simulations to generate `power_data.rds`.
    *   `06_power_results.R`: Generates power curves and sample size calculations.
    *   `07_power_curves.R`: Power by sample size, archetype, cohort and setting -> `data/power_curves.csv`.
    *   `08_manuscript_figures.R`: Reproduces the manuscript result figures (2, 3, S2-S9) and Supplemental Table S1.
    *   `09_interactive_plot_data.R`: Figure 2 draws as a flat `.csv` for the journal's Interactive Plot Viewer.
    *   `schematics/make_schematics.py`: Draws the model schematics (Figures 1 and S1) as PDF/PNG/SVG.
    *   `schematics/audit.py`: Checks that every arrow drawn in the schematics is a transition the model implements.
*   **`figures/`**: Output directory for generated plots.

## Usage

### 1. Setup
Open the project in RStudio (or set working directory to `TB_VE_Model/`).
Install required packages:
```r
install.packages(c("tidyverse", "rateratio.test", "ggsci", "ggpubr", "mgcv"))
```

### 2. Run Analysis
You can run the analysis scripts individually:

**Base Model Results:**
```r
source("analysis/02_results_base.R")
# Outputs summary statistics and figures.
```

**FastProg Model Results:**
```r
source("analysis/03_results_fastprog.R")
```

**Power Analysis:**
First, generate the data (computationally intensive):
```r
source("analysis/05_power_generation.R")
```
Then, generate the power curves and sample size estimates:
```r
source("analysis/06_power_results.R")
```

### 3. Reproduce the manuscript figures

Requires `matplotlib` for the schematics; the rest is R. From the repository root:

```bash
Rscript analysis/07_power_curves.R          # ~25 min; writes data/power_curves.csv
Rscript analysis/08_manuscript_figures.R    # ~8 min;  Figures 2, 3, S2-S9 + Table S1
Rscript analysis/09_interactive_plot_data.R # Figure 2 draws as .csv
python analysis/schematics/make_schematics.py   # Figures 1 and S1
python analysis/schematics/audit.py             # verifies the schematics against the model
```

Output goes to `figures/manuscript/` by default (`FIG_OUT` to change it). Set
`FIG_FMT=pdf` for vector output instead of PNG:

```bash
FIG_FMT=pdf Rscript analysis/08_manuscript_figures.R
```

`07_power_curves.R` must run before `08_manuscript_figures.R`, which reads its output
for Figure 3. `data/power_curves.csv` is committed, so you can skip step 1.

Every draw is seeded, so Figures 2, S2-S9 and Supplemental Table S1 reproduce exactly.
Figure 3 is redrawn from `data/power_curves.csv`; the original sample-size grid was not
retained, so its curves differ from the published figure by Monte Carlo noise. The
sample sizes at 90% power (`data/power_sample_sizes.csv`) agree with those reported in
the paper to within a few participants per arm.

## Model Description

The model is a stochastic SEIR-type compartmental model simulated in discrete time steps. It accounts for different vaccine archetypes ("Leaky" vs "All-or-None"), reinfection dynamics, and distinct disease pathways (Fast Progression).
