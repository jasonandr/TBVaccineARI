# TB Vaccine Efficacy Bias Model

This repository contains the R code and data for analyzing bias in Tuberculosis Vaccine Efficacy (VE) estimates in high transmission settings.

## Project Structure

*   **`R/`**
    *   `model.R`: Core stochastic model functions (`VaccineModel`). Unified to support Base, Fast Progression, and Variable POD variants.
    *   `utils.R`: Helper functions for sampling and VE calculation.
*   **`data/`**
    *   `epi_parameters.csv`: Epidemiological parameters (ARI, LTBI).
    *   `calibration_base.csv` / `calibration_fastprog.csv`: Calibrated progression ratios.
    *   `power_data.rds`: Generated data for power analysis.
*   **`analysis/`**
    *   `01_calibrate.R`: Script to calibrate progression ratios.
    *   `02_results_base.R`: Runs the **Base Model Analysis** (VE plots for all cohorts).
    *   `03_results_fastprog.R`: Runs the **Fast Progression Model Analysis**.
    *   `04_results_variablePOD.R`: Runs the **Variable POD Model Analysis**.
    *   `05_power_generation.R`: Runs simulations to generate `power_data.rds`.
    *   `06_power_results.R`: Generates power curves and sample size calculations.
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

## Model Description

The model is a stochastic SEIR-type compartmental model simulated in discrete time steps. It accounts for vaccine heterogeneity (Leaky vs AoN), reinfection dynamics, and distinct disease pathways (Fast Progression).
