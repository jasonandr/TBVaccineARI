# TB Vaccine Efficacy Bias Model

This repository contains the R code and data for the analysis of **bias in Tuberculosis Vaccine Efficacy (VE) estimates in high transmission settings**.

## Project Structure

*   `R/`
    *   `model.R`: Core stochastic model functions (`VaccineModel`). Unified to support Base, Fast Progression, and Variable POD variants.
    *   `utils.R`: Helper functions for sampling and VE calculation.
*   `data/`
    *   `epi_parameters.csv`: Epidemiological parameters (ARI, LTBI).
    *   `calibration_base.csv` / `calibration_fastprog.csv`: Calibrated progression ratios.
    *   `power_data.rds`: Generated data for power analysis (output of `05_power_generation.R`).
*   `analysis/`
    *   `01_calibrate.R`: Calibration script.
    *   `02_results_base.R`: **Base Model Analysis**. Generates VE plots for QFT+, QFT-, and All cohorts.
    *   `03_results_fastprog.R`: **Fast Progression Model Analysis**. Generates VE plots for all cohorts.
    *   `04_results_variablePOD.R`: **Variable POD Model Analysis**.
    *   `05_power_generation.R`: **Power Analysis (Simulation)**. Generates `power_data.rds` (computationally intensive).
    *   `06_power_results.R`: **Power Analysis (Results)**. Generates power curves and sample size/power text summaries.
*   `figures/`: Generated plots (10+ PDFs).

## Reproducibility

### 1. Setup
Open the project in RStudio (or set working directory to `TB_VE_Model/`).
Install packages:
```r
install.packages(c("tidyverse", "rateratio.test", "ggsci", "ggpubr", "mgcv"))
```

### 2. Run Analysis
You can run the scripts individually:

**Base Model Results:**
```r
source("analysis/02_results_base.R")
# Outputs: figures/02_Base_VE_onlypos.pdf, _all.pdf, _onlyneg.pdf + Text Summaries
```

**FastProg Model Results:**
```r
source("analysis/03_results_fastprog.R")
# Outputs: figures/03_FastProg_VE_onlypos.pdf, ...
```

**Power Analysis:**
First, generate the data (warning: takes time):
```r
source("analysis/05_power_generation.R")
```
Then, visualize and calculate sample sizes:
```r
source("analysis/06_power_results.R")
# Outputs: figures/05_Power_Curves.pdf + Text Results for Sample Size
```

## Model Description

The model is a stochastic SEIR-type compartmental model simulated in discrete time steps. It accounts for vaccine heterogeneity (Leaky vs AoN), reinfection dynamics, and distinct disease pathways (Fast Progression).
