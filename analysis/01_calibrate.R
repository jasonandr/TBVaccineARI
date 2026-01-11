# Calibration Script
# Corresponds to: ModelCalibration_2025_10_08.R and ModelCalibration_FastProg...

library(readr)
library(tidyverse)
# library(stats) # optim is in stats, loaded by default

source("R/model.R")
source("R/utils.R")

# Load Target Data
# Note: Since the original EpiData CSV is missing, we use our recreated parameters to DRIVE the model,
# but calibration usually targets *Incidence* to find the *progression ratio*.
# We need target incidence.
# In original code: target.incidence.low, etc. came from EpiData.
# I will define them approximately to match the calibrated values if I were to run this.
# For now, this script documents HOW to calibration, but keeps the hardcoded "optimal" parameters as defaults unless run.

epi_data <- read_csv("data/epi_parameters.csv")

# Helper wrapper for optimization
IncidenceOut <- function(sims, popsize, epi_setting, progratio, fastprogon) {
    # Create temporary prog_data row
    temp_prog <- data.frame(setting = epi_setting, progratio = progratio)

    # Run model (no vaccine)
    res <- VaccineModel(
        sims = sims, popsize = popsize,
        epi.params = epi_data[epi_data$setting == epi_setting, ],
        prog.ratio = progratio,
        QFTenrol = "all", fastprogon = fastprogon,
        vaxon = 0, LeakyPOI = 0, LeakyPOD = 0, AoN.I = 0, AoN.D = 0
    )

    # Calculate Incidence (per 100k PY)
    # 1e5 * mean(res_final / PY)
    # PY approx = popsize * duration?
    # Refined PY: popsize * rowSums(timestep * (1 - res/popsize))
    timestep <- 0.1
    py <- popsize * rowSums(timestep * (1 - res / popsize))
    cases <- res[, ncol(res)]
    return(1e5 * mean(cases / py))
}

# Example Usage (Commented out to prevent long runtime during verification)
# calibrate_setting <- function(setting, target_inc, fastprogon) {
#   optim(par=2.5, fn=function(p) (IncidenceOut(1000, 5000, setting, p, fastprogon) - target_inc)^2, method="L-BFGS-B", lower=0.1, upper=10)$par
# }
