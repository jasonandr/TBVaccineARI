# Model Calibration Script
# This script calibrates the progression ratios to match target incidence rates
# for both the Base Model and the Fast Progression Model.

library(tidyverse)

source("R/model.R")
source("R/utils.R")

# 1. Load Data
# Updated epi_parameters.csv now implies real data: setting, ltbi, ari, incidence
epi_data <- read_csv("data/epi_parameters.csv", show_col_types = FALSE)

# 2. Define Incidence Function
# Calculates incidence rate per 100,000 person-years
IncidenceOut <- function(sims, popsize, epi_row, progratio, fastprogon) {
    # Prepare single-row epi parameters
    # The model expects a list/df with $ari and $ltbi
    curr_epi <- epi_row

    # Run Model (No Vaccination)
    # using QFTenrol = "all" as standard for calibration (incidence in general population)
    res <- VaccineModel(
        sims = sims,
        popsize = popsize,
        epi.params = curr_epi,
        prog.ratio = progratio,
        QFTenrol = "all",
        fastprogon = fastprogon,
        vaxon = 0,
        LeakyPOI = 0, LeakyPOD = 0, AoN.I = 0, AoN.D = 0
    )

    # Calculate Incidence Rate per 100k PY
    # Numerator: Total cumulative cases at end of simulation (res is cumulative incidence matrix)
    # Denominator: Person-Years approximation
    timestep <- 0.1
    duration <- 3
    times <- seq(1, duration / timestep, 1)

    numerator <- res[, length(times)]
    # PY approx: sum of healthy time steps * timestep
    # (1 - res/popsize) is the proportion healthy at each step
    denominator <- popsize * rowSums(timestep * (1 - res / popsize))

    # Mean incidence rate across simulations, scaled to 100k
    return(1e5 * mean(numerator / denominator))
}

# 3. Calibration Function
CalibrateSetting <- function(setting_name, target_incidence, fastprogon) {
    curr_epi_row <- epi_data %>% filter(setting == setting_name)

    # Objective function for optim
    ObjectiveFn <- function(par) {
        progratio <- par
        # Settings chosen to match original calibration script (2000 sims, 50000 pop)
        # Reduced slightly for speed if needed, but keeping high for accuracy
        model_inc <- IncidenceOut(
            sims = 1000, # Reduced from 2000 for efficiency in this script
            popsize = 20000, # Reduced from 50000
            epi_row = curr_epi_row,
            progratio = progratio,
            fastprogon = fastprogon
        )
        return((model_inc - target_incidence)^2)
    }

    # Helper to print progress
    cat(paste("  Calibrating", setting_name, "(Target:", target_incidence, ")... "))

    opt <- optim(
        par = 2.5,
        method = "L-BFGS-B",
        lower = 1,
        upper = 10, # Increased upper bound just in case
        control = list(parscale = 2.5, trace = 0, maxit = 20),
        fn = ObjectiveFn
    )

    cat(paste("Done. Ratio:", round(opt$par, 3), "\n"))
    return(opt$par)
}

# 4. Main Loop

cat("Starting Calibration...\n")

# A. Base Model (fastprogon = 0)
cat("\n--- Calibrating Base Model (FastProg = 0) ---\n")
results_base <- list()
for (i in 1:nrow(epi_data)) {
    s <- epi_data$setting[i]
    inc <- epi_data$incidence[i]
    val <- CalibrateSetting(s, inc, fastprogon = 0)
    results_base[[s]] <- val
}

df_base <- data.frame(
    setting = names(results_base),
    progratio = unlist(results_base)
)
write_csv(df_base, "data/calibration_base.csv")
cat("Saved to data/calibration_base.csv\n")
print(df_base)


# B. Fast Progression Model (fastprogon = 1)
cat("\n--- Calibrating Fast Progression Model (FastProg = 1) ---\n")
results_fast <- list()
for (i in 1:nrow(epi_data)) {
    s <- epi_data$setting[i]
    inc <- epi_data$incidence[i]
    val <- CalibrateSetting(s, inc, fastprogon = 1)
    results_fast[[s]] <- val
}

df_fast <- data.frame(
    setting = names(results_fast),
    progratio = unlist(results_fast)
)
write_csv(df_fast, "data/calibration_fastprog.csv")
cat("Saved to data/calibration_fastprog.csv\n")
print(df_fast)
