# Model Calibration Script
# This script calibrates the progression ratios to match target incidence rates
# for both the Base Model and the Fast Progression Model.
#
# Simulated incidence is a stochastic function of the progression ratio, so a
# general-purpose optimiser working on a fresh set of random draws at each
# evaluation is unreliable (its numerical gradients are dominated by Monte Carlo
# noise). We instead fix the random seed inside the objective. That makes the
# objective deterministic given the seed -- it is a single stochastic realisation,
# not a smooth function, and monotonicity is empirical rather than guaranteed --
# and solve for the target by bisection. Seed sensitivity is checked below.

library(tidyverse)

source("R/model.R")
source("R/utils.R")

# 1. Load Data
# setting, ltbi, ari, incidence
epi_data <- read_csv("data/epi_parameters.csv", show_col_types = FALSE)

# Calibration precision. 200 simulations of 50,000 people gives 10 million
# person-simulations per evaluation, which is ample.
CAL_SIMS <- 200
CAL_POP <- 50000
CAL_SEED <- 1

# 2. Define Incidence Function
# Calculates incidence rate per 100,000 person-years
IncidenceOut <- function(sims, popsize, epi_row, progratio, fastprogon, seed = CAL_SEED,
                         init = "ari") {
    # Common random numbers: the same draws are used at every progression ratio,
    # so the objective is deterministic given the seed.
    set.seed(seed)

    # Run Model (No Vaccination)
    # using QFTenrol = "all" as standard for calibration (incidence in general population)
    res <- VaccineModel(
        sims = sims,
        popsize = popsize,
        epi.params = epi_row,
        prog.ratio = progratio,
        QFTenrol = "all",
        fastprogon = fastprogon,
        vaxon = 0,
        LeakyPOI = 0, LeakyPOD = 0, AoN.I = 0, AoN.D = 0,
        init = init
    )

    timestep <- 0.1
    nt <- ncol(res)

    # Numerator: cumulative cases at the end of follow-up
    numerator <- res[, nt]
    # Denominator: disease-free person-time, integrated over the nt - 1 intervals
    # between stored states so it spans exactly the same window as the numerator.
    left <- res[, 1:(nt - 1), drop = FALSE]
    right <- res[, 2:nt, drop = FALSE]
    denominator <- popsize * rowSums(timestep * (1 - (left + right) / (2 * popsize)))

    # Pooled population incidence rate (total cases / total person-time), which is
    # the estimand the target incidence refers to, rather than a mean of
    # replicate-specific rates.
    return(1e5 * sum(numerator) / sum(denominator))
}

# 3. Calibration Function
CalibrateSetting <- function(setting_name, target_incidence, fastprogon, init = "ari") {
    curr_epi_row <- epi_data %>% filter(setting == setting_name)

    ObjectiveFn <- function(progratio) {
        IncidenceOut(CAL_SIMS, CAL_POP, curr_epi_row, progratio, fastprogon, init = init) -
            target_incidence
    }

    cat(paste("  Calibrating", setting_name, "(Target:", target_incidence, ")... "))

    lower <- 0.05
    upper <- 20
    if (ObjectiveFn(lower) > 0) {
        warning(paste("Target incidence below achievable range for", setting_name))
        ratio <- lower
    } else if (ObjectiveFn(upper) < 0) {
        warning(paste("Target incidence above achievable range for", setting_name))
        ratio <- upper
    } else {
        ratio <- uniroot(ObjectiveFn, c(lower, upper), tol = 1e-4)$root
    }

    achieved <- IncidenceOut(CAL_SIMS, CAL_POP, curr_epi_row, ratio, fastprogon, init = init)
    cat(paste0("Done. Ratio: ", round(ratio, 3), " (incidence ", round(achieved), ")\n"))
    return(ratio)
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


# C. Seed sensitivity check
# The objective is conditioned on a single RNG stream. Confirm the root does not
# depend materially on which stream is used.
cat("\n--- Seed sensitivity of the calibrated ratio (base model) ---\n")
for (i in 1:nrow(epi_data)) {
    row <- epi_data[i, ]
    roots <- sapply(1:8, function(sd) {
        f <- function(pr) IncidenceOut(CAL_SIMS, CAL_POP, row, pr, 0, seed = sd) - row$incidence
        uniroot(f, c(0.05, 20), tol = 1e-4)$root
    })
    cat(sprintf("  %-7s mean %.3f  sd %.4f  range %.3f-%.3f\n",
                row$setting, mean(roots), sd(roots), min(roots), max(roots)))
}

# D. Calibration for the stationary-initialisation sensitivity analysis
cat("\n--- Calibrating Base Model, stationary E/L initialisation ---\n")
results_stat <- list()
for (i in 1:nrow(epi_data)) {
    s <- epi_data$setting[i]
    results_stat[[s]] <- CalibrateSetting(s, epi_data$incidence[i], fastprogon = 0,
                                          init = "stationary")
}
df_stat <- data.frame(setting = names(results_stat), progratio = unlist(results_stat))
write_csv(df_stat, "data/calibration_base_stationary.csv")
cat("Saved to data/calibration_base_stationary.csv\n")
print(df_stat)
