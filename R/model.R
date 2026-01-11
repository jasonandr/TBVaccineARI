#' TB Vaccine Efficacy Model (Stochastic) - Unified
#'
#' @param sims Number of simulations to run.
#' @param popsize Population size per simulation (cohort size).
#' @param epi.params One-row data.frame or list containing `ari` and `ltbi`.
#' @param prog.ratio Progression scalar (from calibration).
#' @param QFTenrol Enrollment strategy: "all", "onlypos", "onlyneg".
#' @param fastprogon Binary toggle (0/1) for fast progression path.
#' @param vaxon Binary toggle (0/1) for vaccination.
#' @param LeakyPOI Vaccine efficacy against Infection (0-1).
#' @param LeakyPOD Vaccine efficacy against Disease (0-1).
#' @param AoN.I All-or-Nothing protection against Infection (proportion protected).
#' @param AoN.D All-or-Nothing protection against Disease (proportion protected).
#' @param AoN.D.pos All-or-Nothing protection against Disease for QFT+ (defaults to AoN.D).
#' @param AoN.D.neg All-or-Nothing protection against Disease for QFT- (defaults to AoN.D).
#' @return A matrix of cumulative incidence (Imat) for each simulation over time.
VaccineModel <- function(sims, popsize, epi.params, prog.ratio, QFTenrol, fastprogon,
                         vaxon, LeakyPOI, LeakyPOD, AoN.I, AoN.D,
                         AoN.D.pos = NULL, AoN.D.neg = NULL) {
    # Handle defaults for Variable POD
    if (is.null(AoN.D.pos)) AoN.D.pos <- AoN.D
    if (is.null(AoN.D.neg)) AoN.D.neg <- AoN.D

    # Extract parameters
    ARI <- epi.params$ari
    ltbiprev <- epi.params$ltbi

    # Adjust LTBI prevalence based on enrollment criteria
    trueltbiprev <- ltbiprev
    if (QFTenrol == "onlypos") {
        ltbiprev <- 1
    }
    if (QFTenrol == "onlyneg") {
        ltbiprev <- 0
    }

    # Fixed Parameters
    timestep <- 0.1
    duration <- 3
    times <- seq(1, duration / timestep, 1)

    lambda <- -log(1 - ARI)
    recentlyinfected <- ltbiprev * (ARI / trueltbiprev)
    distantlyinfected <- ltbiprev * (1 - ARI / trueltbiprev)
    r1 <- prog.ratio * 0.05
    r2 <- prog.ratio * 0.0015
    e <- 1 / 0.89
    omega <- 0.21

    ### Create the baseline states
    Smat <- matrix(NA, nrow = sims, ncol = length(times))
    Emat <- Smat
    Lmat <- Smat
    Imat <- Smat
    Pmat <- Smat
    EPmat <- Smat
    LPmat <- Smat

    # Initial conditions:
    # Use AoN.D.neg for S (QFT-negative / Uninfected at baseline)
    # Use AoN.D.pos for E/L (QFT-positive / Infected at baseline)
    # Smat = ... (1-vaxon*AoN.D.neg) ...
    # Emat = ... (1-vaxon*AoN.D.pos) ...
    # This implies AoN.D.neg applies to naive recruits, AoN.D.pos applies to prevalent infected recruits.

    Smat[, 1] <- round((1 - vaxon * AoN.D.neg) * (1 - vaxon * AoN.I) * (1 - ltbiprev) * popsize, 0)
    Emat[, 1] <- round(popsize * ((1 - vaxon * AoN.D.pos) * (1 - vaxon * AoN.I)) * recentlyinfected, 0)
    Lmat[, 1] <- round(popsize * ((1 - vaxon * AoN.D.pos) * (1 - vaxon * AoN.I)) * distantlyinfected, 0)
    Imat[, 1] <- 0
    Pmat[, 1] <- round(vaxon * (ltbiprev * vaxon * AoN.D.pos + (1 - ltbiprev) * (1 - (1 - AoN.I) * (1 - AoN.D.neg))) * popsize, 0)
    EPmat[, 1] <- round(vaxon * popsize * (1 - AoN.D.pos) * AoN.I * recentlyinfected, 0)
    LPmat[, 1] <- round(vaxon * popsize * (1 - AoN.D.pos) * AoN.I * distantlyinfected, 0)

    ### Define the rates
    R.infectionrate_S <- lambda * (1 - vaxon * LeakyPOI)
    R.primaryprogression <- r1 * (1 - LeakyPOD * vaxon)
    R.EtoL <- e * timestep
    R.fastprogression <- r1 * lambda * fastprogon * (1 - LeakyPOD * vaxon) * (1 - LeakyPOI * vaxon)
    R.reactivation <- r2 * (1 - LeakyPOD * vaxon)
    R.reinfection <- lambda * omega * (1 - vaxon * LeakyPOI)

    for (i in 2:max(times)) {
        # S events
        infectionrate_S <- sample_competing_events(Smat[, i - 1], R.infectionrate_S, timestep)[, 1]

        # E events
        E_events <- sample_competing_events(Emat[, i - 1], c(R.primaryprogression, R.EtoL, R.fastprogression), timestep)
        primaryprogression <- E_events[, 1]
        EtoL <- E_events[, 2]
        fastprogression <- E_events[, 3]

        # L events
        L_events <- sample_competing_events(Lmat[, i - 1], c(R.reactivation, R.reinfection), timestep)
        reactivation <- L_events[, 1]
        reinfection <- L_events[, 2]

        # EP events
        primaryprogressionEP <- ifelse(EPmat[, i - 1] == 0, 0, rbinom(sims, EPmat[, i - 1], prob = 1 - exp(-R.primaryprogression * timestep)))

        # LP events
        reactivationEP <- ifelse(LPmat[, i - 1] == 0, 0, rbinom(sims, LPmat[, i - 1], prob = 1 - exp(-R.reactivation * timestep)))

        # state updates
        Smat[, i] <- Smat[, i - 1] - infectionrate_S
        Emat[, i] <- Emat[, i - 1] + infectionrate_S - primaryprogression - EtoL + reinfection - fastprogression
        Lmat[, i] <- Lmat[, i - 1] + EtoL - reactivation - reinfection
        Imat[, i] <- Imat[, i - 1] + primaryprogression + reactivation + fastprogression + primaryprogressionEP + reactivationEP
        EPmat[, i] <- EPmat[, i - 1] - primaryprogressionEP
        LPmat[, i] <- LPmat[, i - 1] - reactivationEP
    }
    return(Imat)
}

#' Wrapper to run Paired Trials and calculate VE
#'
#' @param sims Number of simulations.
#' @param popsize Population per arm.
#' @param epi_setting "low", "medium", or "high".
#' @param epi_data Data.frame of ARI/LTBI parameters.
#' @param prog_data Data.frame of calibration parameters.
#' @param QFTenrol "all", "onlypos", "onlyneg".
#' @param fastprogon Binary.
#' @param LeakyPOI ...
#' @param LeakyPOD ...
#' @param AoN.I ...
#' @param AoN.D ...
#' @param AoN.D.pos ...
#' @param AoN.D.neg ...
#' @return List with `vax.eff` vector and `p.values` vector.
VEfromModel <- function(sims, popsize, epi_setting, epi_data, prog_data,
                        QFTenrol, fastprogon,
                        LeakyPOI, LeakyPOD, AoN.I, AoN.D = NULL,
                        AoN.D.pos = NULL, AoN.D.neg = NULL) {
    # Lookup parameters for the specific setting
    curr_epi_params <- epi_data[epi_data$setting == epi_setting, ]
    curr_progratio <- prog_data$progratio[prog_data$setting == epi_setting]

    if (nrow(curr_epi_params) == 0) stop(paste("Epi setting not found:", epi_setting))

    res.novax <- VaccineModel(
        sims = sims, popsize = popsize, epi.params = curr_epi_params, prog.ratio = curr_progratio,
        QFTenrol = QFTenrol, fastprogon = fastprogon,
        vaxon = 0, LeakyPOI = LeakyPOI, LeakyPOD = LeakyPOD, AoN.I = AoN.I,
        AoN.D = AoN.D, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg
    )

    res.vax <- VaccineModel(
        sims = sims, popsize = popsize, epi.params = curr_epi_params, prog.ratio = curr_progratio,
        QFTenrol = QFTenrol, fastprogon = fastprogon,
        vaxon = 1, LeakyPOI = LeakyPOI, LeakyPOD = LeakyPOD, AoN.I = AoN.I,
        AoN.D = AoN.D, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg
    )

    duration <- 3
    timestep <- 0.1
    times <- seq(1, duration / timestep, 1)

    out <- VE(res.vax, res.novax, popsize, times, timestep)
    return(out)
}
