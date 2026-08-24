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
#' @param init Early/late split of prevalent infection: "stationary" (the
#'   specification used for all published results) or "ari" (superseded; retained
#'   for reference only -- see the note below).
#' @return A matrix of cumulative incidence (Imat), one column per stored state
#'   (t = 0, 0.1, ..., duration), for each simulation.
VaccineModel <- function(sims, popsize, epi.params, prog.ratio, QFTenrol, fastprogon,
                         vaxon, LeakyPOI, LeakyPOD, AoN.I, AoN.D,
                         AoN.D.pos = NULL, AoN.D.neg = NULL, init = "stationary") {
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
    # nsteps transitions across nsteps + 1 stored states, i.e. t = 0, 0.1, ..., duration.
    nsteps <- duration / timestep
    times <- seq(1, nsteps + 1, 1)

    lambda <- -log(1 - ARI)
    r1 <- prog.ratio * 0.05
    r2 <- prog.ratio * 0.0015
    e <- 1 / 0.89
    # omega is the RELATIVE RISK of infection among the previously infected: the
    # reinfection hazard is lambda * omega, i.e. 21% of the hazard faced by a
    # susceptible person (79% partial protection from prior infection, ref 19).
    # It is a natural history parameter applied identically in both trial arms;
    # vaccine protection against (re)infection is separate, via LeakyPOI / AoN.I.
    omega <- 0.21

    # Split of prevalent infection between early and late states.
    #  "stationary" - base specification: the split implied by a stationary open
    #                 population sustaining this force of infection and infection
    #                 prevalence (see stationary_EL_split()).
    #  "ari"        - superseded. Takes the fraction of infected individuals
    #                 infected within the past year to be ARI / prevalence. This
    #                 counts new infections against the whole population rather
    #                 than against susceptibles, so it overstates the early
    #                 fraction where prevalence is high (71% vs 26% at 50% ARI
    #                 and 70% prevalence). It is not the split implied by the
    #                 model's own natural history and is not used for any
    #                 published result; kept only so the calibration fallback
    #                 below stays defined and for reference.
    recent_frac <- ARI / trueltbiprev
    if (init == "stationary") {
        sp <- stationary_EL_split(lambda, trueltbiprev, r1, r2, e, omega)
        # No stationary solution exists for extreme progression rates (disease
        # removes infected individuals faster than the force of infection can
        # sustain the target prevalence). Those values are outside the calibrated
        # range; fall back on the ari split so the objective stays defined while
        # the calibration search passes through them.
        if (!is.na(sp$recent)) recent_frac <- sp$recent
    }
    recent_frac <- min(max(recent_frac, 0), 1)
    recentlyinfected <- ltbiprev * recent_frac
    distantlyinfected <- ltbiprev * (1 - recent_frac)

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

    # Largest-remainder apportionment so the six compartments sum exactly to popsize.
    # Rounding each independently could leave the arm up to two people short or over.
    frac <- c(
        S  = (1 - vaxon * AoN.D.neg) * (1 - vaxon * AoN.I) * (1 - ltbiprev),
        E  = ((1 - vaxon * AoN.D.pos) * (1 - vaxon * AoN.I)) * recentlyinfected,
        L  = ((1 - vaxon * AoN.D.pos) * (1 - vaxon * AoN.I)) * distantlyinfected,
        P  = vaxon * (ltbiprev * vaxon * AoN.D.pos + (1 - ltbiprev) * (1 - (1 - AoN.I) * (1 - AoN.D.neg))),
        EP = vaxon * (1 - AoN.D.pos) * AoN.I * recentlyinfected,
        LP = vaxon * (1 - AoN.D.pos) * AoN.I * distantlyinfected
    )
    frac[frac < 0] <- 0
    raw <- frac * popsize
    cnt <- floor(raw)
    short <- popsize - sum(cnt)
    if (short > 0) {
        ord <- order(raw - cnt, decreasing = TRUE)
        cnt[ord[seq_len(short)]] <- cnt[ord[seq_len(short)]] + 1
    }
    Smat[, 1] <- cnt[["S"]]
    Emat[, 1] <- cnt[["E"]]
    Lmat[, 1] <- cnt[["L"]]
    Imat[, 1] <- 0
    Pmat[, 1] <- cnt[["P"]]
    EPmat[, 1] <- cnt[["EP"]]
    LPmat[, 1] <- cnt[["LP"]]

    ### Define the rates
    R.infectionrate_S <- lambda * (1 - vaxon * LeakyPOI)
    R.primaryprogression <- r1 * (1 - LeakyPOD * vaxon)
    # Rates are annual; sample_competing_events() applies the timestep.
    R.EtoL <- e
    R.fastprogression <- r1 * lambda * fastprogon * (1 - LeakyPOD * vaxon) * (1 - LeakyPOI * vaxon)
    R.reactivation <- r2 * (1 - LeakyPOD * vaxon)
    R.reinfection <- lambda * omega * (1 - vaxon * LeakyPOI)

    for (i in 2:(nsteps + 1)) {
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

        # EP events: vaccinated early-infected either progress to disease or move to
        # late infection, at the same rates as their unvaccinated counterparts in E.
        EP_events <- sample_competing_events(EPmat[, i - 1], c(R.primaryprogression, R.EtoL), timestep)
        primaryprogressionEP <- EP_events[, 1]
        EPtoLP <- EP_events[, 2]

        # LP events
        reactivationEP <- ifelse(LPmat[, i - 1] == 0, 0, rbinom(sims, LPmat[, i - 1], prob = 1 - exp(-R.reactivation * timestep)))

        # state updates
        Smat[, i] <- Smat[, i - 1] - infectionrate_S
        Emat[, i] <- Emat[, i - 1] + infectionrate_S - primaryprogression - EtoL + reinfection - fastprogression
        Lmat[, i] <- Lmat[, i - 1] + EtoL - reactivation - reinfection
        Imat[, i] <- Imat[, i - 1] + primaryprogression + reactivation + fastprogression + primaryprogressionEP + reactivationEP
        EPmat[, i] <- EPmat[, i - 1] - primaryprogressionEP - EPtoLP
        LPmat[, i] <- LPmat[, i - 1] + EPtoLP - reactivationEP
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
                        AoN.D.pos = NULL, AoN.D.neg = NULL, init = "stationary") {
    # Lookup parameters for the specific setting
    curr_epi_params <- epi_data[epi_data$setting == epi_setting, ]
    curr_progratio <- prog_data$progratio[prog_data$setting == epi_setting]

    if (nrow(curr_epi_params) == 0) stop(paste("Epi setting not found:", epi_setting))

    res.novax <- VaccineModel(
        sims = sims, popsize = popsize, epi.params = curr_epi_params, prog.ratio = curr_progratio,
        QFTenrol = QFTenrol, fastprogon = fastprogon,
        vaxon = 0, LeakyPOI = LeakyPOI, LeakyPOD = LeakyPOD, AoN.I = AoN.I,
        AoN.D = AoN.D, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg, init = init
    )

    res.vax <- VaccineModel(
        sims = sims, popsize = popsize, epi.params = curr_epi_params, prog.ratio = curr_progratio,
        QFTenrol = QFTenrol, fastprogon = fastprogon,
        vaxon = 1, LeakyPOI = LeakyPOI, LeakyPOD = LeakyPOD, AoN.I = AoN.I,
        AoN.D = AoN.D, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg, init = init
    )

    timestep <- 0.1
    out <- VE(res.vax, res.novax, popsize, timestep = timestep)
    return(out)
}
