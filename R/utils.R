#' Helper Functions for TB VE Model
#'
#' This file contains utility functions for the TB Vaccine Efficacy model,
#' including multinomial sampling helpers and vaccine efficacy calculation logic.

library(rateratio.test)

#' Multinomial sampling from a matrix of probabilities
#'
#' @param mat A matrix where each column is a vector of simulations (one row per sim, one column per event type).
#' @return A matrix of event counts.
multinom.mat <- function(mat) {
    # mat needs to be formatted: cbind(v1, v2, v3)
    res <- apply(mat, 1, function(x) {
        total <- sum(x)
        if (total == 0) {
            return(rep(0, length(x))) # return all zeros
        } else {
            probs <- x / total
            return(as.vector(rmultinom(1, size = total, prob = probs)))
        }
    })
    return(t(res)) # returns a matrix where each column is vector of events, one row per sim
}

#' Sample competing events
#'
#' Multinomial model that takes a vector of states, a vector of rates, and a time-step
#' and returns the event counts.
#'
#' @param n_vec Vector of population sizes for each simulation.
#' @param rate_vec Vector of rates for the competing events.
#' @param dt Time step size.
#' @return A matrix of sampled events (rows=sims, cols=event types).
sample_competing_events <- function(n_vec, rate_vec, dt) {
    n_sims <- length(n_vec)

    # Replicate rate_vec across all simulations
    rate_mat <- matrix(rate_vec, nrow = n_sims, ncol = length(rate_vec), byrow = TRUE)

    # Convert rates to probabilities using discrete-time formula
    prob_mat <- 1 - exp(-rate_mat * dt)
    prob_mat[is.na(prob_mat) | prob_mat < 0] <- 0

    # Normalize if row sum > 1
    row_totals <- rowSums(prob_mat)
    over_1 <- row_totals > 1
    if (any(over_1)) {
        prob_mat[over_1, ] <- prob_mat[over_1, ] / row_totals[over_1]
    }

    # Add "no event" probability
    prob_remainder <- pmax(0, 1 - rowSums(prob_mat))
    prob_mat_full <- cbind(prob_mat, prob_remainder)

    # Preallocate result
    n_events <- ncol(prob_mat_full)
    draws <- matrix(0, nrow = n_sims, ncol = n_events)

    # Multinomial sampling
    for (i in seq_len(n_sims)) {
        draws[i, ] <- as.vector(rmultinom(1, size = n_vec[i], prob = prob_mat_full[i, ]))
    }

    return(draws) # Each row = one simulation; columns = event counts (+ no event)
}

#' Competing probability helper
#' @param p1 Probability 1
#' @param p2 Probability 2
p_compete <- function(p1, p2) {
    (
        return(ifelse(p1 == 0 & p2 == 0, 0, (p1 / (p1 + p2)) * (1 - (1 - p1) * (1 - p2))))
    )
}

#' Calculate Vaccine Efficacy
#'
#' @param res1 Matrix of Incidence/Outcomes for Vaccine arm.
#' @param res2 Matrix of Incidence/Outcomes for Control arm.
#' @param popsize Population size per arm.
#' @param times Time sequence (used for person-year calculation).
#' @param timestep Time step size.
#' @return A list containing VE estimate and p-values.
VE <- function(res1, res2, popsize, times, timestep = 0.1) {
    p.values <- rep(NA, length(res1[, 1]))

    # Cumulative Person-Years calculation
    # Note: logic adapted from original code: popsize * rowSums(timestep * (1 - res/popsize))
    # This approximates PY assuming prevalence is low or using prevalence-adjusted remaining time?
    # Actually, res contains cumulative incidence Imat.
    # The original logic: rowSums(timestep * (1 - res/popsize)) roughly sums up healthy time steps if res is prevalence?
    # In VaccineModel, res is Imat (cumulative incidence).
    # So (1 - I/N) is proportion healthy. Summing this * timestep gives PY at risk.

    PYlist1 <- popsize * rowSums(timestep * (1 - res1 / popsize))
    PYlist2 <- popsize * rowSums(timestep * (1 - res2 / popsize))

    Ilist1 <- res1[, max(times)]
    Ilist2 <- res2[, max(times)]

    for (i in 1:length(res1[, 1])) {
        I1 <- Ilist1[i]
        I2 <- Ilist2[i]
        PY1 <- PYlist1[i]
        PY2 <- PYlist2[i]
        # Handle potential zeros or NAs if needed, but rateratio.test handles basic cases
        p.values[i] <- as.numeric(rateratio.test(c(I1, I2), c(PY1, PY2))[1])
    }

    return(list(
        vax.eff = 1 - (Ilist1 / PYlist1) / (Ilist2 / PYlist2),
        p.values = p.values
    ))
}
