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
    r <- rate_vec
    r[is.na(r) | r < 0] <- 0
    k <- length(r)
    R <- sum(r)

    if (R <= 0) {
        return(matrix(0, nrow = n_sims, ncol = k + 1))
    }

    # Exact continuous-time competing risks: the probability that ANY event occurs
    # in dt is 1 - exp(-sum(r) * dt), split between causes in proportion to their
    # hazards. Converting each cause separately with 1 - exp(-r_j * dt) overstates
    # both the total exit probability and the share of the faster causes.
    p_any <- 1 - exp(-R * dt)
    prob_full <- c(p_any * r / R, 1 - p_any)

    draws <- matrix(0, nrow = n_sims, ncol = k + 1)
    for (i in seq_len(n_sims)) {
        if (n_vec[i] > 0) {
            draws[i, ] <- as.vector(rmultinom(1, size = n_vec[i], prob = prob_full))
        }
    }

    return(draws) # Each row = one simulation; columns = event counts (+ no event)
}

#' Stationary early/late split of prevalent infection
#'
#' The base specification takes the fraction of infected individuals who are in the
#' early state to be ARI / prevalence. That is not the split implied by the
#' model's own natural history: it counts new infections against the whole
#' population rather than against susceptibles, and it ignores the fact that the
#' early state drains with mean duration 1/e.
#'
#' This helper returns the split implied by a stationary open population (entry
#' and exit rate mu, chosen so that infection prevalence matches the setting):
#'   dS/dt = mu - (lambda + mu) S
#'   dE/dt = lambda S + omega lambda L - (r1 + e + mu) E
#'   dL/dt = e E - (r2 + omega lambda + mu) L
#'
#' @return list(recent = fraction of prevalent infection in the early state,
#'              mu = implied population turnover rate per year)
stationary_EL_split <- function(lambda, prevalence, r1, r2, e, omega) {
    states <- function(mu) {
        S <- mu / (lambda + mu)
        L_per_E <- e / (r2 + omega * lambda + mu)
        E <- lambda * S / ((r1 + e + mu) - omega * lambda * L_per_E)
        c(S = S, E = E, L = E * L_per_E)
    }
    prev_gap <- function(mu) {
        x <- states(mu)
        x[["E"]] + x[["L"]] - prevalence * sum(x)
    }
    lo <- 1e-6
    hi <- 50
    if (prev_gap(lo) * prev_gap(hi) > 0) {
        # No turnover rate reproduces this prevalence; fall back on the ari split.
        return(list(recent = NA_real_, mu = NA_real_))
    }
    mu <- uniroot(prev_gap, c(lo, hi), tol = 1e-10)$root
    x <- states(mu)
    list(recent = x[["E"]] / (x[["E"]] + x[["L"]]), mu = mu)
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
VE <- function(res1, res2, popsize, times = NULL, timestep = 0.1) {
    nt <- ncol(res1)

    # res is cumulative incidence at t = 0, dt, ..., duration. Disease-free
    # person-time is integrated over the nt - 1 intervals between those states
    # (trapezoid), so the denominator spans exactly the same window as the
    # numerator. Summing all nt states would credit one interval too many.
    person_years <- function(res) {
        left <- res[, 1:(nt - 1), drop = FALSE]
        right <- res[, 2:nt, drop = FALSE]
        popsize * rowSums(timestep * (1 - (left + right) / (2 * popsize)))
    }

    PYlist1 <- person_years(res1)
    PYlist2 <- person_years(res2)

    Ilist1 <- res1[, nt]
    Ilist2 <- res2[, nt]

    p.values <- rep(NA_real_, nrow(res1))
    for (i in seq_len(nrow(res1))) {
        p.values[i] <- rateratio.test(
            c(Ilist1[i], Ilist2[i]), c(PYlist1[i], PYlist2[i])
        )$p.value
    }

    return(list(
        vax.eff = 1 - (Ilist1 / PYlist1) / (Ilist2 / PYlist2),
        p.values = p.values
    ))
}
