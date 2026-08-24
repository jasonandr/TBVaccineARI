# 07_power_curves.R -- power to detect a difference in incidence between arms,
# across sample size, vaccine archetype, cohort and epidemiologic setting.
#
# Produces data/power_curves.csv, which 08_manuscript_figures.R plots as Figure 3,
# plus data/power_sample_sizes.csv giving the sample size per arm at 90% power
# (linear interpolation between adjacent grid points).
#
# Run from the repository root:  Rscript analysis/07_power_curves.R
suppressMessages({library(tidyverse); library(rateratio.test)})
source("R/model.R"); source("R/utils.R")

epi   <- read_csv("data/epi_parameters.csv", show_col_types = FALSE)
progB <- read_csv("data/calibration_base.csv", show_col_types = FALSE)

SET  <- c("low", "medium", "high")          # labelled Medium / High / Very high
COH  <- c("onlypos", "all", "onlyneg")
SIMS <- 5000
SEED <- 11
ALPHA <- 0.05

POPS <- c(50, 75, 100, 125, 150, 175, 200, 250, 325, 400, 500, 750, 1000,
          1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)

# AoN.I, AoN.D, leakyPOD, leakyPOI
ARCH <- list("AoN POD"       = c(0, .5, 0, 0),
             "Leaky POD"     = c(0, 0, .5, 0),
             "AoN POD/POI"   = c(.5, .5, 0, 0),
             "Leaky POD/POI" = c(0, 0, .5, .5),
             "Leaky POI"     = c(0, 0, 0, .5))

power_at <- function(n, s, q, p) {
    set.seed(SEED)
    r <- VEfromModel(SIMS, n, s, epi, progB, q, 0, p[4], p[3], p[1], p[2])
    mean(r$p.values < ALPHA, na.rm = TRUE)
}

grid <- expand_grid(archetype = names(ARCH), cohort = COH, setting = SET, n = POPS)
cat("running", nrow(grid), "power cells at", SIMS, "simulations each...\n")
t0 <- Sys.time()
pw <- grid %>%
    mutate(power = pmap_dbl(list(archetype, cohort, setting, n),
                            function(a, q, s, N) power_at(N, s, q, ARCH[[a]])))
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

pw$config <- "base"
write_csv(pw %>% select(config, archetype, cohort, setting, n, power),
          "data/power_curves.csv")
cat("wrote data/power_curves.csv\n")

# Sample size per arm at 90% power, interpolated between bracketing grid points.
n90 <- pw %>%
    arrange(archetype, cohort, setting, n) %>%
    group_by(archetype, cohort, setting) %>%
    summarise(n_90 = {
        i <- which(power >= 0.90)[1]
        if (is.na(i)) NA_real_
        else if (i == 1) n[1]
        else {
            x0 <- n[i - 1]; x1 <- n[i]; y0 <- power[i - 1]; y1 <- power[i]
            round(x0 + (0.90 - y0) * (x1 - x0) / (y1 - y0))
        }
    }, .groups = "drop")
write_csv(n90, "data/power_sample_sizes.csv")
cat("wrote data/power_sample_sizes.csv\n\n")

SETLAB <- c(low = "Medium", medium = "High", high = "Very high")
COHLAB <- c(onlypos = "QFT+", all = "Mixed", onlyneg = "QFT-")
n90 %>%
    mutate(setting = SETLAB[setting], cohort = COHLAB[cohort]) %>%
    pivot_wider(names_from = setting, values_from = n_90) %>%
    as.data.frame() %>% print(row.names = FALSE)
