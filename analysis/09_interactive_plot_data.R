# 09_interactive_plot_data.R -- the simulated vaccine efficacy distributions
# underlying Figure 2, as a flat .csv for the journal's Interactive Plot Viewer.
#
# Same draws as Figure 2 in 08_manuscript_figures.R: base model, mixed cohort,
# 1500 participants per arm, 5000 simulations per archetype and setting.
#
# Run from the repository root:  Rscript analysis/09_interactive_plot_data.R
# Environment:
#   FIG_OUT  output directory (default figures/manuscript/)
suppressMessages({library(tidyverse); library(rateratio.test)})
source("R/model.R"); source("R/utils.R")

OUT <- Sys.getenv("FIG_OUT", unset = "figures/manuscript/")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

epi   <- read_csv("data/epi_parameters.csv", show_col_types = FALSE)
progB <- read_csv("data/calibration_base.csv", show_col_types = FALSE)

SET    <- c("low", "medium", "high")
SETLAB <- c(low = "Medium", medium = "High", high = "Very high")
SIMS   <- 5000

# AoN.I, AoN.D, leakyPOD, leakyPOI -- the four archetypes shown in Figure 2
ARCH <- list("AoN POD"       = c(0, .5, 0, 0),
             "Leaky POD"     = c(0, 0, .5, 0),
             "AoN POD/POI"   = c(.5, .5, 0, 0),
             "Leaky POD/POI" = c(0, 0, .5, .5))

draws <- function(n, s, q, p, prog, sims = SIMS, seed = 7) {
    set.seed(seed)
    r <- VEfromModel(sims, n, s, epi, prog, q, 0, p[4], p[3], p[1], p[2])
    r$vax.eff[is.finite(r$vax.eff)]
}

cat("collecting Figure 2 draws (mixed cohort, 1500 per arm)...\n")
out <- map_dfr(names(ARCH), function(a)
    map_dfr(SET, function(s) {
        v <- draws(1500, s, "all", ARCH[[a]], progB)
        tibble(archetype = a,
               ari_setting = SETLAB[[s]],
               simulation = seq_along(v),
               vaccine_efficacy_percent = round(100 * v, 3))
    }))

out$ari_setting <- factor(out$ari_setting, levels = c("Medium", "High", "Very high"))
out <- out %>% arrange(match(archetype, names(ARCH)), ari_setting, simulation)

write_csv(out, paste0(OUT, "Figure2_interactive_plot_data.csv"))
cat("wrote ", OUT, "Figure2_interactive_plot_data.csv (", nrow(out), " rows)\n", sep = "")

out %>% group_by(archetype, ari_setting) %>%
    summarise(median = round(median(vaccine_efficacy_percent), 1), .groups = "drop") %>%
    pivot_wider(names_from = ari_setting, values_from = median) %>%
    as.data.frame() %>% print(row.names = FALSE)
