# 08_manuscript_figures.R -- reproduces the result figures in the manuscript:
#   Figure 2      observed VE, base model, mixed cohort
#   Figure 3      power curves (needs data/power_curves.csv from 07_power_curves.R)
#   Figure S2/S3  observed VE, base model, QFT+ / QFT-
#   Figure S4     power with reduced efficacy in QFT-negative participants
#   Figure S5     POI-only vaccine across the three cohorts
#   Figure S6     trial-size sensitivity at very high ARI
#   Figure S7-S9  observed VE, alternative model, mixed / QFT- / QFT+
# and Supplemental Table S1 (expected efficacy by archetype, cohort and setting).
#
# Figures 1 and S1 are model schematics, drawn by analysis/make_schematics.py.
#
# Run from the repository root, after 07_power_curves.R:
#   Rscript analysis/08_manuscript_figures.R
# Environment:
#   FIG_OUT  output directory              (default figures/manuscript/)
#   FIG_FMT  "png" or "pdf"; pdf is vector (default png)
suppressMessages({library(tidyverse); library(rateratio.test); library(ggsci); library(ggpubr)})
source("R/model.R"); source("R/utils.R")

OUT <- Sys.getenv("FIG_OUT", unset = "figures/manuscript/")
FMT <- Sys.getenv("FIG_FMT", unset = "png")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fig <- function(stem) paste0(OUT, stem, ".", FMT)

epi   <- read_csv("data/epi_parameters.csv", show_col_types = FALSE)
progB <- read_csv("data/calibration_base.csv", show_col_types = FALSE)
progF <- read_csv("data/calibration_fastprog.csv", show_col_types = FALSE)

SET    <- c("low", "medium", "high")
SETLAB <- c(low = "Medium", medium = "High", high = "Very high")
COH    <- c(onlypos = "QFT+", all = "Mixed", onlyneg = "QFT-")
SIMS   <- 5000

# AoN.I, AoN.D, leakyPOD, leakyPOI
ARCH <- list("AoN POD"       = c(0, .5, 0, 0),
             "Leaky POD"     = c(0, 0, .5, 0),
             "AoN POD/POI"   = c(.5, .5, 0, 0),
             "Leaky POD/POI" = c(0, 0, .5, .5),
             "Leaky POI"     = c(0, 0, 0, .5))

# The early/late split of prevalent infection is the stationary solution
# (VaccineModel's default); see R/model.R.
draws <- function(n, s, q, p, fp, prog, sims = SIMS, seed = 7) {
    set.seed(seed)
    r <- VEfromModel(sims, n, s, epi, prog, q, fp, p[4], p[3], p[1], p[2])
    r$vax.eff[is.finite(r$vax.eff)]
}

# ---------------------------------------------------------------- collect data
# base = natural history as published; alt = additional reinfection-driven
# progression out of the early infection state (fastprogon = 1).
CONF <- list(base = list(progB, 0), alt = list(progF, 1))

cat("collecting VE distributions...\n")
ve <- map_dfr(names(CONF), function(cf) {
  z <- CONF[[cf]]
  map_dfr(names(ARCH), function(a)
    map_dfr(names(COH), function(q)
      map_dfr(SET, function(s)
        tibble(config = cf, archetype = a, cohort = q, setting = s,
               ve = draws(1500, s, q, ARCH[[a]], z[[2]], z[[1]])))))
})

cat("collecting trial-size distributions...\n")
ss <- map_dfr(names(ARCH), function(a)
  map_dfr(names(COH), function(q)
    map_dfr(c(1500, 450, 300), function(n)
      tibble(archetype = a, cohort = q, n = n,
             ve = draws(n, "high", q, ARCH[[a]], 0, progB)))))

# Expected ("asymptotic") efficacy from large cohorts, so that sampling
# variability and the finite-sample behaviour of 1 - IRR are removed. The
# medium-ARI value is the low-saturation benchmark drawn as a dashed line.
cat("computing expected-efficacy reference...\n")
ref <- map_dfr(names(CONF), function(cf) {
  z <- CONF[[cf]]
  map_dfr(names(ARCH), function(a)
    map_dfr(names(COH), function(q)
      map_dfr(SET, function(s)
        tibble(config = cf, archetype = a, cohort = q, setting = s,
               asym = 100 * mean(draws(50000, s, q, ARCH[[a]], z[[2]], z[[1]],
                                       sims = 400))))))
})
write_csv(ref, paste0(OUT, "TableS1_expected_efficacy.csv"))

# ------------------------------------------------------------------- plotting
theme_ms <- function() {
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "plain"),
          axis.title = element_text(size = 9),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
}

box_panel <- function(d, title, refline = NA) {
    p <- d %>% mutate(setting = factor(setting, levels = SET)) %>%
        ggplot(aes(x = setting, y = 100 * ve, fill = setting))
    if (!is.na(refline)) {
        p <- p + geom_hline(yintercept = refline, linetype = "dashed",
                            colour = "grey35", linewidth = 0.4)
    }
    p + geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25, linewidth = 0.3) +
        scale_x_discrete(labels = SETLAB) +
        scale_y_continuous(limits = c(-20, 100), breaks = seq(-20, 100, 20)) +
        scale_fill_nejm() + theme_ms() +
        labs(title = title, y = "Observed vaccine efficacy (%)",
             x = "Annual risk of infection")
}

ve_figure <- function(cohort, config, archs, stem, width = 9, height = 6, tag = TRUE) {
    d <- ve %>% filter(cohort == !!cohort, config == !!config, archetype %in% archs)
    rr <- ref %>% filter(cohort == !!cohort, config == !!config, setting == "low")
    panels <- imap(archs, function(a, i) {
        rl <- rr$asym[rr$archetype == a]
        box_panel(d %>% filter(archetype == a),
                  paste0(if (tag) paste0(LETTERS[i], ". ") else "", a),
                  if (length(rl)) rl else NA)
    })
    g <- ggarrange(plotlist = panels, nrow = ceiling(length(archs) / 2),
                   ncol = min(2, length(archs)))
    ggsave(fig(stem), g, width = width, height = height, dpi = 300)
    cat("  wrote", basename(fig(stem)), "\n")
}

MAIN4 <- c("AoN POD", "Leaky POD", "AoN POD/POI", "Leaky POD/POI")
cat("drawing figures...\n")
ve_figure("all",     "base", MAIN4, "Figure2")
ve_figure("onlypos", "base", MAIN4, "FigureS2")
ve_figure("onlyneg", "base", MAIN4, "FigureS3")
ve_figure("all",     "alt",  MAIN4, "FigureS7")
ve_figure("onlyneg", "alt",  MAIN4, "FigureS8")
ve_figure("onlypos", "alt",  MAIN4, "FigureS9")

# S5: POI-only across the three cohorts
poi <- ve %>% filter(config == "base", archetype == "Leaky POI")
pr <- ref %>% filter(config == "base", archetype == "Leaky POI", setting == "low")
panels <- imap(names(COH), function(q, i)
    box_panel(poi %>% filter(cohort == q),
              paste0(LETTERS[i], ". ", COH[[q]]), pr$asym[pr$cohort == q]))
ggsave(fig("FigureS5"), ggarrange(plotlist = panels, nrow = 1, ncol = 3),
       width = 11, height = 3.4, dpi = 300)
cat("  wrote", basename(fig("FigureS5")), "\n")

# S6: trial-size sensitivity at very high ARI
s6 <- ss %>% filter(archetype %in% MAIN4) %>%
    mutate(n = factor(n, levels = c(1500, 450, 300)),
           cohort = factor(COH[cohort], levels = c("QFT+", "Mixed", "QFT-")))
g6 <- ggplot(s6, aes(x = n, y = 100 * ve, fill = n)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, linewidth = 0.3) +
    facet_grid(cohort ~ archetype) +
    scale_y_continuous(limits = c(-40, 100), breaks = seq(-40, 100, 20)) +
    scale_fill_nejm() + theme_ms() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9)) +
    labs(y = "Observed vaccine efficacy (%)", x = "Participants per arm")
ggsave(fig("FigureS6"), g6, width = 10, height = 6.5, dpi = 300)
cat("  wrote", basename(fig("FigureS6")), "\n")

# Figure 3: power curves
pwfile <- "data/power_curves.csv"
if (!file.exists(pwfile)) {
    stop("missing ", pwfile, " -- run analysis/07_power_curves.R first")
}
pw <- read_csv(pwfile, show_col_types = FALSE) %>%
      filter(config == "base", archetype %in% MAIN4)
pw$cohort <- factor(COH[pw$cohort], levels = c("QFT+", "Mixed", "QFT-"))
pw$archetype <- factor(pw$archetype, levels = MAIN4)
pw$setting <- factor(pw$setting, levels = SET, labels = SETLAB)
g3 <- ggplot(pw, aes(x = n, y = 100 * power, colour = setting)) +
    geom_hline(yintercept = 90, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    facet_grid(cohort ~ archetype) +
    scale_x_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 2500)) +
    scale_colour_nejm(name = "Annual risk\nof infection") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right", strip.background = element_blank(),
          strip.text = element_text(size = 9),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)) +
    labs(x = "Sample size per arm", y = "Power (%)")
ggsave(fig("Figure3"), g3, width = 10, height = 6.5, dpi = 300)
cat("  wrote", basename(fig("Figure3")), "\n")

# S4: reduced efficacy in QFT-negative participants
cat("running variable-POD power...\n")
POPS <- c(50, 100, 200, 325, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000)
vp <- expand_grid(neg = c(0.5, 0.335, 0.25), cohort = c("onlypos", "all"),
                  setting = SET, n = POPS) %>%
    mutate(power = pmap_dbl(list(neg, cohort, setting, n), function(nn, q, s, N) {
        set.seed(11)
        r <- VEfromModel(2000, N, s, epi, progB, q, 0, 0, 0, 0, NULL,
                         AoN.D.pos = 0.5, AoN.D.neg = nn)
        mean(r$p.values < 0.05, na.rm = TRUE)
    }))
vp$cohort <- factor(COH[vp$cohort], levels = c("QFT+", "Mixed"))
vp$setting <- factor(vp$setting, levels = SET, labels = SETLAB)
vp$neg <- factor(vp$neg, levels = c(0.5, 0.335, 0.25),
                 labels = c("Equal (50%)", "33% lower", "50% lower"))
g4 <- ggplot(vp, aes(x = n, y = 100 * power, colour = setting)) +
    geom_hline(yintercept = 90, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + facet_grid(neg ~ cohort) +
    scale_x_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 2500)) +
    scale_colour_nejm(name = "Annual risk\nof infection") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right", strip.background = element_blank(),
          strip.text = element_text(size = 9),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)) +
    labs(x = "Sample size per arm", y = "Power (%)")
ggsave(fig("FigureS4"), g4, width = 8, height = 6.5, dpi = 300)
cat("  wrote", basename(fig("FigureS4")), "\nDONE\n")
