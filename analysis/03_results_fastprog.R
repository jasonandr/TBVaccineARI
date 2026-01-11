# TB Vaccine Efficacy Analysis: Fast Progression Model
# This script runs the analysis for the model variant including fast progression.

library(tidyverse)
library(rateratio.test)
library(ggsci)
library(ggpubr)

source("R/model.R")
source("R/utils.R")

# Load Data
epi_data <- read_csv("data/epi_parameters.csv")
prog_data <- read_csv("data/calibration_fastprog.csv")

if (!dir.exists("figures")) dir.create("figures")

# Settings
sims <- 1000
fastprogon <- 1
popsize <- 1500

# Helper Scenarios Runner
run_scenarios_df <- function(QFTenrol, AoN.I, AoN.D, LeakyPOD, LeakyPOI) {
    out_high <- VEfromModel(sims, popsize, "high", epi_data, prog_data, QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I, AoN.D)
    out_med <- VEfromModel(sims, popsize, "medium", epi_data, prog_data, QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I, AoN.D)
    out_low <- VEfromModel(sims, popsize, "low", epi_data, prog_data, QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I, AoN.D)

    bind_rows(
        data.frame(ve = out_high$vax.eff, epi = "high"),
        data.frame(ve = out_med$vax.eff, epi = "medium"),
        data.frame(ve = out_low$vax.eff, epi = "low")
    )
}

# Plotting Function
plot_from_df <- function(df, title) {
    df %>%
        mutate(epi = factor(epi, levels = c("low", "medium", "high"))) %>%
        ggplot(aes(x = epi, y = ve, fill = epi)) +
        geom_boxplot() +
        scale_x_discrete(labels = c("Medium", "High", "Very High")) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(title = title, y = "Vaccine Efficacy", x = "Annual Risk of Infection") +
        scale_y_continuous(limits = c(-0.2, 1)) +
        scale_fill_nejm()
}

# Main Loop for Cohorts
cohorts <- c("onlypos", "all", "onlyneg")
names(cohorts) <- c("QFT+", "All", "QFT-")

for (cohort_name in names(cohorts)) {
    qft_code <- cohorts[cohort_name]
    cat(paste0("Running FastProg Model Scenarios for ", cohort_name, "...\n"))

    # 1. Scenarios
    df1 <- run_scenarios_df(qft_code, 0, 0.5, 0, 0) # AoN POD
    df2 <- run_scenarios_df(qft_code, 0, 0, 0.5, 0) # Leaky POD
    df3 <- run_scenarios_df(qft_code, 0.5, 0.5, 0, 0) # AoN POD/POI
    df4 <- run_scenarios_df(qft_code, 0, 0, 0.5, 0.5) # Leaky POD/POI

    # 2. Figures
    p1 <- plot_from_df(df1, "A. AoN POD")
    p2 <- plot_from_df(df2, "B. Leaky POD")
    p3 <- plot_from_df(df3, "C. AoN POD/POI")
    p4 <- plot_from_df(df4, "D. Leaky POD/POI")

    comb.fig <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
    filename <- paste0("figures/03_FastProg_VE_", qft_code, ".pdf")
    ggsave(filename, comb.fig, width = 10, height = 6)
    cat(paste("Saved", filename, "\n"))

    # 3. Text Results (Summaries)
    cat(paste0("\n--- Text Summary for ", cohort_name, " ---\n"))

    cat("Scenario: Leaky POD/POI\n")
    print(df4 %>% group_by(epi) %>%
        summarise(
            median = median(ve, na.rm = TRUE),
            lower = quantile(ve, 0.025, na.rm = TRUE),
            upper = quantile(ve, 0.975, na.rm = TRUE)
        ))
    cat("\n")
}
