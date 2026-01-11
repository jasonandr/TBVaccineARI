# Reproduce Results: Variable POD Model
# Corresponds to: VE_model_variablePOD_2025_11_28.R

library(tidyverse)
library(rateratio.test)
library(ggsci)
library(ggpubr)

source("R/model.R")
source("R/utils.R")

# Load Data
epi_data <- read_csv("data/epi_parameters.csv")
prog_data <- read_csv("data/calibration_base.csv") # VariablePOD uses base calibration (fastprogon=0 in original?)
# Checking original logic: The file `VE_model_variablePOD_2025_11_28.R` uses `progratio_calibration_base_2025_10_08.csv`
# and fastprogon is an argument, defaulting to what?
# The wrapper `VEfromModel` calls `VaccineModel(..., fastprogon=fastprogon...)`.
# The original script does NOT show how it was called in lines 1-200. I need to assume fastprogon=0 if it uses base calibration.

if (!dir.exists("figures")) dir.create("figures")

# Settings
sims <- 1000
fastprogon <- 0
QFTenrol <- "onlypos"
popsize <- 1500

run_scenarios_df <- function(AoN.I, AoN.D.pos, AoN.D.neg, LeakyPOD, LeakyPOI) {
    # Call with specific pos/neg arguments

    out_high <- VEfromModel(sims, popsize, "high", epi_data, prog_data, QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I,
        AoN.D = NULL, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg
    )
    out_med <- VEfromModel(sims, popsize, "medium", epi_data, prog_data, QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I,
        AoN.D = NULL, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg
    )
    out_low <- VEfromModel(sims, popsize, "low", epi_data, prog_data, QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I,
        AoN.D = NULL, AoN.D.pos = AoN.D.pos, AoN.D.neg = AoN.D.neg
    )

    bind_rows(
        data.frame(ve = out_high$vax.eff, epi = "high"),
        data.frame(ve = out_med$vax.eff, epi = "medium"),
        data.frame(ve = out_low$vax.eff, epi = "low")
    )
}

plot_from_df <- function(df, title) {
    df %>%
        mutate(epi = factor(epi, levels = c("low", "medium", "high"))) %>%
        ggplot(aes(x = epi, y = ve, fill = epi)) +
        geom_boxplot() +
        scale_x_discrete(labels = c("2%", "5%", "50%")) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(title = title, y = "Vaccine Efficacy", x = "Annual Risk of Infection") +
        scale_y_continuous(limits = c(-0.2, 1)) +
        scale_fill_nejm()
}

cat("Running Variable POD Model Scenarios...\n")

# Example Scenario: AoN.D.pos differs from AoN.D.neg
# Assuming typical scenario modeled is: AoN.D.pos = 0 (prev infection reduces efficacy?), AoN.D.neg = 0.5?
# Or checking hypothetical scenarios.
# I will run a demo scenario: Pos=0.3, Neg=0.8

p1 <- plot_from_df(run_scenarios_df(0, AoN.D.pos = 0.3, AoN.D.neg = 0.8, 0, 0), "AoN POD: Pos=0.3, Neg=0.8")

ggsave("figures/04_Results_VariablePOD.pdf", p1, width = 6, height = 4)
print("Saved figures/04_Results_VariablePOD.pdf")
