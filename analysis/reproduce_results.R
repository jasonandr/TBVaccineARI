# Reproduce Vaccine Efficacy Bias Results

library(tidyverse)
library(rateratio.test)
library(ggsci)
library(ggpubr)

# Source model logic
source("R/model.R")
source("R/utils.R")

# Load data
epi_data <- read_csv("data/epi_parameters.csv")
prog_data <- read_csv("data/calibration_results.csv")

# Ensure output directory exists
if (!dir.exists("figures")) dir.create("figures")

# Simulation Settings
# Note: Reduced sims to 1000 for demonstration. Increase to 5000 for full paper reproduction.
sims <- 1000
popsize <- 1500
fastprogon <- 1
QFTenrol <- "onlypos" # "Result 1" mainly focuses on QFT+ (or check original script)
# Original script had sections for "onlypos" and "all".
# We will focus on reproducing the "onlypos" figure (Supplemental Figure equivalent)

run_scenarios <- function(AoN.I, AoN.D, LeakyPOD, LeakyPOI, title_suffix) {
    results_list <- list()

    for (epi_setting in c("high", "medium", "low")) {
        cat(paste("Running:", title_suffix, "-", epi_setting, "...\n"))

        out <- VEfromModel(
            sims = sims, popsize = popsize,
            epi_setting = epi_setting, epi_data = epi_data, prog_data = prog_data,
            QFTenrol = QFTenrol, fastprogon = fastprogon,
            LeakyPOI = LeakyPOI, LeakyPOD = LeakyPOD, AoN.I = AoN.I, AoN.D = AoN.D
        )

        df <- data.frame(ve = out$vax.eff, epi = epi_setting)
        results_list[[epi_setting]] <- df
    }

    bind_rows(results_list)
}

plot_results <- function(data, title) {
    data %>%
        mutate(epi = factor(epi, levels = c("low", "medium", "high"))) %>%
        ggplot(aes(x = epi, y = ve, fill = epi)) +
        geom_boxplot() +
        scale_x_discrete(labels = c("2%", "5%", "50%")) +
        xlab("Annual Risk of Infection") +
        theme_classic() +
        theme(legend.position = "none") +
        ylab("Vaccine Efficacy") +
        scale_y_continuous(limits = c(-0.2, 1)) + # Extended lower limit slightly for visibility
        scale_fill_nejm() +
        labs(title = title)
}

# 1. POD Vaccine with All-or-None Protection
res1 <- run_scenarios(AoN.I = 0, AoN.D = 0.5, LeakyPOD = 0, LeakyPOI = 0, "AoN POD")
p1 <- plot_results(res1, "A. POD Vaccine (AoN)")

# 2. POD Vaccine with Leaky Protection
res2 <- run_scenarios(AoN.I = 0, AoN.D = 0, LeakyPOD = 0.5, LeakyPOI = 0, "Leaky POD")
p2 <- plot_results(res2, "B. POD Vaccine (Leaky)")

# 3. POD/POI Vaccine with All-or-None Protection
res3 <- run_scenarios(AoN.I = 0.5, AoN.D = 0.5, LeakyPOD = 0, LeakyPOI = 0, "AoN POD/POI")
p3 <- plot_results(res3, "C. POD/POI Vaccine (AoN)")

# 4. POD/POI Vaccine with Leaky Protection
res4 <- run_scenarios(AoN.I = 0, AoN.D = 0, LeakyPOD = 0.5, LeakyPOI = 0.5, "Leaky POD/POI")
p4 <- plot_results(res4, "D. POD/POI Vaccine (Leaky)")

# Combine and Save
comb.fig <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggexport(comb.fig, width = 10, height = 6, filename = "figures/VE_Bias_HighTransmission.pdf")

cat("Plot saved to figures/VE_Bias_HighTransmission.pdf\n")
