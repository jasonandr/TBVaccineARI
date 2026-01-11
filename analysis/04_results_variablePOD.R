# TB Vaccine Efficacy Analysis: Variable POD Model
# This script runs the analysis for the model variant with variable protection against disease.
# It reproduces the Power vs Sample Size analysis for different QFT cohorts and VE assumptions.

library(tidyverse)
library(rateratio.test)
library(ggsci)
library(gridExtra)

source("R/model.R")
source("R/utils.R")

# Load Data
epi_data <- read_csv("data/epi_parameters.csv", show_col_types = FALSE)
prog_data <- read_csv("data/calibration_base.csv", show_col_types = FALSE)

# Check output directory
if (!dir.exists("figures")) dir.create("figures")

# --- Helper Function ---
# Calculates Power (proportion of significant results) from the model
PowerFromModel <- function(sims, popsize, epi, QFTenrol, fastprogon,
                           LeakyPOI, LeakyPOD, AoN.I, AoN.D.pos, AoN.D.neg) {
    # Call the unified VEfromModel (which handles parameter lookup using strict data passing)
    out <- VEfromModel(
        sims = sims,
        popsize = popsize,
        epi_setting = epi, # "low", "medium", "high"
        epi_data = epi_data,
        prog_data = prog_data,
        QFTenrol = QFTenrol,
        fastprogon = fastprogon,
        LeakyPOI = LeakyPOI,
        LeakyPOD = LeakyPOD,
        AoN.I = AoN.I,
        AoN.D.pos = AoN.D.pos,
        AoN.D.neg = AoN.D.neg
    )

    # Calculate power: proportion of p-values < 0.05
    power <- length(which(out$p.values < 0.05)) / length(out$p.values)
    return(power)
}

# --- Settings ---
sims <- 5000 # High precision as per user script
# p.s. This takes a long time. For testing, you might want to reduce this,
# but for reproduction we use 5000.

# Population Sizes (Full list from user script)
pop.list <- c(
    50, 75, 100, 125, 150, 175, 200, 250, 325, 400, 500, 750,
    1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000
)

# Parameters
AoN.D.pos.VE <- 0.5
AoN.D.neg.VE <- c(0.5 * 0.5, (2 / 3) * 0.5, 0.5) # 0.25, 0.333..., 0.5
epi <- c("low", "medium", "high")
QFTenrol <- c("all", "onlypos")

# Create Grid
df <- expand.grid(pop.list, AoN.D.pos.VE, AoN.D.neg.VE, epi, QFTenrol)
names(df) <- c("pop.list", "AoN.D.pos.VE", "AoN.D.neg.VE", "epi", "QFTenrol")
df$power <- NA

# --- Run Scenarios ---
cat("Running Variable POD Power Analysis. Total scenarios:", nrow(df), "\n")
cat("Note: This may take a significant amount of time due to 5000 simulations per scenario.\n")

for (i in seq_len(nrow(df))) {
    if (i %% 10 == 0) cat(paste(i, "of", nrow(df), "\n"))

    df$power[i] <- PowerFromModel(
        sims = sims,
        popsize = df$pop.list[i],
        epi = as.character(df$epi[i]),
        QFTenrol = as.character(df$QFTenrol[i]),
        fastprogon = 0,
        LeakyPOD = 0,
        LeakyPOI = 0,
        AoN.D.pos = df$AoN.D.pos.VE[i],
        AoN.D.neg = df$AoN.D.neg.VE[i],
        AoN.I = 0
    )
}

# Save intermediate results (optional but good practice)
# saveRDS(df, file = "data/power_variable_pod.rds")


# --- Plotting ---
titlesize <- 10

# Helper to keep code DRYer? The user script separates them explicitly.
# We will follow the user's structure to match the output exactly.

# Plot 1: QFT+, QFT- reduced 50%
p1 <- df %>%
    filter(QFTenrol == "onlypos") %>%
    filter(abs(AoN.D.neg.VE - 0.25) < 0.01) %>% # Use tolerance for float comparison
    ggplot(aes(x = pop.list, y = power, color = epi)) +
    geom_line() +
    theme_classic() +
    xlab("Sample Size per Arm") +
    scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
    theme(legend.position = c(0.8, 0.4)) +
    labs(title = "AoN POD Vaccine (VE reduced 50% in QFT-), QFT+") +
    theme(
        plot.title = element_text(size = titlesize),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.spacing = unit(0.25, "lines")
    ) +
    theme(legend.position = "none")

# Plot 2: QFT+, QFT- reduced 33% (approx 0.333...)
p2 <- df %>%
    filter(QFTenrol == "onlypos") %>%
    filter(AoN.D.neg.VE > 0.3 & AoN.D.neg.VE < 0.4) %>%
    ggplot(aes(x = pop.list, y = power, color = epi)) +
    geom_line() +
    theme_classic() +
    xlab("Sample Size per Arm") +
    scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
    theme(legend.position = c(0.8, 0.4)) +
    labs(title = "AoN POD Vaccine (VE reduced 33% in QFT-), QFT+") +
    theme(
        plot.title = element_text(size = titlesize),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.spacing = unit(0.25, "lines")
    ) +
    theme(legend.position = "none")

# Plot 3: QFT+, QFT- No reduction (0.5)
p3 <- df %>%
    filter(QFTenrol == "onlypos") %>%
    filter(abs(AoN.D.neg.VE - 0.5) < 0.01) %>%
    ggplot(aes(x = pop.list, y = power, color = epi)) +
    geom_line() +
    theme_classic() +
    xlab("Sample Size per Arm") +
    scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
    theme(legend.position = c(0.8, 0.4)) +
    labs(title = "AoN POD Vaccine (VE not reduced in QFT-), QFT+") +
    theme(
        plot.title = element_text(size = titlesize),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.spacing = unit(0.25, "lines")
    )

# Plot 4: All, QFT- reduced 50%
p4 <- df %>%
    filter(QFTenrol == "all") %>%
    filter(abs(AoN.D.neg.VE - 0.25) < 0.01) %>%
    ggplot(aes(x = pop.list, y = power, color = epi)) +
    geom_line() +
    theme_classic() +
    xlab("Sample Size per Arm") +
    scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
    theme(legend.position = c(0.8, 0.4)) +
    labs(title = "AoN POD Vaccine (VE reduced 50% in QFT-), All QFT") +
    theme(
        plot.title = element_text(size = titlesize),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.spacing = unit(0.25, "lines")
    ) +
    theme(legend.position = "none")

# Plot 5: All, QFT- reduced 33%
p5 <- df %>%
    filter(QFTenrol == "all") %>%
    filter(AoN.D.neg.VE > 0.3 & AoN.D.neg.VE < 0.4) %>%
    ggplot(aes(x = pop.list, y = power, color = epi)) +
    geom_line() +
    theme_classic() +
    xlab("Sample Size per Arm") +
    scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
    theme(legend.position = c(0.8, 0.4)) +
    labs(title = "AoN POD Vaccine (VE reduced 33% in QFT-), All QFT") +
    theme(
        plot.title = element_text(size = titlesize),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.spacing = unit(0.25, "lines")
    ) +
    theme(legend.position = "none")

# Plot 6: All, QFT- No reduction
p6 <- df %>%
    filter(QFTenrol == "all") %>%
    filter(abs(AoN.D.neg.VE - 0.5) < 0.01) %>%
    ggplot(aes(x = pop.list, y = power, color = epi)) +
    geom_line() +
    theme_classic() +
    xlab("Sample Size per Arm") +
    scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
    theme(legend.position = c(0.8, 0.4)) +
    labs(title = "AoN POD Vaccine (VE not reduced in QFT-), All QFT") +
    theme(
        plot.title = element_text(size = titlesize),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.spacing = unit(0.25, "lines")
    ) +
    theme(legend.position = "none")

# Combine Logic (Rows: 3)
# Row 1: p3, p6 (No reduction)
# Row 2: p2, p5 (33% reduction)
# Row 3: p1, p4 (50% reduction)
comb.fig <- grid.arrange(p3, p6, p2, p5, p1, p4, nrow = 3)

# Save
filename <- "figures/04_Variable_VE_Power.pdf"
ggsave(filename, comb.fig, width = 10, height = 6)
cat(paste("Saved", filename, "\n"))
