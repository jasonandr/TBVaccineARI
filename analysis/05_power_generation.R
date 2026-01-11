# Generate Power Analysis Data
# Corresponds to: SampleSizeVE_2025_11_05.R (Generation part)

library(tidyverse)
library(rateratio.test)

source("R/model.R")
source("R/utils.R")

# Load Data
epi_data <- read_csv("data/epi_parameters.csv")
prog_data <- read_csv("data/calibration_base.csv")

# Internal Power Calc Helper
PowerFromModel <- function(sims, popsize, epi_setting, QFTenrol, fastprogon,
                           LeakyPOI, LeakyPOD, AoN.I, AoN.D) {
    out <- VEfromModel(
        sims, popsize, epi_setting, epi_data, prog_data,
        QFTenrol, fastprogon, LeakyPOI, LeakyPOD, AoN.I, AoN.D
    )
    return(length(which(out$p.values < 0.05)) / length(out$p.values))
}

# Settings
sims <- 5000 # Increased for smoother results as requested.
# Full list for proper curves:
pop.list <- c(50, 75, 100, 125, 150, 175, 200, 250, 325, 400, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)

leakyPOD.VE <- c(0, 0.5)
leakyPOI.VE <- c(0, 0.5)
AoN.D.VE <- c(0, 0.5)
AoN.I.VE <- c(0, 0.5)
epi <- c("low", "medium", "high")
QFTenrol <- c("all", "onlyneg", "onlypos")

# Create Grid
df <- expand.grid(pop.list, leakyPOD.VE, leakyPOI.VE, AoN.D.VE, AoN.I.VE, epi, QFTenrol)
names(df) <- c("pop.list", "leakyPOD.VE", "leakyPOI.VE", "AoN.D.VE", "AoN.I.VE", "epi", "QFTenrol")
df$power <- NA

# Filter illogical combinations (as per original code)
df <- df %>%
    filter(!((leakyPOD.VE == 0.5 | leakyPOI.VE == 0.5) & (AoN.D.VE == 0.5 | AoN.I.VE == 0.5))) %>%
    filter(!(leakyPOI.VE == 0.5 & leakyPOD.VE == 0)) %>%
    filter(!(AoN.D.VE == 0 & AoN.I.VE == 0.5))

cat(paste("Generating Power Data. Total rows:", nrow(df), "\n"))

# Loop
for (i in 1:nrow(df)) {
    if (i %% 10 == 0) cat(paste("Processing row", i, "of", nrow(df), "...\n"))

    df$power[i] <- PowerFromModel(
        sims = sims, popsize = df$pop.list[i], epi = as.character(df$epi[i]),
        QFTenrol = as.character(df$QFTenrol[i]),
        fastprogon = 0,
        LeakyPOD = df$leakyPOD.VE[i],
        LeakyPOI = df$leakyPOI.VE[i],
        AoN.D = df$AoN.D.VE[i],
        AoN.I = df$AoN.I.VE[i]
    )
}

if (!dir.exists("data")) dir.create("data")
saveRDS(df, file = "data/power_data.rds")
cat("Saved data/power_data.rds\n")
