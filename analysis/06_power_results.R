# Analyze Power Results
# This script generates power curves and calculates required sample sizes based on the simulated data.

library(tidyverse)
library(ggsci)
library(gridExtra)
library(mgcv) # For GAM output in PowerToSample

# Load Data
if (!file.exists("data/power_data.rds")) {
    stop("data/power_data.rds used not found. Please run 'analysis/05_power_generation.R' first.")
}
df <- readRDS("data/power_data.rds")

# Anchor curves at 0 (Power = 0 at Sample Size = 0)
# Create a 0-sample row for every scenario
scenarios <- df %>%
    select(-pop.list, -power) %>%
    distinct() %>%
    mutate(pop.list = 0, power = 0)

df <- bind_rows(df, scenarios) %>%
    arrange(pop.list)

if (!dir.exists("figures")) dir.create("figures")

# 1. Figures (Power Curves)
titlesize <- 10

# Definition of Plots (Mapped from original code p1-p12)
plot_power <- function(data, title, show_legend = FALSE) {
    p <- data %>%
        ggplot(aes(x = pop.list, y = power, color = epi)) +
        geom_line() +
        theme_classic() +
        xlab("Sample Size per Arm") +
        scale_x_continuous(limits = c(0, 5000)) +
        scale_color_nejm(name = "ARI", labels = c("Medium", "High", "Very high")) +
        labs(title = title) +
        theme(plot.title = element_text(size = titlesize))

    if (show_legend) {
        p <- p + theme(
            legend.position = c(0.75, 0.4), # Inside Top-Rightish
            legend.background = element_rect(fill = NA)
        )
    } else {
        p <- p + theme(legend.position = "none")
    }
    return(p)
}

# Helper to filter scenarios
# AoN POD: D=0.5, I=0, L_POD=0, L_POI=0
filter_aon_pod <- function(d) filter(d, AoN.D.VE == 0.5 & AoN.I.VE == 0 & leakyPOD.VE == 0 & leakyPOI.VE == 0)
filter_aon_pp <- function(d) filter(d, AoN.D.VE == 0.5 & AoN.I.VE == 0.5 & leakyPOD.VE == 0 & leakyPOI.VE == 0)
filter_lky_pod <- function(d) filter(d, AoN.D.VE == 0 & AoN.I.VE == 0 & leakyPOD.VE == 0.5 & leakyPOI.VE == 0)
filter_lky_pp <- function(d) filter(d, AoN.D.VE == 0 & AoN.I.VE == 0 & leakyPOD.VE == 0.5 & leakyPOI.VE == 0.5)

# Generate plots for each QFT cohort
plots <- list()
counter <- 1
first_plot <- TRUE # Only show legend on the very first plot

for (qft in c("onlypos", "all", "onlyneg")) {
    # Mapped title for display
    qft_disp <- dplyr::case_when(
        qft == "onlypos" ~ "QFT+",
        qft == "all" ~ "All QFT",
        qft == "onlyneg" ~ "QFT-"
    )

    sub <- df %>% filter(QFTenrol == qft)

    plots[[counter]] <- plot_power(filter_aon_pod(sub), paste("AoN POD Vaccine,", qft_disp), show_legend = first_plot)
    first_plot <- FALSE
    counter <- counter + 1
    plots[[counter]] <- plot_power(filter_aon_pp(sub), paste("AoN POD/POI Vaccine,", qft_disp))
    counter <- counter + 1
    plots[[counter]] <- plot_power(filter_lky_pod(sub), paste("Leaky POD Vaccine,", qft_disp))
    counter <- counter + 1
    plots[[counter]] <- plot_power(filter_lky_pp(sub), paste("Leaky POD/POI Vaccine,", qft_disp))
    counter <- counter + 1
}

comb.fig <- do.call(grid.arrange, c(plots, nrow = 3))
ggsave("figures/05_Power_Curves.pdf", comb.fig, width = 12, height = 9)
cat("Saved figures/05_Power_Curves.pdf\n")


# 2. Text Results (Sample Size Calculation)
# Function from PowerAnalysisTextResults...
PowerToSample <- function(dfnew, powertarget, qft, epitype, vaccinetype) {
    df.2 <- dfnew

    # Filter based on type
    if (vaccinetype == "AoN.POD") {
        df.3 <- df.2 %>% filter(AoN.D.VE == 0.5 & leakyPOD.VE == 0 & leakyPOI.VE == 0 & AoN.I.VE == 0 & epi == epitype & QFTenrol == qft)
    } else if (vaccinetype == "AoN.PODPOI") {
        df.3 <- df.2 %>% filter(AoN.D.VE == 0.5 & leakyPOD.VE == 0 & leakyPOI.VE == 0 & AoN.I.VE == 0.5 & epi == epitype & QFTenrol == qft)
    } else if (vaccinetype == "Leaky.POD") {
        df.3 <- df.2 %>% filter(AoN.D.VE == 0 & leakyPOD.VE == 0.5 & leakyPOI.VE == 0 & AoN.I.VE == 0 & epi == epitype & QFTenrol == qft)
    } else if (vaccinetype == "Leaky.PODPOI") {
        df.3 <- df.2 %>% filter(AoN.D.VE == 0 & leakyPOD.VE == 0.5 & leakyPOI.VE == 0.5 & AoN.I.VE == 0 & epi == epitype & QFTenrol == qft)
    } else {
        return(NA)
    }

    if (nrow(df.3) < 4) {
        return(NA)
    } # Not enough points for smooth

    # Bound probabilities for logit
    df.3 <- df.3 %>% mutate(power_p = pmin(pmax(power, 1e-6), 1 - 1e-6))

    m <- gam(power_p ~ s(pop.list, k = min(10, nrow(df.3) - 1)),
        family = quasibinomial(link = "logit"), data = df.3
    )

    grid <- tibble(pop.list = seq(min(df.3$pop.list), max(df.3$pop.list) + 1000, by = 1))
    pred <- predict(m, newdata = grid, type = "response")

    n90 <- grid$pop.list[which(pred >= powertarget)[1]]
    return(n90)
}

cat("\n--- Sample Sizes for 90% Power ---\n")
# Example calls from original script
cat("AoN POD (QFT+, Low):", PowerToSample(df, 0.9, "onlypos", "low", "AoN.POD"), "\n")
cat("AoN POD (QFT+, High):", PowerToSample(df, 0.9, "onlypos", "high", "AoN.POD"), "\n")
cat("Leaky POD (QFT+, Low):", PowerToSample(df, 0.9, "onlypos", "low", "Leaky.POD"), "\n")
cat("Leaky POD/POI (QFT+, High):", PowerToSample(df, 0.9, "onlypos", "high", "Leaky.PODPOI"), "\n")
