# ============================================================================
# TMA DUPLICATE CORE REPRODUCIBILITY ANALYSIS
# BCA Cohort — Clustering Index (CI) Agreement Between Duplicate Cores
#
# Outputs:
#   - Intraclass correlation coefficient (ICC)
#   - Spearman rank correlation
#   - Spearman-Brown corrected ICC for averaged cores
#   - Cohen's Kappa (median-split binary agreement)
#   - Chi-square test of concordance
#   - Heatmap of CI per core per patient
#   - Scatter plot of Core 1 vs Core 2 CI
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(irr)        # ICC and Kappa
library(patchwork)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
# bca_raw should contain columns: Name, Mean (CI per core), Patient
# Load your Geyer output CSV for the BCA cohort

bca_raw <- read.csv("BCA_Mp2p6n_Geyer_modified_bestAIC.csv")

# ============================================================================
# 2. PREPARE HEATMAP DATA
# ============================================================================
# Keep only patients with more than one core
# Order cores within each patient and compute mean CI per patient

bca_heatmap <- bca_raw %>%
  mutate(Mean = as.numeric(Mean)) %>%
  filter(!is.na(Mean), !is.na(Patient)) %>%
  group_by(Patient) %>%
  mutate(
    n_cores  = n(),
    mean_CI  = mean(Mean, na.rm = TRUE),
    core_num = paste0("Core ", row_number())
  ) %>%
  ungroup() %>%
  filter(n_cores > 1)   # Retain only patients with duplicate cores

# ============================================================================
# 3. ICC — INTRACLASS CORRELATION COEFFICIENT
# ============================================================================
# Pivot to wide format: one row per patient, one column per core
# Only patients with complete data across both cores are retained

bca_icc_data <- bca_heatmap %>%
  select(Patient, core_num, Mean) %>%
  pivot_wider(names_from = core_num, values_from = Mean) %>%
  select(-Patient) %>%
  filter(complete.cases(.))

# Two-way mixed ICC, agreement type
icc_result <- icc(bca_icc_data, model = "twoway", type = "agreement")
print(icc_result)

# Spearman-Brown prophecy formula:
# Estimates the reliability gain from averaging two cores
# ICC_averaged = (2 * ICC_single) / (1 + ICC_single)
icc_single  <- icc_result$value
icc_average <- (2 * icc_single) / (1 + icc_single)

cat("ICC single core:              ", round(icc_single,  3), "\n")
cat("ICC averaged (Spearman-Brown):", round(icc_average, 3), "\n")

# ============================================================================
# 4. SCATTER PLOT — CORE 1 vs CORE 2 CI
# ============================================================================
bca_icc_plot <- bca_heatmap %>%
  select(Patient, core_num, Mean) %>%
  pivot_wider(names_from = core_num, values_from = Mean) %>%
  filter(complete.cases(.))

ax_min <- floor(min(bca_icc_plot$`Core 1`, bca_icc_plot$`Core 2`, na.rm = TRUE))
ax_max <- ceiling(max(bca_icc_plot$`Core 1`, bca_icc_plot$`Core 2`, na.rm = TRUE))

p_scatter <- ggplot(bca_icc_plot, aes(x = `Core 1`, y = `Core 2`)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  geom_point(colour = "maroon", alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
  stat_cor(method = "spearman", size = 4,
           label.x.npc = "left", label.y.npc = "top") +
  annotate("text", x = ax_max, y = ax_min,
           label = paste0("ICC = ", round(icc_single, 2)),
           hjust = 1, vjust = -1, size = 4, colour = "grey30") +
  scale_x_continuous(limits = c(ax_min, ax_max)) +
  scale_y_continuous(limits = c(ax_min, ax_max)) +
  coord_fixed(ratio = 1) +
  labs(x = "Core 1 CI", y = "Core 2 CI",
       title = "Core 1 vs Core 2 CI agreement") +
  theme_classic(base_size = 12)

print(p_scatter)

# ============================================================================
# 5. HEATMAP — CI PER CORE PER PATIENT
# ============================================================================
# Patients ordered by mean CI across cores

p_heatmap <- ggplot(bca_heatmap,
                    aes(x    = core_num,
                        y    = reorder(Patient, mean_CI),
                        fill = Mean)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "maroon",
    mid      = "white",
    high     = "navy",
    midpoint = median(bca_heatmap$Mean, na.rm = TRUE),
    name     = "CI"
  ) +
  labs(x     = "Core",
       y     = "Patient (ordered by mean CI)",
       title = "CI per core per patient — BCA cohort") +
  theme_classic(base_size = 12) +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(size = 11))

print(p_heatmap)

# ============================================================================
# 6. BINARY CONCORDANCE — MEDIAN SPLIT
# ============================================================================
# Dichotomise each core at the median CI and assess agreement
# between Core 1 and Core 2 group assignments

bca_binary <- bca_icc_plot %>%
  mutate(
    group_core1 = ifelse(`Core 1` <= median(`Core 1`, na.rm = TRUE),
                         "Low CI", "High CI"),
    group_core2 = ifelse(`Core 2` <= median(`Core 2`, na.rm = TRUE),
                         "Low CI", "High CI"),
    group_core1 = factor(group_core1, levels = c("Low CI", "High CI")),
    group_core2 = factor(group_core2, levels = c("Low CI", "High CI"))
  )

# Cohen's Kappa — agreement beyond chance
kappa_result <- kappa2(
  bca_binary %>%
    select(group_core1, group_core2) %>%
    mutate(across(everything(), as.character))
)
print(kappa_result)

# Chi-square test of concordance
confusion_mat <- table(bca_binary$group_core1, bca_binary$group_core2)
print(chisq.test(confusion_mat))

# ============================================================================
# 7. CORRELATION OF EACH CORE WITH MEAN CI
# ============================================================================
# Validates that averaging cores improves reliability

bca_icc_plot2 <- bca_icc_plot %>%
  mutate(mean_CI = (`Core 1` + `Core 2`) / 2)

cat("Spearman correlation — Core 1 vs mean CI:",
    round(cor(bca_icc_plot2$`Core 1`, bca_icc_plot2$mean_CI,
              method = "spearman"), 3), "\n")
cat("Spearman correlation — Core 2 vs mean CI:",
    round(cor(bca_icc_plot2$`Core 2`, bca_icc_plot2$mean_CI,
              method = "spearman"), 3), "\n")

# ============================================================================
# 8. COMBINED FIGURE
# ============================================================================
p_heatmap | p_scatter
