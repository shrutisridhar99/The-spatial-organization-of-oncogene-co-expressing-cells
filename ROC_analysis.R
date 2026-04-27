# ============================================================================
# ROC / AUC ANALYSIS — 2-YEAR OS AND PFS
# Comparing predictive performance of:
#   Model 1: M26 DE score + Cell of Origin (COO)
#   Model 2: M26 DE score + Mean Clustering Index (CI)
#
# Statistical comparison via DeLong test (pROC::roc.test)
#
# Input:
#   df_bca — BCA cohort dataframe with columns:
#     M26          : continuous DE score
#     COO_new      : cell of origin ("ABC", "GCB", or other)
#     mean         : mean clustering index (CI) per patient
#     OS, OS_status: overall survival time and event indicator
#     PFS, PFS_status: progression-free survival time and event indicator
#
# Output:
#   - Smoothed ROC curves for 2-year OS and PFS
#   - AUC values for each model
#   - DeLong p-value comparing model pairs
# ============================================================================

library(pROC)
library(ggplot2)
library(dplyr)
library(patchwork)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
# df_bca <- read.csv("BCA_with_clinical.csv")

# ============================================================================
# 2. DATA PREPARATION
# ============================================================================
# Define CI quartile thresholds for spatial group assignment
q25 <- quantile(df_bca$mean, 0.25, na.rm = TRUE)
q75 <- quantile(df_bca$mean, 0.75, na.rm = TRUE)

df <- df_bca %>%
  mutate(
    # Clean M26 DE score — coerce non-numeric to NA
    M26_bin    = ifelse(M26 %in% c("#N/A", "NA", ""),
                        NA, as.numeric(as.character(M26))),

    # Ensure survival variables are numeric
    OS_num     = as.numeric(as.character(OS)),
    PFS_num    = suppressWarnings(as.numeric(as.character(PFS))),
    OS_status  = as.numeric(as.character(OS_status)),
    PFS_status = as.numeric(as.character(PFS_status)),

    # ── 2-year binary outcomes ──────────────────────────────────────────────
    # Event   : died/progressed before 2 years        → 1
    # Success : survived/no event beyond 2 years      → 0
    # Censored before 2 years                         → NA (excluded)
    OS_2yr = case_when(
      OS_status == 1 & OS_num <= 2 ~ 1,
      OS_num > 2                   ~ 0,
      TRUE                         ~ NA_real_
    ),
    PFS_2yr = case_when(
      PFS_status == 1 & PFS_num <= 2 ~ 1,
      PFS_num > 2                    ~ 0,
      TRUE                           ~ NA_real_
    ),

    # ── COO binary encoding ─────────────────────────────────────────────────
    COO_clean = ifelse(COO_new %in% c("GCB", "ABC"), COO_new, NA),
    COO_bin   = ifelse(COO_clean == "ABC", 1,
                       ifelse(COO_clean == "GCB", 0, NA)),

    # ── Spatial group (Q1 = Dispersed, Q4 = Clustered) ─────────────────────
    Spatial_Group = case_when(
      mean <= q25 ~ "Dispersed",
      mean >= q75 ~ "Clustered",
      TRUE        ~ NA_character_
    )
  )

# Check 2-year outcome counts
cat("2yr OS  — events:", sum(df$OS_2yr  == 1, na.rm = TRUE),
    " non-events:", sum(df$OS_2yr  == 0, na.rm = TRUE),
    " excluded:", sum(is.na(df$OS_2yr)), "\n")
cat("2yr PFS — events:", sum(df$PFS_2yr == 1, na.rm = TRUE),
    " non-events:", sum(df$PFS_2yr == 0, na.rm = TRUE),
    " excluded:", sum(is.na(df$PFS_2yr)), "\n")

# ============================================================================
# 3. RESTRICT TO COMPLETE CASES FOR ROC ANALYSIS
# ============================================================================
df_both <- df %>%
  filter(!is.na(M26_bin), !is.na(COO_bin), !is.na(mean),
         !is.na(OS_2yr),  !is.na(PFS_2yr))

cat("N for ROC analysis:", nrow(df_both), "\n")

# ============================================================================
# 4. FIT LOGISTIC REGRESSION MODELS
# ============================================================================
# Model 1: M26 DE score + COO (clinical surrogate for spatial topology)
# Model 2: M26 DE score + mean CI (direct spatial measure)

# 2-year OS
glm_coo_os  <- glm(OS_2yr ~ M26_bin + COO_bin, data = df_both, family = binomial)
glm_mean_os <- glm(OS_2yr ~ M26_bin + mean,    data = df_both, family = binomial)

# 2-year PFS
glm_coo_pfs  <- glm(PFS_2yr ~ M26_bin + COO_bin, data = df_both, family = binomial)
glm_mean_pfs <- glm(PFS_2yr ~ M26_bin + mean,    data = df_both, family = binomial)

# ============================================================================
# 5. COMPUTE ROC CURVES
# ============================================================================
roc_coo_os  <- roc(df_both$OS_2yr,  predict(glm_coo_os,  type = "response"), quiet = TRUE)
roc_mean_os <- roc(df_both$OS_2yr,  predict(glm_mean_os, type = "response"), quiet = TRUE)

roc_coo_pfs  <- roc(df_both$PFS_2yr, predict(glm_coo_pfs,  type = "response"), quiet = TRUE)
roc_mean_pfs <- roc(df_both$PFS_2yr, predict(glm_mean_pfs, type = "response"), quiet = TRUE)

# ============================================================================
# 6. PLOT FUNCTION — SMOOTHED ROC WITH DELONG TEST
# ============================================================================
col_coo  <- "#C1666B"   # Model 1: M26 + COO
col_mean <- "#4A7C99"   # Model 2: M26 + CI

make_roc_plot <- function(roc_coo, roc_mean, outcome_label) {

  # Smooth ROC curves using density-based method
  roc_coo_sm  <- smooth(roc_coo,  method = "density")
  roc_mean_sm <- smooth(roc_mean, method = "density")

  # DeLong test comparing AUCs between the two models
  test <- roc.test(roc_coo, roc_mean)

  # Build plotting dataframe
  roc_df <- bind_rows(
    data.frame(
      FPR   = 1 - roc_coo_sm$specificities,
      TPR   = roc_coo_sm$sensitivities,
      Model = paste0("M26 + COO  (AUC = ", round(auc(roc_coo),  3), ")")
    ),
    data.frame(
      FPR   = 1 - roc_mean_sm$specificities,
      TPR   = roc_mean_sm$sensitivities,
      Model = paste0("M26 + Mean CI (AUC = ", round(auc(roc_mean), 3), ")")
    )
  )
  roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))

  ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
    geom_abline(linetype = "dashed", color = "#CCCCCC", linewidth = 0.8) +
    geom_line(linewidth = 1.3) +
    annotate("text", x = 0.97, y = 0.03,
             label  = paste0("DeLong p = ", round(test$p.value, 3)),
             size   = 3.2, hjust = 1, color = "grey40") +
    scale_color_manual(values = c(col_coo, col_mean)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(
      title = paste0("2-year ", outcome_label, " ROC"),
      x     = "1 - Specificity",
      y     = "Sensitivity",
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position   = c(0.65, 0.12),
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.text       = element_text(size = 8.5),
      axis.line         = element_line(color = "grey40"),
      panel.grid.major  = element_line(color = "grey94", linewidth = 0.4)
    )
}

# ============================================================================
# 7. GENERATE AND COMBINE PLOTS
# ============================================================================
p_os  <- make_roc_plot(roc_coo_os,  roc_mean_os,  "OS")
p_pfs <- make_roc_plot(roc_coo_pfs, roc_mean_pfs, "PFS")

p_os | p_pfs
