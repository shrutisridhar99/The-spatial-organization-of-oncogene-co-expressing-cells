# ============================================================================
# DIGITAL SPATIAL PROFILING (DSP) — DIFFERENTIAL EXPRESSION ANALYSIS
# CD3+ T Cell Compartment: Dispersed vs Clustered DE Cases
#
# Workflow:
#   1. Load Q3-normalised GeoMx DSP data
#   2. Separate CD3+ and CD20+ ROI compartments
#   3. Score CD20+ compartment for dispersed-DE signature using GSVA
#   4. Stratify patients by signature score (Q1 vs Q4)
#   5. Run limma-voom DEG analysis in matched CD3+ compartment
#
# Input:
#   - Q3-normalised DSP expression matrix (Excel, sheet 3)
#   - Dispersed-DE upregulated gene list (CSV)
#
# Output:
#   - Differentially expressed genes in CD3+ compartment
#     between dispersed-DE (Q4) and clustered-DE (Q1) cases
# ============================================================================

library(dplyr)
library(limma)
library(GSVA)

# ============================================================================
# 1. LOAD Q3-NORMALISED DSP DATA
# ============================================================================
# Sheet 3 contains the Q3-normalised expression matrix
# Rows = genes (TargetName), Columns = ROIs

dds <- readxl::read_xlsx(
  "Q3_normalization_DLBCL.xlsx",
  sheet = 3
)
dds           <- as.data.frame(dds)
rownames(dds) <- dds$TargetName
dds           <- dds[, -1]   # Remove TargetName column after setting rownames

# ============================================================================
# 2. REMOVE DUPLICATE ROIs
# ============================================================================
# These CD68 ROIs have matched CD20 and CD3 counterparts that are retained.
# Removing CD68 ROIs to avoid duplication in downstream analysis.

cols_to_remove <- c(
  "DLBCL A3 | 009 | CD68", "DLBCL A3 | 007 | CD68",
  "DLBCL A3 | 013 | CD68", "DLBCL A3 | 010 | CD68",
  "DLBCL A3 | 011 | CD68", "DLBCL A3 | 004 | CD68",
  "DLBCL A3 | 003 | CD68", "DLBCL B2 | 001 | CD68",
  "DLBCL B2 | 007 | CD68", "DLBCL B2 | 008 | CD68",
  "DLBCL B2 | 005 | CD68", "DLBCL B2 | 004 | CD68",
  "DLBCL B2 | 006 | CD68", "DLBCL B2 | 009 | CD68",
  "DLBCL B2 | 010 | CD68", "DLBCL B2 | 013 | CD68",
  "DLBCL C1 | 001 | CD68", "DLBCL C1 | 002 | CD68",
  "DLBCL C1 | 014 | CD68", "DLBCL C1 | 003 | CD68",
  "DLBCL C1 | 016 | CD68", "DLBCL C1 | 005 | CD68",
  "DLBCL C1 | 008 | CD68", "DLBCL C1 | 010 | CD68"
)

# Derive matched CD20 and CD3 column names from CD68 list
new_names_cd20 <- sub("CD68", "CD20", cols_to_remove)
new_names_cd3  <- sub("CD68", "CD3",  cols_to_remove)

# ============================================================================
# 3. SEPARATE CD3+ AND CD20+ COMPARTMENTS
# ============================================================================
# Remove flagged duplicate ROIs and retain only the compartment of interest.
# CD3 column names are renamed to match CD20 for paired sample alignment.

# CD3+ compartment
dds_cd3_1     <- dds[, !(colnames(dds) %in% new_names_cd3)]
cd3_indices   <- grep("CD3", colnames(dds_cd3_1))
dds_cd3_final <- dds_cd3_1[, cd3_indices]

# CD20+ compartment
dds_cd20_1     <- dds[, !(colnames(dds) %in% new_names_cd20)]
cd20_indices   <- grep("CD20", colnames(dds_cd20_1))
dds_cd20_final <- dds_cd20_1[, cd20_indices]

# Rename CD3 columns to CD20 for alignment, then retain only shared samples
colnames(dds_cd3_final) <- gsub("CD3", "CD20", colnames(dds_cd3_final))
common_cols          <- intersect(colnames(dds_cd3_final), colnames(dds_cd20_final))
dds_cd3_final_1      <- dds_cd3_final[, common_cols]
dds_cd20_final_1     <- dds_cd20_final[, common_cols]

# ============================================================================
# 4. SCORE CD20+ COMPARTMENT FOR DISPERSED-DE SIGNATURE USING GSVA
# ============================================================================
# Load the dispersed-DE upregulated gene list derived from BCA bulk RNA-seq
# (274-gene signature, Q4 CI vs Q1 CI differential expression)

deg <- read.csv("results_final_Quartile_4_v_1.csv")

# Invert fold change direction: signature was derived as dispersed upregulated
# relative to clustered, so positive FC = dispersed-enriched
deg$log2FoldChange <- -1 * deg$log2FoldChange

# Retain upregulated genes (dispersed-enriched, padj < 0.05)
firma_up <- deg %>% filter(log2FoldChange > 0 & padj < 0.05)

# Score each CD20+ ROI for dispersed-DE signature enrichment using z-score GSVA
gbmPar <- zscoreParam(
  data.matrix(dds_cd20_final_1),
  list(firma_up$X)
)
gbm_es <- gsva(gbmPar)
gbm_es <- t(as.data.frame(gbm_es))

# ============================================================================
# 5. STRATIFY SAMPLES BY DISPERSED-DE SIGNATURE SCORE
# ============================================================================
# Assign quartiles based on GSVA enrichment score in CD20+ compartment
# Q4 = highest dispersed-DE signature = dispersed-like cases
# Q1 = lowest dispersed-DE signature  = clustered-like cases

gbm_es           <- as.data.frame(gbm_es)
gbm_es$quartile  <- ntile(gbm_es$V1, 4)

# Retain only Q1 and Q4 for comparison (exclude intermediate cases)
dds_m26_q4 <- gbm_es %>% filter(quartile %in% c(1, 4))

# Subset matched CD3+ compartment to Q1 and Q4 samples only
dds_cd3_m26 <- dds_cd3_final_1[, colnames(dds_cd3_final_1) %in% rownames(dds_m26_q4)]
dds_m26_q4  <- dds_m26_q4[colnames(dds_cd3_m26), , drop = FALSE]

# ============================================================================
# 6. DIFFERENTIAL EXPRESSION — limma-voom
# ============================================================================
# Compare CD3+ compartment gene expression between dispersed-DE (Q4)
# and clustered-DE (Q1) cases

# Log2 transform Q3-normalised counts
log_exprs <- log2(dds_cd3_m26 + 1)

# Design matrix: Q1 vs Q4 (no intercept)
dds_m26_q4$quartile <- as.factor(dds_m26_q4$quartile)
design <- model.matrix(~ 0 + dds_m26_q4$quartile)
colnames(design) <- gsub("dds_m26_q4\\$quartile", "quartile", colnames(design))

# Contrast: Q4 (dispersed-like) vs Q1 (clustered-like)
contrast.matrix <- makeContrasts(
  dispersed_vs_clustered = quartile4 - quartile1,
  levels = design
)

# voom transformation — estimates mean-variance relationship
v    <- voom(log_exprs, design)
fit  <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract all results with Benjamini-Hochberg FDR correction
res     <- topTable(fit2, number = Inf, adjust = "BH")
res_sig <- res %>% filter(adj.P.Val < 0.05)

cat("Significant DEGs (adj.P.Val < 0.05):", nrow(res_sig), "\n")

# Export results
write.csv(res,     "DSP_DEG_all.csv",        row.names = TRUE)
write.csv(res_sig, "DSP_DEG_significant.csv", row.names = TRUE)
