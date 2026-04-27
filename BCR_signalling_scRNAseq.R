# ============================================================================
# BCR SIGNALLING & CD79B ANALYSIS IN scRNA-seq DATA
# Dispersed vs Clustered M+2+6- Malignant B Cells — DLBCL
#
# Inputs:
#   - Seurat object (dds) with metadata columns:
#       final         : broad cell type (e.g. "Malignant B cells")
#       cell_type     : fine cell type including spatial labels
#                       ("Dispersed M+2+6- cell", "Clustered M+2+6- cell")
#       RNA.subtype   : COO classification ("ABC", "GCB")
#
# Outputs:
#   - Violin plots: Chronic BCR score and CD79B expression
#       (i)  M+2+6- vs other malignant B cells
#       (ii) Dispersed vs Clustered M+2+6- cells
#   - DotPlot: individual BCR pathway genes
#   - UMAP: spatial phenotype highlighting
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ============================================================================
# 1. LOAD SEURAT OBJECT
# ============================================================================
# dds <- readRDS("path/to/your/seurat_object.rds")
# Ensure metadata contains: final, cell_type, RNA.subtype

# ============================================================================
# 2. SUBSET TO DLBCL MALIGNANT B CELLS
# ============================================================================
# Subset to malignant B cells across all subtypes
dds_malig_all <- subset(dds, final == "Malignant B cells")

# Assign group labels:
#   - M+2+6- cells identified by spatial cell_type label
#     (RNA.subtype is NA for these cells — do NOT use RNA.subtype to identify them)
#   - All other ABC/GCB malignant B cells = "Other malignant B cells"
dds_malig_all@meta.data <- dds_malig_all@meta.data %>%
  mutate(group = case_when(
    cell_type == "Dispersed M+2+6- cell" ~ "M+2+6-",
    cell_type == "Clustered M+2+6- cell" ~ "M+2+6-",
    RNA.subtype %in% c("ABC", "GCB")    ~ "Other malignant B cells",
    TRUE                                 ~ NA_character_
  ))

# Retain only labelled cells (M+2+6- and other DLBCL malignant B cells)
dds_dlbcl_malig <- subset(
  dds_malig_all,
  cells = rownames(dds_malig_all@meta.data)[
    !is.na(dds_malig_all@meta.data$group)
  ]
)

# ============================================================================
# 3. CHRONIC BCR SIGNALLING GENE SET
# ============================================================================
# Canonical BCR pathway components covering receptor proximal signalling,
# PI3K/AKT activation, NF-kB signalling, and downstream transcription factors

chronic_bcr_genes <- list(Chronic_BCR = c(
  "CD79A", "CD79B", "CD19",
  "PIK3CD", "PIK3R1", "AKT1", "AKT2",
  "SYK", "LYN", "BLK",
  "NFKB1", "RELA", "BCL10", "CARD11",
  "IRF4", "MYC", "BCL2"
))

# Compute module score across all DLBCL malignant B cells
dds_dlbcl_malig <- AddModuleScore(
  dds_dlbcl_malig,
  features = chronic_bcr_genes,
  name     = "ChronicBCR"
)

# Fetch CD79B expression per cell
dds_dlbcl_malig@meta.data$CD79B_expr <- FetchData(
  dds_dlbcl_malig,
  vars = "CD79B"
)[, 1]

# Prepare metadata for plotting
meta_dlbcl <- dds_dlbcl_malig@meta.data %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group,
                        levels = c("M+2+6-", "Other malignant B cells")))

# ============================================================================
# 4. VIOLIN PLOTS — M+2+6- vs OTHER MALIGNANT B CELLS
# ============================================================================

# Plot 1: Chronic BCR score
p1 <- ggplot(meta_dlbcl,
             aes(x = group, y = ChronicBCR1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format", size = 4.5) +
  scale_fill_manual(values = c("M+2+6-"                  = "maroon",
                               "Other malignant B cells" = "grey70")) +
  labs(x        = NULL,
       y        = "Chronic BCR Score",
       title    = "Chronic BCR signalling",
       subtitle = "M+2+6- vs other malignant B cells (DLBCL)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Plot 2: CD79B expression
p2 <- ggplot(meta_dlbcl,
             aes(x = group, y = CD79B_expr, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format", size = 4.5) +
  scale_fill_manual(values = c("M+2+6-"                  = "maroon",
                               "Other malignant B cells" = "grey70")) +
  labs(x        = NULL,
       y        = "CD79B expression",
       title    = "CD79B expression",
       subtitle = "M+2+6- vs other malignant B cells (DLBCL)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# ============================================================================
# 5. DISPERSED vs CLUSTERED M+2+6- CELLS
# ============================================================================
# Subset to M+2+6- cells only and recompute BCR score on this subset
# (recomputing on the subset ensures the background gene distribution
#  is specific to M+2+6- cells rather than all malignant B cells)

dds_m26_dlbcl <- subset(
  dds_dlbcl_malig,
  cell_type %in% c("Dispersed M+2+6- cell", "Clustered M+2+6- cell")
)

dds_m26_dlbcl <- AddModuleScore(
  dds_m26_dlbcl,
  features = chronic_bcr_genes,
  name     = "ChronicBCR_m26"
)

meta_spatial <- dds_m26_dlbcl@meta.data %>%
  mutate(spatial_group = factor(
    ifelse(cell_type == "Dispersed M+2+6- cell",
           "Dispersed (Q1)", "Clustered (Q4)"),
    levels = c("Dispersed (Q1)", "Clustered (Q4)")
  ))

# Plot 3: Chronic BCR score — Dispersed vs Clustered
p3 <- ggplot(meta_spatial,
             aes(x = spatial_group, y = ChronicBCR_m261, fill = spatial_group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format", size = 4.5) +
  scale_fill_manual(values = c("Dispersed (Q1)" = "maroon",
                               "Clustered (Q4)" = "navy")) +
  labs(x        = NULL,
       y        = "Chronic BCR Score",
       title    = "Chronic BCR signalling",
       subtitle = "Dispersed vs Clustered M+2+6- (DLBCL)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Combined violin figure
p1 | p2 | p3

# ============================================================================
# 6. DOTPLOT — INDIVIDUAL BCR PATHWAY GENES
# ============================================================================
# Shows mean expression and percent expressed for key BCR genes
# across Dispersed vs Clustered M+2+6- cells

bcr_genes_show <- c(
  "CD79A", "CD79B", "SYK", "BTK", "LYN",
  "BLNK", "PLCG2", "BCL10", "CARD11", "NFKB1"
)

DotPlot(
  dds_m26_dlbcl,
  features = rev(bcr_genes_show),
  group.by = "cell_type",
  cols     = c("lightgrey", "maroon"),
  dot.scale = 6
) +
  coord_flip() +
  labs(title = "BCR signalling genes: Dispersed vs Clustered M+2+6-",
       x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x    = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# ============================================================================
# 7. UMAP — HIGHLIGHT SPATIAL PHENOTYPES
# ============================================================================
# Overlay Dispersed and Clustered M+2+6- cells on the full UMAP
# Other cells shown in grey to provide context

dds@meta.data$highlight <- case_when(
  dds@meta.data$cell_type == "Dispersed M+2+6- cell" ~ "Dispersed (Q1)",
  dds@meta.data$cell_type == "Clustered M+2+6- cell" ~ "Clustered (Q4)",
  TRUE                                                ~ "Other"
)

DimPlot(
  dds,
  group.by = "highlight",
  cols     = c("Dispersed (Q1)" = "maroon",
               "Clustered (Q4)" = "navy",
               "Other"          = "grey85"),
  pt.size  = 0.3,
  order    = c("Dispersed (Q1)", "Clustered (Q4)")
) +
  theme_classic(base_size = 13) +
  labs(title = "M+2+6- spatial phenotypes on UMAP")
