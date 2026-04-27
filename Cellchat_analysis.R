# ============================================================================
# CELLCHAT CELL-CELL COMMUNICATION ANALYSIS
# Dispersed vs Clustered M+2+6- Malignant B Cells — DLBCL
#
# Key methodological note:
#   Communication probabilities aggregated using truncatedMean (trim = 0)
#   instead of the CellChat default trimean, to reduce outlier influence
#   while retaining sensitivity to lowly expressed ligands.
#
# Input:
#   cells.use — Seurat object subsetted to cells of interest
#               with cell type identity in metadata
#
# Output:
#   - CellChat object with computed communication probabilities
#   - Circle plots of interaction counts and weights
#   - Pathway-level communication visualisations
# ============================================================================

library(CellChat)

# ============================================================================
# 1. CREATE CELLCHAT OBJECT
# ============================================================================
cellchat <- createCellChat(
  object = cells.use,
  meta   = cells.use@meta.data
)

# ============================================================================
# 2. SET LIGAND-RECEPTOR DATABASE
# ============================================================================
# Use full human CellChat database (all signalling categories)
# To restrict to secreted signalling only, uncomment the subsetDB line below

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

# ============================================================================
# 3. PREPROCESS EXPRESSION DATA
# ============================================================================
# Subset to genes in the database to reduce computation
cellchat <- subsetData(cellchat)

# Parallelise using 4 workers
future::plan("multisession", workers = 4)

# Identify overexpressed genes and ligand-receptor interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression onto protein-protein interaction network
# to smooth sparse single-cell data
cellchat <- projectData(cellchat, PPI.human)

# ============================================================================
# 4. COMPUTE COMMUNICATION PROBABILITIES
# ============================================================================
# type = "truncatedMean": uses truncated mean instead of default trimean
#   trim = 0: no trimming of extreme values (full truncatedMean)
#   population.size = FALSE: do not correct for cell population size
#   k.min = 2: minimum number of cells per group
#   nboot = 100: permutation tests for significance
#   raw.use = TRUE: use raw count data

cellchat <- computeCommunProb(
  cellchat,
  type            = "truncatedMean",
  trim            = 0,
  population.size = FALSE,
  k.min           = 2,
  nboot           = 100,
  raw.use         = TRUE
)

# Compute pathway-level communication probabilities
# thresh = 100: retain all pathways regardless of p-value threshold
cellchat <- computeCommunProbPathway(cellchat, thresh = 100)

# Aggregate cell-cell interaction network
cellchat <- aggregateNet(cellchat)

# ============================================================================
# 5. VISUALISE INTERACTION NETWORK
# ============================================================================
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1, 2), xpd = TRUE)

# Number of interactions between cell types
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale  = TRUE,
  label.edge    = TRUE,
  title.name    = "Number of interactions"
)

# Interaction strength between cell types
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale  = TRUE,
  label.edge    = TRUE,
  title.name    = "Interaction weights/strength"
)

# ============================================================================
# 6. PER CELL TYPE OUTGOING INTERACTION STRENGTH
# ============================================================================
mat <- cellchat@net$weight

par(mfrow = c(1, 2), xpd = TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,
                 nrow     = nrow(mat),
                 ncol     = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight  = groupSize,
    weight.scale   = TRUE,
    edge.weight.max = max(mat),
    title.name     = rownames(mat)[i]
  )
}

# ============================================================================
# 7. PATHWAY-SPECIFIC VISUALISATION
# ============================================================================
# Example: visualise prostaglandin signalling pathway
# Replace "Prostaglandin" with any pathway of interest

pathways.show <- c("Prostaglandin")

# Hierarchy plot
vertex.receiver <- seq(1, 4)
netVisual_aggregate(
  cellchat,
  signaling       = pathways.show,
  vertex.receiver = vertex.receiver,
  idents.use      = cellchat@idents
)

# Circle plot
par(mfrow = c(1, 1))
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  layout    = "circle"
)
