print("=== CellChat module: start ===")
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(future)
  library(jsonlite)
})

options(future.globals.maxSize = 1e10)

## =========================================================
## INPUT / OUTPUT
## =========================================================
in_rds <- file.path(Sys.getenv("HOME"), "results/visium/objects/visium_crc5_query_annotated_rctd_COARSEref_QUERYsubsample.rds")

# Save ALL CellChat plots in the SAME folder used by previous scripts
## -------------------- GLOBAL PLOT SETTINGS --------------------
## Keep fonts and sizes consistent across all images
## -------------------- Output folder for plots --------------------
plot_dir <- file.path(Sys.getenv("HOME"), "results/visium/figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 14))

plot_theme <- function() {
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )
}

save_png <- function(plot, filename, width = 14.5, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    device = "png",
    dpi = dpi,
    limitsize = FALSE
  )
}

# Where to save CellChat RDS + tables
outdir <- file.path(Sys.getenv("HOME"), "results/visium/figures/cellchat_crc5_rctd_coarse")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## LOAD SEURAT OBJECT
## =========================================================
print("Loading Seurat object...")
obj <- readRDS(in_rds)

stopifnot("full_first_type" %in% colnames(obj@meta.data))
obj$full_first_type <- as.character(obj$full_first_type)
obj$full_first_type[is.na(obj$full_first_type)] <- "Unknown"
Idents(obj) <- "full_first_type"

print("Cell counts per group (full_first_type):")
print(table(Idents(obj)))

## =========================================================
## ASSAY SETUP (Visium HD bins)
## CellChat expects normalized data -> ensure Spatial.008um has 'data'
## =========================================================
DefaultAssay(obj) <- "Spatial.008um"
Idents(obj) <- "full_first_type"

## ---------------------------------------------------------
## Ensure normalized data layer exists
## ---------------------------------------------------------
# Seurat v5 layer check (kept simple and robust)

has_data_layer <- FALSE
try({
  has_data_layer <- ("data" %in% Layers(obj[["Spatial.008um"]]))
}, silent = TRUE)

if (!has_data_layer) {
  print("No 'data' layer found in Spatial.008um -> running NormalizeData(LogNormalize)...")
  obj <- NormalizeData(
    obj,
    assay = "Spatial.008um",
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    verbose = FALSE
  )
}

data.input <- GetAssayData(obj, assay = "Spatial.008um", slot = "data")
meta <- obj@meta.data %>% dplyr::select(full_first_type)

## =========================================================
## SPATIAL COORDINATES (Visium HD)
## =========================================================
spatial_locs <- GetTissueCoordinates(
  obj,
  cols = c("x", "y")
)
spatial_locs <- spatial_locs[, 1:2]

## =========================================================
## VISIUM HD SCALING FACTORS (8 μm bins)
## =========================================================

scale_factors <- fromJSON(
  txt = file.path(Sys.getenv("HOME"), "data/visium/GEO_CRC5_VHD/outs/binned_outputs/square_008um/spatial/scalefactors_json.json")
)

# Approximate full-res pixel diameter for 8 μm bin
spot_size <- 64
conversion_factor <- spot_size / scale_factors$spot_diameter_fullres

spatial_factors <- data.frame(
  ratio = conversion_factor,
  tol   = spot_size / 2
)

## =========================================================
## COMPUTE SPATIAL DISTANCES
## =========================================================
library(future)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)

spatial <- computeCellDistance(
  coordinates = spatial_locs,
  ratio = spatial_factors$ratio,
  tol   = spatial_factors$tol
)

## =========================================================
## CREATE CELLCHAT (HUMAN DATABASE)
## =========================================================
print("Creating CellChat object...")
cellchat <- createCellChat(
  object          = data.input,
  meta            = meta,
  group.by        = "full_first_type",
  datatype        = "spatial",
  coordinates     = spatial_locs,
  spatial.factors = spatial_factors
)

print("Setting CellChat database: HUMAN")
cellchat@DB <- CellChatDB.human

## =========================================================
## PREPROCESS + PARALLEL
## =========================================================
print("Subsetting signaling genes...")
cellchat <- subsetData(cellchat)

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
print(paste0("Using cores (SLURM_CPUS_PER_TASK): ", ncores))
plan(multisession, workers = ncores)

print("Identifying overexpressed genes/interactions...")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## =========================================================
## COMPUTE COMMUNICATION
## =========================================================
print("Computing communication probabilities...")
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, contact.range = 50)

# For Visium bins, filter more strictly to reduce noisy interactions
print("Filtering communications (min.cells = 5)...")
cellchat <- filterCommunication(cellchat, min.cells = 5)

print("Computing pathway probabilities + aggregating network...")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

## =========================================================
## SAVE CHECKPOINTS
## =========================================================
saveRDS(cellchat, file = file.path(outdir, "cellchat_crc5_rctd_coarse.rds"))
print(paste0("Saved: ", file.path(outdir, "cellchat_crc5_rctd_coarse.rds")))

print("Exporting interaction tables...")
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = file.path(outdir, "cellchat_interactions_all.csv"), row.names = FALSE)

df.path <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.path, file = file.path(outdir, "cellchat_pathways_all.csv"), row.names = FALSE)

## =========================================================
## PLOTS (Saved in SAME plot_dir as previous scripts)
## =========================================================

# 1) Circle: number of interactions
print("Plotting: circle (number of interactions)...")
p1 <- netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
# netVisual_circle returns a base plot; still works with ggsave if captured as recordedplot
# Safer: open device manually
png(file.path(plot_dir, "17-CellChat_Circle_Number.png"), width = 2000, height = 1800, res = 200)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
dev.off()

# 2) Circle: interaction strength
print("Plotting: circle (interaction strength)...")
png(file.path(plot_dir, "18-CellChat_Circle_Weight.png"), width = 2000, height = 1800, res = 200)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength"
)
dev.off()

# 3) Heatmap: interaction strength
print("Plotting: heatmap (interaction strength)...")
png(file.path(plot_dir, "19-CellChat_Heatmap_Weight.png"), width = 2200, height = 1800, res = 200)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# 4) RankNet: senders/receivers

# ---- Compute centrality (REQUIRED) ----
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

print("Plotting: rankNet...")

# ---- Outgoing signals ----
png(file.path(plot_dir, "20-CellChat_SignalingRole_Outgoing.png"),
    width = 12, height = 12, units = "in", res = 300)

par(cex = 1.6, mar = c(10, 10, 4, 4))

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()

# ---- Incoming signals ----
png(file.path(plot_dir, "21-CellChat_SignalingRole_Incoming.png"),
    width = 3600, height = 3600, res = 300)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()


## =========================================================
## CLEANUP
## =========================================================
plan(sequential)

print("=== CellChat module: finished successfully ===")
print(paste0("Plots saved to: ", plot_dir))
print(paste0("CellChat outputs saved to: ", outdir))
