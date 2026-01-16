############################################################
## 03_clustering_banksy.R
## Clustering + CLUSTREE + BANKSY on Visium HD SUBSET
############################################################

set.seed(1234)


suppressPackageStartupMessages({
  library(hdf5r)
  library(Seurat)
  library(clustree)
  library(Banksy)
  library(SeuratWrappers)
  library(spacexr)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
})

options(future.globals.maxSize = 1e10)

## -------------------- PATHS --------------------

plot_dir <- file.path(Sys.getenv("HOME"), "results/visium/figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(plot, filename, width = 12, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE
  )
}

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

## -------------------- LOAD SUBSET OBJECT --------------------

message("Loading subset object")
bc_visium_subset <- readRDS(file.path(Sys.getenv("HOME"), "results/visium/objects/bc_visium_subset.rds"))
DefaultAssay(bc_visium_subset) <- "Spatial"

## -------------------- STANDARD PREPROCESSING --------------------

message("Normalize / HVG / Scale / PCA")

bc_visium_subset <- NormalizeData(bc_visium_subset)
bc_visium_subset <- FindVariableFeatures(bc_visium_subset, nfeatures = 2000)
bc_visium_subset <- ScaleData(bc_visium_subset, features = VariableFeatures(bc_visium_subset), verbose = TRUE)
bc_visium_subset <- RunPCA(bc_visium_subset, npcs = 30)

## -------------------- PCA ELBOW --------------------

p_elbow <- ElbowPlot(bc_visium_subset, ndims = 20) + plot_theme()
save_png(p_elbow, "07-SUBSET_PCA_Elbow.png", 6.5, 5.5)
rm(p_elbow); gc()

## -------------------- CLUSTERING (PRE-BANKSY) --------------------

message("Pre-BANKSY clustering")

bc_visium_subset <- FindNeighbors(bc_visium_subset, dims = 1:30)

bc_visium_subset <- FindClusters(
  bc_visium_subset,
  resolution = seq(0.1, 2, by = 0.1),
  algorithm  = 4,
  random.seed = 1234
)

## -------------------- CLUSTREE (PRE-BANKSY) --------------------

p_clustree_pre <- clustree(
  bc_visium_subset,
  prefix = "Spatial_snn_res."
) + plot_theme()


save_png(p_clustree_pre, "08-SUBSET_Clustree_PreBANKSY.png", 10, 12)
rm(p_clustree_pre); gc()

## -------------------- UMAP (PRE-BANKSY) --------------------

bc_visium_subset <- RunUMAP(
  bc_visium_subset,
  dims = 1:20,
  reduction.name = "umap.prebanksy"
)

p_umap_pre <- DimPlot(
  bc_visium_subset,
  reduction = "umap.prebanksy",
  group.by = "Spatial_snn_res.0.8",   # adjust resolution
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = FALSE
) + ggtitle("Pre-BANKSY clustering") + plot_theme()

save_png(p_umap_pre, "09-SUBSET_UMAP_PreBANKSY.png", 7, 6)
rm(p_umap_pre); gc()

## -------------------- HARD MEMORY CLEANUP BEFORE BANKSY ---------------
bc_visium_subset@graphs <- list()
bc_visium_subset@reductions <- bc_visium_subset@reductions["pca"]
gc()

## -------------------- BANKSY --------------------

DefaultAssay(bc_visium_subset) <- "Spatial"

bc_visium_subset <- RunBanksy(
  bc_visium_subset,
  lambda   = 0.6,
  verbose  = TRUE,
  assay    = "Spatial",
  slot     = "data",
  features = "variable",
  k_geom   = 15
)

DefaultAssay(bc_visium_subset) <- "BANKSY"

bc_visium_subset <- RunPCA(
  bc_visium_subset,
  assay = "BANKSY",
  reduction.name = "pca.banksy",
  features = rownames(bc_visium_subset),
  npcs = 30
)

p_banksy_elbow <- ElbowPlot(
  bc_visium_subset,
  reduction = "pca.banksy",
  ndims = 15
) + plot_theme()

save_png(p_banksy_elbow, "10-SUBSET_Banksy_Elbow.png", 6.5, 5.5)
rm(p_banksy_elbow); gc()

## -------------------- CLUSTERING (POST-BANKSY) --------------------

message("Post-BANKSY clustering")

bc_visium_subset <- FindNeighbors(
  bc_visium_subset,
  reduction  = "pca.banksy",
  dims       = 1:20,
  graph.name = c("banksy_nn", "banksy_snn")
)

print(names(bc_visium_subset@graphs))

res_vec <- seq(0.1, 1.1, by = 0.2)

for (res in res_vec) {
  message("Clustering at resolution = ", res)

  bc_visium_subset <- FindClusters(
    bc_visium_subset,
    graph.name  = "banksy_snn",
    resolution  = res,
    algorithm   = 4,
    random.seed = 1234,
    verbose     = FALSE
  )

  bc_visium_subset[[paste0("banksy_cluster_", res)]] <- Idents(bc_visium_subset)
}

## -------------------- CLUSTREE (POST-BANKSY) --------------------

p_clustree <- clustree(bc_visium_subset, prefix = "banksy_cluster_") + plot_theme()
save_png(p_clustree, "11-SUBSET_Clustree_PostBANKSY.png", 10, 12)



## -------------------- FINAL BANKSY UMAP --------------------

RES_USE <- "banksy_cluster_0.3"
stopifnot(RES_USE %in% colnames(bc_visium_subset@meta.data))

message("Running UMAP for selected resolution: ", RES_USE)

bc_visium_subset <- RunUMAP(
  bc_visium_subset,
  reduction = "pca.banksy",
  reduction.name = "umap.banksy",
  dims = 1:30
)

Idents(bc_visium_subset) <- RES_USE

p_umap <- DimPlot(
  bc_visium_subset,
  reduction = "umap.banksy",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = FALSE
) + ggtitle(paste0("BANKSY (", RES_USE, ")")) + plot_theme()

p_spatial <- SpatialDimPlot(
  bc_visium_subset,
  group.by = RES_USE,
  pt.size.factor = 6,
  label = TRUE,
  repel = TRUE,
  label.size = 4
) + plot_theme()

combined <- p_umap | p_spatial
save_png(combined, paste0("12-SUBSET_Banksy_UMAP_Spatial_", RES_USE, ".png"), 13, 6)
rm(p_umap, p_spatial, combined); gc()

## -------------------- SAVE OUTPUT --------------------
message("Saving object with multi-resolution clusters + UMAP")
saveRDS(
  bc_visium_subset,
  file = file.path(Sys.getenv("HOME"), "results/visium/objects/bc_visium_banksy_subset.rds")
)

message("Done")

