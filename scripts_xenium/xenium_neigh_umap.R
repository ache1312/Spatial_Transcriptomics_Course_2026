#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Neighbors, Clustering, and UMAP
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)
library(clustree)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/neigh_umap_plots"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("QC plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc_banksy <- readRDS(file = "/home/courses/student75/xen_bc_banksy.rds")

DefaultAssay(xen_bc_banksy) = "Xenium"

# -------------------- Find Neighbors ---------------------

xen_bc_banksy <- FindNeighbors(xen_bc_banksy, 
                               dims = 1:20)

# -------------------- Find Clusters and clustree plot ----

xen_bc_banksy <- FindClusters(xen_bc_banksy, 
                              resolution = seq(0.1, 2, by = 0.1),
                              algorithm = 4,
                              random.seed = 1234)

cluster_tree_plot <- clustree(xen_bc_banksy, 
                              prefix = "Xenium_snn_res.")

# -------------------- UMAP -------------------------------

xen_bc_banksy <- RunUMAP(xen_bc_banksy, dims = 1:20)

# -------------------- Dim Plots --------------------------

umap_clusters <- DimPlot(xen_bc_banksy, group.by = "Xenium_snn_res.0.7", label = TRUE)

image_clusters <- ImageDimPlot(xen_bc_banksy, group.by = "Xenium_snn_res.0.7", dark.background = F,
                              size = 0.5)

# -------------------- Banksy Neighbors -------------------

DefaultAssay(xen_bc_banksy) = "BANKSY"

xen_bc_banksy <- FindNeighbors(xen_bc_banksy, 
                               dims = 1:20, reduction = "pca.banksy")

# -------------------- Banksy Clusters and clustree plot --

xen_bc_banksy <- FindClusters(xen_bc_banksy, 
                              resolution = seq(0.1, 2, by = 0.1),
                              algorithm = 4,
                              random.seed = 1234)

banksy_tree_plot <- clustree(xen_bc_banksy, 
                              prefix = "BANKSY_snn_res.")

# ---------------- Banksy UMAP-----------------------------

xen_bc_banksy <- RunUMAP(xen_bc_banksy, dims = 1:20, reduction = "pca.banksy", reduction.name = "umap.banksy")

# ------------------ Banksy Dim Plots ---------------------

umap_clusters_banksy <- DimPlot(xen_bc_banksy, group.by = "BANKSY_snn_res.0.7", label = TRUE)

image_clusters_banksy <- ImageDimPlot(xen_bc_banksy, group.by = "BANKSY_snn_res.0.7", dark.background = F,
                              size = 0.5)

# -------------------- Save plots ------------------------- 

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(qc_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    device   = "jpeg",
    dpi = dpi
  )
}

save_plot_clustree <- function(plot, filename, width = 12, height = 8, dpi = 300) {
  ggsave(
    filename = file.path(qc_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    device   = "jpeg",
    dpi = dpi
  )
}

# Clustree Plot

save_plot_clustree(cluster_tree_plot, "plot38_clustree_Xenium.jpeg")
save_plot_clustree(banksy_tree_plot, "plot39_banksy_clustree_Xenium.jpeg")

# UMAP

save_plot(umap_clusters, "plot40_umap_clusters_Xenium.jpeg")
save_plot(umap_clusters_banksy, "plot41_umap_clusters_banksy_Xenium.jpeg")

save_plot(image_clusters, "plot42_image_clusters_Xenium.jpeg")
save_plot(image_clusters_banksy, "plot43_image_clusters_banksy_Xenium.jpeg")

cat("Done. Neighbors-UMAP  plots saved in:", qc_dir, "\n")
saveRDS(xen_bc_banksy, file = "xen_bc_umap.rds")
