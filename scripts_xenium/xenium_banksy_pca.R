#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Banksy - PCA
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(ggplot2)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/bansky_plots"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("Banksy-PCA plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc_pca <- readRDS(file = "/home/courses/student75/xen_bc_pca.rds")

# -------------------- Banksy -----------------------------

xen_bc_banksy <- RunBanksy(xen_bc_pca, lambda = 0.2, verbose = TRUE, assay = "Xenium", slot = "data", features = "variable", k_geom = 15)

DefaultAssay(xen_bc_banksy) <- "BANKSY"

xen_bc_banksy <- RunPCA(xen_bc_banksy, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(xen_bc_banksy))

banksy_vizdim_plot <- VizDimLoadings(xen_bc_banksy, dims = 1:2, reduction = "pca.banksy")

banksy_pca_plot <- DimPlot(xen_bc_banksy, reduction = "pca.banksy", dims = 1:2) + NoLegend()

banksy_pca_feat1 <- FeaturePlot(xen_bc_banksy,
                         reduction = "pca.banksy",
                         dims = 1:2,
                         cols = c("white", "red"),
                         features = "XBP1")

banksy_pca_feat2 <- FeaturePlot(xen_bc_banksy,
                         reduction = "pca.banksy",
                         dims = 1:2,
                         cols = c("white", "red"),
                         features = "CTSL")

banksy_pca_heatmap <- DimHeatmap(xen_bc_banksy, 
                          cells = 500, 
                          reduction = "pca.banksy",
                          assay = "BANKSY",
                          balanced = TRUE, 
                          dims = 1,
                          fast = FALSE) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) + 
  ggtitle("Heatmap PC1 Banksy")

# ------------------- Elbow Plots -------------------------

banksy_elbow_plot <- ElbowPlot(xen_bc_banksy, reduction = "pca.banksy", ndims = 50)

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

# Viz Dim Plot

save_plot(banksy_vizdim_plot, "plot32_banksy_vizdim_plot_Xenium.jpeg")

# PCA Plot

save_plot(banksy_pca_plot, "plot33_banksy_pca_plot_Xenium.jpeg")
save_plot(banksy_pca_feat1, "plot34_banksy_pca_feat1_Xenium.jpeg")
save_plot(banksy_pca_feat2, "plot35_banksy_pca_feat2_Xenium.jpeg")

# PCA heatmap

save_plot(banksy_pca_heatmap, "plot36_banksy_pc1_heatmap_Xenium.jpeg")

# Elbow Plot

save_plot(banksy_elbow_plot, "plot37_banksy_elbow_plot_Xenium.jpeg")

cat("Done. Banksy-PCA plots saved in:", qc_dir, "\n")
saveRDS(xen_bc_banksy, file = "xen_bc_banksy.rds")


