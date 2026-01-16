#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Normalization, Data Scaling, Principal
# Component Analysis (PCA), and Plots
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/norm_pca_plots"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("Norm-PCA plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc_sub <- readRDS(file = "/home/courses/student75/xen_bc_sub.rds")

# -------------------- Normalize Data ---------------------

xen_bc_sub <- NormalizeData(xen_bc_sub, scale.factor = median(xen_bc_sub$nCount_Xenium))

# ------------------- Find Variable Features --------------

xen_bc_sub <- FindVariableFeatures(xen_bc_sub,
                                   nfeatures = 1500)

top10 <- head(VariableFeatures(xen_bc_sub), 10)

var_plot <- VariableFeaturePlot(xen_bc_sub)

var_plot_final <- LabelPoints(plot = var_plot, points = top10, repel = TRUE)

# ------------------- Scale Data --------------------------

xen_bc_sub <- ScaleData(xen_bc_sub)

# ------------------- Principal Component Analysis --------

xen_bc_sub <- RunPCA(xen_bc_sub)

vizdim_plot <- VizDimLoadings(xen_bc_sub, dims = 1:2, reduction = "pca", balanced = TRUE)

pca_plot <- DimPlot(xen_bc_sub, reduction = "pca", dims = 1:2) + NoLegend()

pca_feat1 <- FeaturePlot(xen_bc_sub,
                         reduction = "pca",
                         dims = 1:2,
                         cols = c("white", "red"),
                         features = "XBP1")

pca_feat2 <- FeaturePlot(xen_bc_sub,
                         reduction = "pca",
                         dims = 1:2,
                         cols = c("white", "red"),
                         features = "CTSL")

pca_heatmap <- DimHeatmap(xen_bc_sub, 
                          cells = 500, 
                          balanced = TRUE, 
                          dims = 1,
                          fast = FALSE) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) + 
  ggtitle("Heatmap PC1")

# ------------------- Elbow Plot --------------------------

elbow_sub <- ElbowPlot(xen_bc_sub, ndims = 50)

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

# Var Plot

save_plot(var_plot_final, "plot25_var_plot_Xenium.jpeg")

# Viz Dim Plot

save_plot(vizdim_plot, "plot26_vizdim_plot_Xenium.jpeg")

# PCA Plot

save_plot(pca_plot, "plot27_pca_plot_Xenium.jpeg")
save_plot(pca_feat1, "plot28_pca_feat1_Xenium.jpeg")
save_plot(pca_feat2, "plot29_pca_feat2_Xenium.jpeg")

# PCA heatmap

save_plot(pca_heatmap, "plot30_pc1_heatmap_Xenium.jpeg")

# Elbow Plot

save_plot(elbow_sub, "plot31_elbow_plot_Xenium.jpeg")

cat("Done. Norm-PCA plots saved in:", qc_dir, "\n")
saveRDS(xen_bc_sub, file = "xen_bc_pca.rds")
