#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Colocalization
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)
library(SpatialExperiment)
library(hoodscanR)
library(pheatmap)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/colocalization_plots"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("Colocalization plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc_umap <- readRDS(file = "/home/courses/student75/xen_bc_umap.rds")

DefaultAssay(xen_bc_umap) = "Xenium"

# -------------------- Spatial experiment object-----------

counts_mat <- GetAssayData(xen_bc_umap, assay = "Xenium",
                           slot = "counts")
data_mat <- GetAssayData(xen_bc_umap, assay = "Xenium",
                         slot = "data")

cell_metadata <- xen_bc_umap@meta.data
gene_metadata <- data.frame(gene = rownames(xen_bc_umap), row.names(xen_bc_umap))

coords <- GetTissueCoordinates(xen_bc_umap)
coords <- coords[, c("x", "y")]
coords$x <- as.numeric(coords$x)
coords$y <- as.numeric(coords$y)
coords <- as.matrix(coords)

spe <- SpatialExperiment(assays = list(counts = counts_mat,
                                       logcounts = data_mat),
                         colData = cell_metadata,
                         rowData = gene_metadata,
                         spatialCoords = coords)

# -------------------- hoodscanR --------------------------

sqe <- readHoodData(spe, anno_col="first_type")
nbs <- findNearCells(sqe, k=100)
mtx <- scanHoods(nbs$distance)      
grp <- mergeByGroup(mtx, nbs$cells) 
sqe <- mergeHoodSpe(sqe, grp)

# -------------------- Generate plots ---------------------

cor <- plotColocal(sqe, pm_cols=colnames(grp), return_matrix=TRUE)
diag(cor) <- NA
pal <- colorRampPalette(rev(hcl.colors(9, "Roma")))(100)
cor_plot <- pheatmap(cor, 
                     cellwidth=15, cellheight=15, 
                     treeheight_row=5, treeheight_col=5,
                     col=pal)

# -------------------- Save plots -------------------------

save_plot <- function(plot, filename, width = 6, height = 4, dpi = 300) {
  ggsave(
    filename = file.path(qc_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    device   = "jpeg",
    dpi = dpi
  )
}

# Heatmaps

save_plot(cor_plot, "plot63_cor_plot_Xenium.jpeg")


cat("Done. Colocalization plots saved in:", qc_dir, "\n")
