#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Annotation - Deconvolution
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)
library(spacexr)
library(RhpcBLASctl)
library(cowplot)
library(tidyverse)
options(future.globals.maxSize = 1e10)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/annot_plots"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("Annotation plots will be saved in:", qc_dir, "\n")

# -------------------- Multi-core Setup --------------------

# Set environment variables telling math libraries to use 10 threads
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  R_MAX_NUM_THREADS = "1"
)

# Set BLAS threads
blas_set_num_threads(1)

# Set OpenMP threads
omp_set_num_threads(1)

# -------------------- Data loading -----------------------

xen_bc_umap <- readRDS(file = "/home/courses/student75/xen_bc_umap.rds")
DefaultAssay(xen_bc_umap) = "Xenium"

atlas_breast <- readRDS("/home/courses/student75/atlas_liu.rds")

# --------------------- DEGs ------------------------------

degs_bc <- FindAllMarkers(xen_bc_umap,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              min.pct = 0.25,
                              group.by = "BANKSY_snn_res.0.7")

degs_bc_sign <- degs_bc %>%
  filter(p_val_adj < 0.05)

write.csv(degs_bc_sign, file = "degs_bc_sign.csv")

top5 <- degs_bc_sign %>%
  group_by(cluster) %>%
  slice_head(n = 3) %>%
  ungroup()

dotplot_degs_top5 <- DotPlot(xen_bc_umap,
                             group.by = "BANKSY_snn_res.0.7",
                             dot.min = 0.1,
                             features = unique(top5$gene)) +
  RotatedAxis()

bc_genes <- c("EPCAM", "KRT8", "KRT18",
              "MKI67", "TOP2A", "UBE2C",
              "COL1A1", "COL1A2", "DCN",
              "ACTA2", "TAGLN", "MYL9",
              "PECAM1", "VWF", "KDR",
              "PROX1", "PDPN", "LYVE1",
              "RGS5", "CSPG4", "MCAM",
              "LYZ", "LST1", "TYROBP",
              "C1QA", "C1QB", "APOE",
              "HLA-DRA","CD74",
              "XCR1", "BATF3", "CLEC9A",
              "FCER1A", "CD1C", "CLEC10A",
              "CD3D", "CD3E", "CD8A", "CD4",
              "IL7R", "CCR7",
              "FOXP3", "IL2RA", "CTLA4",
              "NKG7", "GNLY", "GZMB", "PRF1",
              "PDCD1", "LAG3", "HAVCR2", "TOX",
              "CD69", "ITGAE",
              "MS4A1", "CD79A", 
              "MZB1", "XBP1", "JCHAIN",
              "TPSB2", "KIT", "MS4A2")

dotplot_markers <- DotPlot(xen_bc_umap,
                           group.by = "BANKSY_snn_res.0.7",
                           dot.min = 0.1,
                           features = bc_genes) +
  RotatedAxis()

img_trac <- ImageFeaturePlot(xen_bc_umap, 
                             features = "TRAC", 
                             max.cutoff = "q95", 
                             dark.background = F,
                             size = 0.5,
                             cols = c("white", "red"))
img_cd8a <- ImageFeaturePlot(xen_bc_umap, 
                             features = "CD8A", 
                             max.cutoff = "q95", 
                             dark.background = F,
                             size = 0.5,
                             cols = c("white", "red"))
img_epcam <- ImageFeaturePlot(xen_bc_umap, 
                              features = "EPCAM", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.5,
                              cols = c("white", "red"))
img_ms4a1 <- ImageFeaturePlot(xen_bc_umap, 
                              features = "MS4A1", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.5,
                              cols = c("white", "red"))
img_col4a1 <- ImageFeaturePlot(xen_bc_umap, 
                               features = "COL4A1", 
                               max.cutoff = "q95", 
                               dark.background = F,
                               size = 0.5,
                               cols = c("white", "red"))

# --------------------- RCTD annotation -------------------

# --------------------- Subset atlas ----------------------

atlas_common_genes <- intersect(rownames(xen_bc_umap), rownames(atlas_breast))
cat("Number of shared genes between Xenium and atlas:")
print(length(atlas_common_genes))

atlas_subset <- CreateSeuratObject(counts = atlas_breast[["RNA"]]$counts[atlas_common_genes,],
                                   meta = atlas_breast@meta.data) %>%
                                   NormalizeData() %>%
                                   FindVariableFeatures() %>%
                                   ScaleData() %>%
                                   RunPCA()
                                           
# -------------------- Prepare Query data (Xenium) --------

counts_xen <- GetAssayData(xen_bc_umap,
                           assay = "Xenium",
                           layer = "counts")

counts_xen <- as.matrix(counts_xen)

cells_xen <- colnames(counts_xen)

coords_all <- GetTissueCoordinates(xen_bc_umap)

coords_all <- as.data.frame(coords_all)

if ("cell" %in% colnames(coords_all)) {
  rownames(coords_all) <- coords_all$cell
} else if ("barcode" %in% colnames(coords_all)) {
  rownames(coords_all) <- coords_all$barcode
}

coords_xen <- coords_all[cells_xen, c("x", "y")]

nUMI <- colSums(counts_xen)

query <- SpatialRNA(coords = coords_xen,
                    counts = counts_xen,
                    nUMI   = nUMI)


# --------------------- Prepare reference data (scRNA) -------

counts_sc <- atlas_subset[["RNA"]]$counts

cluster_sc <- as.factor(atlas_subset$Labels)

nUMI_sc <- atlas_subset$nCount_RNA

cluster_sc <- droplevels(cluster_sc)

reference <- Reference(counts_sc, cluster_sc, nUMI_sc)

# --------------------- Run RCTD -------------------------------

RCTD <- create.RCTD(query, reference, max_cores = 20, test_mode = FALSE, CELL_MIN_INSTANCE = 5)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# --------------------- Add RCTD results to Xenium data --------

xen_bc_umap <- AddMetaData(xen_bc_umap, metadata = RCTD@results$results_df)
xen_bc_umap$first_type <- as.character(xen_bc_umap$first_type)
xen_bc_umap$first_type[is.na(xen_bc_umap$first_type)] <- "Unknown"

cat("Number of annotated cells:")
table(xen_bc_umap$first_type)

umap_clusters_rctd <- DimPlot(xen_bc_umap, 
                              group.by = "first_type", 
                              label = TRUE,
                              reduction = "umap.banksy", repel = TRUE)

image_clusters_rctd <- ImageDimPlot(xen_bc_umap, 
                                    group.by = "first_type",
                                    dark.background = FALSE,
                                    size = 0.5)

# --------------------- Plots ----------------------------------

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

save_plot_dotplot <- function(plot, filename, width = 11, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(qc_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    device   = "jpeg",
    dpi = dpi
  )
}

# Dot Plots

save_plot_dotplot(dotplot_degs_top5, "plot44_dotplot_degs_top5_Xenium.jpeg")
save_plot_dotplot(dotplot_markers, "plot45_dotplot_markers_Xenium.jpeg")

# UMAP

save_plot(umap_clusters_rctd, "plot46_umap_clusters_rctd_Xenium.jpeg")

# Image Plot

save_plot(image_clusters_rctd, "plot47_image_clusters_rctd_Xenium.jpeg")
save_plot(img_trac, "plot48_img_trac_Xenium.jpeg")
save_plot(img_cd8a, "plot49_img_cd8a_Xenium.jpeg")
save_plot(img_epcam, "plot50_img_epcam_Xenium.jpeg")
save_plot(img_ms4a1, "plot51_img_ms4a1_Xenium.jpeg")
save_plot(img_col4a1, "plot52_img_col4a1_Xenium.jpeg")

cat("Done. Annotation plots saved in:", qc_dir, "\n")
saveRDS(xen_bc_umap, file = "xen_bc_umap.rds")
