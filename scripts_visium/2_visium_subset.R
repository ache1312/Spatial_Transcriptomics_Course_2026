#### Visium HD ####

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

## -------------------- GLOBAL PLOT SETTINGS --------------------
## Keep fonts and sizes consistent across all images
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

plot_dir <- file.path(Sys.getenv("HOME"),"results/visium/figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(plot, filename, width, height, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE
  )
}


## -------------------- DATA LOADING -------------------------

bc_visium <- readRDS(file = file.path(Sys.getenv("HOME"), "results/visium/objects/bc_visium.rds"))

## -------------------- Data SUBSETTING ----------------------

bc_visium_subset <- subset(bc_visium, subset = nCount_Spatial> 100 & nFeature_Spatial > 100)

subset_cells <- read.csv(file.path(Sys.getenv("HOME"), "data/visium/roi_mama_vis.csv"), header = TRUE)

bc_visium_subset <- subset(bc_visium_subset, cells = subset_cells$Barcode)

bc_visium
bc_visium_subset

## -------------------- QC PLOTS (Post-QC) --------------------

## Plot Counts (UMIs)
plot_hd1 <- VlnPlot(bc_visium_subset, features = "nCount_Spatial", pt.size = 0) + xlab("") + plot_theme()
plot_hd2 <- SpatialFeaturePlot(bc_visium_subset, features = "nCount_Spatial", pt.size.factor = 8) + plot_theme()
plot_hd3 <- bc_visium_subset@meta.data %>%
  ggplot(aes(x = nCount_Spatial)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  ylab("Cell density") +
  xlab("Number of Counts (UMIs) / bin") +
  ggtitle("Pre-QC UMIs/Bin BC") +
  geom_vline(xintercept = 100, linetype = "dashed") +
  plot_theme()

combined_plot <- plot_hd1 | plot_hd2 | plot_hd3
save_png(combined_plot, "04-BC_PosQC_UMIs.png", width = 14, height = 4.5)
rm(plot_hd1, plot_hd2, plot_hd3, combined_plot); gc()

## Plot Features (genes)
plot_hd4 <- VlnPlot(bc_visium_subset, features = "nFeature_Spatial", pt.size = 0) + xlab("") + plot_theme()
plot_hd5 <- SpatialFeaturePlot(bc_visium_subset, features = "nFeature_Spatial", pt.size.factor = 8) + plot_theme()
plot_hd6 <- bc_visium_subset@meta.data %>%
  ggplot(aes(x = nFeature_Spatial)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  ylab("Cell density") +
  xlab("Number of Genes / bin") +
  ggtitle("Pre-QC Genes/Bin BC") +
  geom_vline(xintercept = 100, linetype = "dashed") +
  plot_theme()

combined_plot <- plot_hd4 | plot_hd5 | plot_hd6
save_png(combined_plot, "05-BC_PosQC_Genes_bin.png", width = 14, height = 4.5)
rm(plot_hd4, plot_hd5, plot_hd6, combined_plot); gc()

## Plot mitochondrial percentage
bc_visium_subset <- PercentageFeatureSet(bc_visium_subset, pattern = "^MT-", col.name = "mito_percent")

plot_hd7 <- VlnPlot(bc_visium_subset, features = "mito_percent", pt.size = 0) + xlab("") + plot_theme()
plot_hd8 <- SpatialFeaturePlot(bc_visium_subset, features = "mito_percent", pt.size.factor = 8) + plot_theme()
plot_hd9 <- bc_visium_subset@meta.data %>%
  ggplot(aes(x = mito_percent)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  ylab("Cell density") +
  xlab("Mito % / bin") +
  ggtitle("Pre-QC Mitochondrial percentage/Bin BC") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  plot_theme()

combined_plot <- plot_hd7 | plot_hd8 | plot_hd9
save_png(combined_plot, "06-BC_PosQC_Mitochondria.png", width = 14, height = 4.5)
rm(plot_hd7, plot_hd8, plot_hd9, combined_plot); gc()

## -------------------- Save file --------------------

saveRDS(bc_visium_subset, file = file.path(Sys.getenv("HOME"), "results/visium/objects/bc_visium_subset.rds"))
