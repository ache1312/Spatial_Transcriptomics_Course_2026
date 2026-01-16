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
## -------------------- Output folder for plots --------------------

plot_dir <- file.path(Sys.getenv("HOME"),"results/visium/figures")
dir.create(plot_dir, recursive = TRUE)

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

## -------------------- DATA LOADING --------------------


bc_visium <- Load10X_Spatial(
  data.dir = file.path(Sys.getenv("HOME"),"breast10x_visium/binned_outputs/square_008um/"))

bc_visium$orig.ident <- "Breast_Cancer"
Idents(bc_visium) <- "orig.ident"

## -------------------- QC PLOTS (PRE-QC) --------------------

## Plot Counts (UMIs)
plot_hd1 <- VlnPlot(bc_visium, features = "nCount_Spatial", pt.size = 0) + xlab("") + plot_theme()
plot_hd2 <- SpatialFeaturePlot(bc_visium, features = "nCount_Spatial", pt.size.factor = 8) + plot_theme()
plot_hd3 <- bc_visium@meta.data %>%
  ggplot(aes(x = nCount_Spatial)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  ylab("Cell density") +
  xlab("Number of Counts (UMIs) / bin") +
  ggtitle("Pre-QC UMIs/Bin BC") +
  geom_vline(xintercept = 100, linetype = "dashed") +
  plot_theme()

combined_plot1 <- plot_hd1 | plot_hd2 | plot_hd3
save_png(combined_plot1, "01-BC_PreQC_UMIs.png")
rm(plot_hd1, plot_hd2, plot_hd3, combined_plot1); gc()

## Plot Features (genes)
plot_hd4 <- VlnPlot(bc_visium, features = "nFeature_Spatial", pt.size = 0) + xlab("") + plot_theme()
plot_hd5 <- SpatialFeaturePlot(bc_visium, features = "nFeature_Spatial", pt.size.factor = 8) + plot_theme()
plot_hd6 <- bc_visium@meta.data %>%
  ggplot(aes(x = nFeature_Spatial)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  ylab("Cell density") +
  xlab("Number of Genes / bin") +
  ggtitle("Pre-QC Genes/Bin BC") +
  geom_vline(xintercept = 100, linetype = "dashed") +
  plot_theme()

combined_plot2 <- plot_hd4 | plot_hd5 | plot_hd6
save_png(combined_plot2, "02-bc_vis_PreQC_Genes_bin.png")
rm(plot_hd4, plot_hd5, plot_hd6, combined_plot2); gc()

## Plot mitochondrial percentage
bc_visium <- PercentageFeatureSet(bc_visium, pattern = "^MT-", col.name = "mito_percent")

plot_hd7 <- VlnPlot(bc_visium, features = "mito_percent", pt.size = 0) + xlab("") + plot_theme()
plot_hd8 <- SpatialFeaturePlot(bc_visium, features = "mito_percent", pt.size.factor = 8) + plot_theme()
plot_hd9 <- bc_visium@meta.data %>%
  ggplot(aes(x = mito_percent)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  ylab("Cell density") +
  xlab("Mito % / bin") +
  ggtitle("Pre-QC Mitochondrial percentage/Bin BC") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  plot_theme()

combined_plot3 <- plot_hd7 | plot_hd8 | plot_hd9
save_png(combined_plot3, "03-bc_vis_PreQC_Mitochondria.png")
rm(plot_hd7, plot_hd8, plot_hd9, combined_plot3); gc()

## -------------------- Save file --------------------

saveRDS(bc_visium, file = file.path(Sys.getenv("HOME"), "results/visium/objects/bc_visium.rds"))
