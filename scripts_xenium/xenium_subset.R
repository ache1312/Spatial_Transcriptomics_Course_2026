#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - QC subset and plots
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/qc_plots_subset"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("Subset plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc <- readRDS(file = "/home/courses/student75/xen_bc.rds")

# -------------------- Subset data ------------------------

xen_bc_sub <- subset(xen_bc, 
                     subset = nCount_Xenium > 10 &
                     nFeature_Xenium > 10 &
                     percent_BlankCodeword < 1 &
                     percent_ControlCodeword < 1 &
                     percent_ControlProbe < 1 &
                     percent_GenomicControl < 1)

# -------------------- Keep only cells from ROI -----------

subset_cells <- read.csv("/home/courses/student75/Xenium_cells_workshop.csv", skip = 2, header = TRUE)

xen_bc_sub <- subset(xen_bc_sub, cells = subset_cells$Cell.ID)

# -------------------- Compare full vs subset datasets ----

xen_bc
xen_bc_sub

# -------------------- Generate plots ---------------------

# Violin plots

vln_count <- VlnPlot(xen_bc_sub, 
                     features = "nCount_Xenium", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_feature <- VlnPlot(xen_bc_sub, 
                     features = "nFeature_Xenium", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_blank <- VlnPlot(xen_bc_sub, 
                     features = "percent_BlankCodeword", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_ctrlcode <- VlnPlot(xen_bc_sub, 
                     features = "percent_ControlCodeword", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_ctrlprobe <- VlnPlot(xen_bc_sub, 
                     features = "percent_ControlProbe", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_genomicctrl <- VlnPlot(xen_bc_sub, 
                     features = "percent_GenomicControl", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# Spatial plots

img_count <- ImageFeaturePlot(xen_bc_sub, 
                              features = "nCount_Xenium", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red"))

img_feature <- ImageFeaturePlot(xen_bc_sub, 
                              features = "nFeature_Xenium", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_blank <- ImageFeaturePlot(xen_bc_sub, 
                              features = "percent_BlankCodeword", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_ctrlcode <- ImageFeaturePlot(xen_bc_sub, 
                              features = "percent_ControlCodeword", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_ctrlprobe <- ImageFeaturePlot(xen_bc_sub, 
                              features = "percent_ControlProbe", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_genomicctrl <- ImageFeaturePlot(xen_bc_sub, 
                              features = "percent_GenomicControl", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 
                                                           
# -------------------- Save plots -------------------------

save_vln_plot <- function(plot, filename, width = 4, height = 4, dpi = 300) {
  ggsave(
    filename = file.path(qc_dir, filename),
    plot = plot,
    width = width,
    height = height,
    device = "jpeg",
    dpi = dpi
  )
}

save_image_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(qc_dir, filename),
    plot = plot,
    width = width,
    height = height,
    device = "jpeg",
    dpi = dpi
  )
}

# Violin plots

save_vln_plot(vln_count, "plot13_vln_nCount_Xenium_qc.jpeg")
save_vln_plot(vln_feature, "plot14_vln_nFeature_Xenium_qc.jpeg")
save_vln_plot(vln_blank, "plot15_vln_percent_blank_Xenium_qc.jpeg")
save_vln_plot(vln_ctrlcode, "plot16_vln_percent_ctrlcode_Xenium_qc.jpeg")
save_vln_plot(vln_ctrlprobe, "plot17_vln_percent_ctrlprobe_Xenium_qc.jpeg")
save_vln_plot(vln_genomicctrl, "plot18_vln_genomicctrl_Xenium_qc.jpeg")

# ImageFeaturePlots

save_image_plot(img_count, "plot19_img_nCount_Xenium_qc.jpeg")
save_image_plot(img_feature, "plot20_img_nFeature_Xenium_qc.jpeg")
save_image_plot(img_blank, "plot21_img_percent_blank_Xenium_qc.jpeg")
save_image_plot(img_ctrlcode, "plot22_img_percent_ctrlcode_Xenium_qc.jpeg")
save_image_plot(img_ctrlprobe, "plot23_img_percent_ctrlprobe_Xenium_qc.jpeg")
save_image_plot(img_genomicctrl, "plot24_img_genomicctrl_Xenium_qc.jpeg")

cat("Done. Subset plots saved in:", qc_dir, "\n")
saveRDS(xen_bc_sub, file = "xen_bc_sub.rds")
