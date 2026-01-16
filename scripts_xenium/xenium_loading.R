#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Data loading and QC Plots
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)
library(arrow)

# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/qc_plots"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("QC plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc <- LoadXenium(data.dir = "/home/courses/student75/breast_10x")
xen_bc$orig.ident = "Xenium_BC"

# -------------------- Quality metrics --------------------

xen_bc$nCount_total_all <- with(xen_bc@meta.data,
                                nCount_Xenium + 
                                  nCount_BlankCodeword + 
                                  nCount_ControlCodeword + 
                                  nCount_ControlProbe + 
                                  nCount_GenomicControl)

xen_bc$percent_BlankCodeword <- 100 * xen_bc$nCount_BlankCodeword / (xen_bc$nCount_total_all + 1)
xen_bc$percent_ControlCodeword <- 100 * xen_bc$nCount_ControlCodeword / (xen_bc$nCount_total_all + 1)
xen_bc$percent_ControlProbe <- 100 * xen_bc$nCount_ControlProbe / (xen_bc$nCount_total_all + 1)
xen_bc$percent_GenomicControl <- 100 * xen_bc$nCount_GenomicControl / (xen_bc$nCount_total_all + 1)

qc_cols <- c("nCount_Xenium",
             "nFeature_Xenium",
             "percent_BlankCodeword",
             "percent_ControlCodeword",
             "percent_ControlProbe",
             "percent_GenomicControl")

qc_table <- xen_bc@meta.data[, qc_cols]

qc_summary <- t(sapply(qc_table, function(x)
  c(min    = min(x, na.rm = TRUE),
    q1     = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    mean   = mean(x, na.rm = TRUE),
    q3     = quantile(x, 0.75, na.rm = TRUE),
    max    = max(x, na.rm = TRUE))))

# -------------------- Print metrics ----------------------

xen_bc

options(scipen = 999)
qc_summary

# -------------------- Generate plots ---------------------

# Violin plots

vln_count <- VlnPlot(xen_bc, 
                     features = "nCount_Xenium", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_feature <- VlnPlot(xen_bc, 
                     features = "nFeature_Xenium", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_blank <- VlnPlot(xen_bc, 
                     features = "percent_BlankCodeword", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_ctrlcode <- VlnPlot(xen_bc, 
                     features = "percent_ControlCodeword", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_ctrlprobe <- VlnPlot(xen_bc, 
                     features = "percent_ControlProbe", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

vln_genomicctrl <- VlnPlot(xen_bc, 
                     features = "percent_GenomicControl", 
                     pt.size = 0.015) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# Spatial plots

img_count <- ImageFeaturePlot(xen_bc, 
                              features = "nCount_Xenium", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red"))

img_feature <- ImageFeaturePlot(xen_bc, 
                              features = "nFeature_Xenium", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_blank <- ImageFeaturePlot(xen_bc, 
                              features = "percent_BlankCodeword", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_ctrlcode <- ImageFeaturePlot(xen_bc, 
                              features = "percent_ControlCodeword", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_ctrlprobe <- ImageFeaturePlot(xen_bc, 
                              features = "percent_ControlProbe", 
                              max.cutoff = "q95", 
                              dark.background = F,
                              size = 0.25,
                              cols = c("white", "red")) 

img_genomicctrl <- ImageFeaturePlot(xen_bc, 
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

save_vln_plot(vln_count, "plot1_vln_nCount_Xenium.jpeg")
save_vln_plot(vln_feature, "plot2_vln_nFeature_Xenium.jpeg")
save_vln_plot(vln_blank, "plot3_vln_percent_blank_Xenium.jpeg")
save_vln_plot(vln_ctrlcode, "plot4_vln_percent_ctrlcode_Xenium.jpeg")
save_vln_plot(vln_ctrlprobe, "plot5_vln_percent_ctrlprobe_Xenium.jpeg")
save_vln_plot(vln_genomicctrl, "plot6_vln_genomicctrl_Xenium.jpeg")

# ImageFeaturePlots

save_image_plot(img_count, "plot7_img_nCount_Xenium.jpeg")
save_image_plot(img_feature, "plot8_img_nFeature_Xenium.jpeg")
save_image_plot(img_blank, "plot9_img_percent_blank_Xenium.jpeg")
save_image_plot(img_ctrlcode, "plot10_img_percent_ctrlcode_Xenium.jpeg")
save_image_plot(img_ctrlprobe, "plot11_img_percent_ctrlprobe_Xenium.jpeg")
save_image_plot(img_genomicctrl, "plot12_img_genomicctrl_Xenium.jpeg")

cat("Done. QC plots saved in:", qc_dir, "\n")
saveRDS(xen_bc, file = "xen_bc.rds")
