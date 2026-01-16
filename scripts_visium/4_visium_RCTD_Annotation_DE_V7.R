print("Loading Libraries")
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
  library(RhpcBLASctl)
})

options(future.globals.maxSize = 1e10)

# -------------------- Threads (safe) --------------------
blas_set_num_threads(1)

# -------------------- Plot output --------------------
plot_dir <- file.path(Sys.getenv("HOME"), "results/visium/figures")
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

plot_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
}

print("Libraries Loaded successfully")

#### =========================================================
####  Cell type annotation with RCTD (SUBSAMPLE QUERY + COARSE REFERENCE)
#### =========================================================

# === Load objects ===
print("Load Objects")

crc5_filt <- readRDS(file.path(Sys.getenv("HOME"), "results/visium/objects/visium_crc5_banksy_subset.rds"))
crc_integrated <- readRDS(file.path(Sys.getenv("HOME"), "data/crc_integrated.rds"))

print("Objects Loaded Successfully")

# -------------------- Add cluster annotations to reference (crc_integrated) --------------------
# Make sure seurat_clusters is numeric for mapping
crc_integrated$seurat_clusters <- as.integer(as.character(crc_integrated$seurat_clusters))

crc_integrated@meta.data <- crc_integrated@meta.data %>%
  mutate(Annotation = case_when(
    seurat_clusters == 1  ~ "Tumor_cells_1",
    seurat_clusters == 7  ~ "Endothelial",
    seurat_clusters == 16 ~ "Lymphatic",
    seurat_clusters == 3  ~ "Fibroblasts_1",
    seurat_clusters == 15 ~ "Fibroblasts_2",
    seurat_clusters == 19 ~ "Glial",
    seurat_clusters == 10 ~ "Smooth muscle",
    seurat_clusters == 11 ~ "Tumor_cells_2",
    seurat_clusters == 14 ~ "B",
    seurat_clusters == 4  ~ "Plasma_1",
    seurat_clusters == 9  ~ "Plasma_2",
    seurat_clusters == 12 ~ "Plasma_3",
    seurat_clusters == 5  ~ "CD4_T",
    seurat_clusters == 6  ~ "CD8_T",
    seurat_clusters == 2  ~ "Macro_Mono_1",
    seurat_clusters == 13 ~ "Macro_Mono_2",
    seurat_clusters == 18 ~ "cDC1",
    seurat_clusters == 20 ~ "mregDC",
    seurat_clusters == 17 ~ "Mast",
    seurat_clusters == 8  ~ "Mixed_cells",
    TRUE ~ "Unknown"
  ))

crc_integrated@meta.data$Annotation <- factor(
  crc_integrated@meta.data$Annotation,
  levels = rev(c(
    "Tumor_cells_1","Endothelial","Lymphatic","Fibroblasts_1","Fibroblasts_2","Glial",
    "Smooth muscle","Tumor_cells_2","B","Plasma_1","Plasma_2","Plasma_3","CD4_T","CD8_T",
    "Macro_Mono_1","Macro_Mono_2","cDC1","mregDC","Mast","Mixed_cells","Unknown"
  ))
)

#### ---------------------------------------------------------
#### Option 2: SUBSAMPLE QUERY BINS (to speed up RCTD)
#### ---------------------------------------------------------
print("Preparing Query Data (Spatial.008um)")

DefaultAssay(crc5_filt) <- "Spatial.008um"

# robust counts extraction
counts_hd_all <- GetAssayData(crc5_filt, assay = "Spatial.008um", slot = "counts")
cells_all <- colnames(counts_hd_all)

# Choose how many bins to use for RCTD
# - Set QUERY_N <- NULL to use FULL query (all bins)
# - Recommended: 10000–20000 for fast runs
QUERY_N <- 10000

if (!is.null(QUERY_N) && QUERY_N < length(cells_all)) {
  set.seed(1234)
  cells_use <- sample(cells_all, size = QUERY_N)
  message("Subsampling query bins: ", QUERY_N, " / ", length(cells_all))
} else {
  cells_use <- cells_all
  message("Using FULL query bins: ", length(cells_use))
}

counts_hd <- counts_hd_all[, cells_use, drop = FALSE]
coords_all <- GetTissueCoordinates(crc5_filt)[, 1:2, drop = FALSE]
coords <- coords_all[cells_use, , drop = FALSE]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# IMPORTANT: If we subsample, we must work downstream on a matching Seurat object
crc5_query <- if (identical(cells_use, cells_all)) crc5_filt else subset(crc5_filt, cells = cells_use)

#### ---------------------------------------------------------
#### Option 3: COARSEN REFERENCE LABELS (to speed up RCTD)
#### ---------------------------------------------------------
print("Preparing Reference Data (crc_integrated; COARSE labels)")

DefaultAssay(crc_integrated) <- "RNA"
counts_ref <- GetAssayData(crc_integrated, slot = "counts")
nUMI_ref <- crc_integrated$nCount_RNA

# Create a coarser annotation (edit mapping as needed)
crc_integrated$Annotation_coarse <- as.character(crc_integrated$Annotation)

crc_integrated$Annotation_coarse <- dplyr::recode(
  crc_integrated$Annotation_coarse,
  "Tumor_cells_1" = "Tumor",
  "Tumor_cells_2" = "Tumor",

  "Plasma_1" = "Plasma",
  "Plasma_2" = "Plasma",
  "Plasma_3" = "Plasma",

  "Macro_Mono_1" = "Macro_Mono",
  "Macro_Mono_2" = "Macro_Mono",

  "Fibroblasts_1" = "Fibroblasts",
  "Fibroblasts_2" = "Fibroblasts",

  "Endothelial" = "Vascular",
  "Lymphatic"   = "Vascular"
)

# keep as-is: CD4_T, CD8_T, B, cDC1, mregDC, Mast, Glial, Smooth muscle, Mixed_cells, Unknown
cluster <- droplevels(as.factor(crc_integrated$Annotation_coarse))

message("Reference (coarse) number of types: ", nlevels(cluster))
print(table(cluster))

reference <- Reference(counts_ref, cluster, nUMI_ref)

# cores (match Slurm if present)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
message("Using ", ncores, " cores (SLURM_CPUS_PER_TASK).")
options(mc.cores = ncores)

# Run RCTD
print("Creating RCTD")
RCTD <- create.RCTD(query, reference, max_cores = ncores, test_mode = FALSE)

print("RCTD Created")
print("Running RCTD")
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

print("RCTD finished successfully")

# -------------------- SAVE CHECKPOINT (after RCTD) --------------------
saveRDS(
  RCTD,
  file = file.path(Sys.getenv("HOME"),"results/visium/objects/RCTD_crc5_subset_doublet_COARSE_ref_QUERYsubsample.rds")
)
print("Saved RCTD: ~/results/visium/objects/RCTD_crc5_subset_doublet_COARSE_ref_QUERYsubsample.rds")

print("Adding Metadata")

# Add results to the matching object (crc5_query)
res_df <- RCTD@results$results_df
stopifnot(all(rownames(res_df) %in% colnames(crc5_query)))
crc5_query <- AddMetaData(crc5_query, metadata = res_df)

crc5_query$first_type <- as.character(crc5_query$first_type)
crc5_query$first_type[is.na(crc5_query$first_type)] <- "Unknown"

# Keep downstream naming consistent
crc5_query$full_first_type <- crc5_query$first_type
Idents(crc5_query) <- "full_first_type"

# ---- Plot 13: Spatial annotation ----
plot_hd13 <- SpatialDimPlot(
  crc5_query,
  group.by = "full_first_type",
  pt.size.factor = 6,
  label = TRUE,
  repel = TRUE,
  label.size = 4
) + ggtitle("RCTD annotation (QUERY subsample + COARSE reference)") + plot_theme()

save_png(plot_hd13, "13-RCTD_Spatial_Annotation.png", 13, 6)
rm(plot_hd13); gc()

#### =========================================================
#### Differential expression analysis (on crc5_query)
#### =========================================================

DefaultAssay(crc5_query) <- "Spatial.008um"
Idents(crc5_query) <- "full_first_type"

print("Differential expression Analysis Started")
print("Finding All Markers for each Annotated Group")

hd_markers <- FindAllMarkers(
  crc5_query,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "Spatial.008um"
)

print("Saving markers table")
write.csv(hd_markers, file = file.path(plot_dir, "16b-RCTD_markers_table.csv"))
print("Markers saved")

# ---- DotPlot (14) ----
top_markers <- hd_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

features_use <- unique(top_markers$gene)

plot_hd14 <- DotPlot(
  crc5_query,
  features = features_use,
  group.by = "full_first_type",
  assay = "Spatial.008um"
) + RotatedAxis() + ggtitle("Top markers (DotPlot)") + plot_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save_png(plot_hd14, "14-RCTD_DotPlot.png", 18, 10)
rm(plot_hd14); gc()

# ---- ViolinPlot (15) ----
plot_hd15 <- VlnPlot(
  crc5_query,
  features = features_use,
  group.by = "full_first_type",
  assay = "Spatial.008um",
  pt.size = 0
) + ggtitle("Top markers (VlnPlot)") + plot_theme()

save_png(plot_hd15, "15-RCTD_ViolinPlot.png", 20, 30)
rm(plot_hd15); gc()

# ---- Heatmap (16) ----
plot_hd16 <- DoHeatmap(
  crc5_query,
  features = features_use,
  group.by = "full_first_type",
  assay = "Spatial.008um"
) + ggtitle("Top markers (Heatmap)") + plot_theme()

save_png(plot_hd16, "16-RCTD_Heatmap.png", 20, 10)
rm(plot_hd16); gc()

print("Pipeline finished successfully")

# Save final object for next steps (CellChat / downstream)
saveRDS(crc5_query, file.path(Sys.getenv("HOME"), "results/visium/objects/visium_crc5_query_annotated_rctd_COARSEref_QUERYsubsample.rds"))
print("Saved: ~/results/visium/objects/visium_crc5_query_annotated_rctd_COARSEref_QUERYsubsample.rds")

