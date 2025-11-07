##################################################################
#====================== VISIUM HD WORKSHOP ======================#
##################################################################

#Is important to set the seed at the beginning to
#secure the reproducibility of the analysis
set.seed(1234)
#Loading Libraries----
library(Seurat)
library(clustree)
library(Banksy)
library(SeuratWrappers)
library(spacexr)
library(tidyverse)

options(future.globals.maxSize = 1e10)

#Data download and setup ----

## Create Data directory (if it doesn't exist)----
#you can change the direction between "" to your own directory

data_dir <- "../../Data" 
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
  cat("Created Data directory\n")
}

## Download and extract data (if not already present)----
zip_file <- file.path(data_dir, "visium_hd_data.zip")
if (!file.exists(zip_file)) {
  cat("Downloading Visium HD data...\n")
  download.file(
    url = "https://apps.cienciavida.org/single_cell/visium_hd_data.zip",
    destfile = zip_file,
    mode = "wb"
  )
  cat("Download complete!\n")
}

## Extract zip file if not already extracted----
crc5_dir <- file.path(data_dir, "GEO_CRC5_VHD")
if (!dir.exists(crc5_dir)) {
  cat("Extracting data...\n")
  unzip(zip_file, exdir = data_dir)
  cat("Extraction complete!\n")
}

#Data loading----
##VIsium HD Data from Patient 5 (crc5)----
crc5 <- Load10X_Spatial(data.dir = file.path(data_dir, "GEO_CRC5_VHD/outs/"),
                        bin.size = 8) 

crc5$orig.ident <- "CRC5" #Replace the information of orig.ident column in metadata
Idents(crc5) <- "orig.ident"

##VIsium HD Data from Patient 2 (crc2)----
crc2 <- Load10X_Spatial(data.dir = file.path(data_dir, "GEO_CRC2_VHD/outs/"),
                        bin.size = 8) # Visium HD data of patient 2

crc2$orig.ident <- "CRC2"
Idents(crc2) <- "orig.ident"

##Load scRNA-seq atlas of both patients (already processed)----
load(file.path(data_dir, "crc_integrated.RData"))

#Data exploration and quality control (QC)----

##CRC5----
###Plot Counts (UMIs)----

plot_hd1 <- VlnPlot(crc5, 
                    features = "nCount_Spatial.008um", 
                    pt.size = 0) +
  xlab("")

plot_hd2 <- SpatialFeaturePlot(crc5, 
                               features = "nCount_Spatial.008um", 
                               pt.size.factor = 8)

plot_hd3 <- crc5@meta.data %>%
  ggplot(aes(x=nCount_Spatial.008um)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Number of Counts (UMIs) / bin") +
  ggtitle('Pre-QC UMIs/Bin CRC5') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 100,
             linetype = "dashed")

print(plot_hd1 | plot_hd2 | plot_hd3)

###Plot Features (genes)----

plot_hd4 <- VlnPlot(crc5, 
                    features = "nFeature_Spatial.008um", 
                    pt.size = 0) +
  xlab("")

plot_hd5 <- SpatialFeaturePlot(crc5, 
                               features = "nFeature_Spatial.008um", 
                               pt.size.factor = 8)

plot_hd6 <- crc5@meta.data %>%
  ggplot(aes(x=nFeature_Spatial.008um)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Number of Genes / bin") +
  ggtitle('Pre-QC Genes/Bin CRC5') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 100,
             linetype = "dashed")

print(plot_hd4 | plot_hd5 | plot_hd6)

###Plot mitochondrial percentage----

crc5 <- PercentageFeatureSet(crc5,
                             pattern = "^MT-",
                             col.name = "mito_percent")

plot_hd7 <- VlnPlot(crc5, 
                    features = "mito_percent", 
                    pt.size = 0) +
  xlab("")

plot_hd8 <- SpatialFeaturePlot(crc5, 
                               features = "mito_percent", 
                               pt.size.factor = 8)

plot_hd9 <- crc5@meta.data %>%
  ggplot(aes(x=mito_percent)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Mito % / bin") +
  ggtitle('Pre-QC Mitochondrial percentage/Bin CRC5') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 20,
             linetype = "dashed")

print(plot_hd7 | plot_hd8 | plot_hd9)

##CRC2----

###Plot Counts (UMIs)----

plot_hd10 <- VlnPlot(crc2, 
                     features = "nCount_Spatial.008um", 
                     pt.size = 0) +
  xlab("")

plot_hd11 <- SpatialFeaturePlot(crc2, 
                                features = "nCount_Spatial.008um", 
                                pt.size.factor = 8)

plot_hd12 <- crc2@meta.data %>%
  ggplot(aes(x=nCount_Spatial.008um)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Number of Counts (UMIs) / bin") +
  ggtitle('Pre-QC UMIs/Bin CRC2') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 100,
             linetype = "dashed")

print(plot_hd10 | plot_hd11 | plot_hd12)

###Plot Features (genes)----

plot_hd13 <- VlnPlot(crc2, 
                     features = "nFeature_Spatial.008um", 
                     pt.size = 0) +
  xlab("")

plot_hd14 <- SpatialFeaturePlot(crc2, 
                                features = "nFeature_Spatial.008um", 
                                pt.size.factor = 8)

plot_hd15 <- crc2@meta.data %>%
  ggplot(aes(x=nFeature_Spatial.008um)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Number of Genes / bin") +
  ggtitle('Pre-QC Genes/Bin CRC5') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 100,
             linetype = "dashed")

print(plot_hd13 | plot_hd14 | plot_hd15)

###Plot mitochondrial percentage----

crc2 <- PercentageFeatureSet(crc2,
                             pattern = "^MT-",
                             col.name = "mito_percent")

plot_hd16 <- VlnPlot(crc2, 
                     features = "mito_percent", 
                     pt.size = 0) +
  xlab("")

plot_hd17 <- SpatialFeaturePlot(crc2, 
                                features = "mito_percent", 
                                pt.size.factor = 8)

plot_hd18 <- crc2@meta.data %>%
  ggplot(aes(x=mito_percent)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Mito % / bin") +
  ggtitle('Pre-QC Mitochondrial percentage/Bin CRC5') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 20,
             linetype = "dashed")

print(plot_hd16 | plot_hd17 | plot_hd18)

##Subset data---- 
#Here we eliminate those bins with low QCs
crc5_filt <- subset(crc5, 
                    subset = nCount_Spatial.008um > 100 & 
                      nFeature_Spatial.008um > 100)

crc2_filt <- subset(crc2,
                    subset = nCount_Spatial.008um > 100 & 
                      nFeature_Spatial.008um > 100)

rm(crc5, crc2)
gc()

#Data sketch----
#An Sketch is a subset of the data that represents all the information,
#it's useful to lower the memory usage in some steps of the analysis
#On the next steps we are going to do the demonstration for only
#Patient 5 (crc5)

##Sketch of Patient 5 (crc5)----
crc5_filt <- NormalizeData(crc5_filt)

crc5_filt <- FindVariableFeatures(crc5_filt)

#With this step we generate a new assay with the sketch
crc5_filt <- SketchData(crc5_filt,
                        assay = "Spatial.008um",
                        ncells = 50000,
                        method = "LeverageScore",
                        sketched.assay = "sketch",
                        features = VariableFeatures(crc5_filt))
#What would happen if we have a lower or bigger number of cells?

DefaultAssay(crc5_filt) <- "sketch" 

crc5_filt <- FindVariableFeatures(crc5_filt)

crc5_filt <- ScaleData(crc5_filt)

crc5_filt <- RunPCA(crc5_filt, 
                    assay = "sketch", 
                    reduction.name = "pca.sketch")

#Visualization of top10 genes in the dimentional reduction
plot_hd19 <- VizDimLoadings(crc5_filt,
                            dims = 1:2,
                            nfeatures = 10,
                            reduction = "pca.sketch",
                            balanced = TRUE)

plot_hd20 <- DimPlot(crc5_filt,
                     dims = 1:2,
                     reduction = "pca.sketch") +
  geom_vline(xintercept = 0, 
             linetype = "dashed") +
  geom_hline(yintercept = 0,
             linetype = "dashed")

plot_hd21 <- FeaturePlot(crc5_filt,
                         features = c("LCN2", 
                                      "COL3A1", 
                                      "ID1", 
                                      "MUC12"), 
                         ncol = 2,
                         reduction = "pca.sketch",
                         dims = 1:2)

print(plot_hd19 | plot_hd20 | plot_hd21)

#See how the standard deviation change across dimensions
plot_hd22 <- ElbowPlot(crc5_filt, 
                       ndims = 50, 
                       reduction = "pca.sketch")

print(plot_hd22) #What will happen if we change the number of dimensions?

crc5_filt <- FindNeighbors(crc5_filt, 
                           assay = "sketch", 
                           reduction = "pca.sketch",
                           dims = 1:30)
###Subset sketch data to visualize resolutions with a clustree----

cells_sketch <- colnames(crc5_filt[["sketch"]])

crc5_sketch <- subset(crc5_filt, 
                      cells = cells_sketch)

crc5_sketch <- FindClusters(crc5_sketch,
                            resolution = seq(0.1, 
                                             2, 
                                             by = 0.1),
                            algorithm = 4,
                            random.seed = 1234)

#Visualization the clusters to know which resolution is the optimal
plot_hd23 <- clustree(crc5_sketch, prefix = "sketch_snn_res.")
#Resolution 0.8 seems Ok!

print(plot_hd23)

###Run FindClusters in the original seurat object----
#Calculate the resolution of interest (be aware of work in the right assay)

DefaultAssay(crc5_filt) <- "sketch"

crc5_filt <- FindClusters(crc5_filt,
                          resolution = 0.8,
                          algorithm = 4,
                          random.seed = 1234)

#Corroborate that we are having the same information from
#Subset of the sketch
table(crc5_sketch$sketch_snn_res.0.8)

#and sketch assay in the og object
table(crc5_filt$sketch_snn_res.0.8)

crc5_filt <- RunUMAP(crc5_filt,
                     reduction = "pca.sketch",
                     reduction.name = "umap.sketch",
                     return.model = TRUE,
                     dims = 1:30)

plot_hd24 <- DimPlot(crc5_filt, 
                     reduction = "umap.sketch",
                     label = TRUE,
                     group.by = "seurat_clusters")

print(plot_hd24)

###Project sketch data into the full dataset----
crc5_filt <- ProjectData(crc5_filt,
                         assay = "Spatial.008um",
                         full.reduction = "full.pca.sketch",
                         sketched.assay = "sketch",
                         sketched.reduction = "pca.sketch",
                         umap.model = "umap.sketch",
                         dims = 1:30,
                         refdata = list(seurat_cluster.projected = "seurat_clusters"))

DefaultAssay(crc5_filt) <- "Spatial.008um"

Idents(crc5_filt) <- "seurat_cluster.projected"

plot_hd25 <- DimPlot(crc5_filt,
                     reduction = "full.umap.sketch",
                     label = TRUE,
                     raster = FALSE)

plot_hd26 <- SpatialDimPlot(crc5_filt,
                            group.by = "seurat_cluster.projected",
                            pt.size.factor = 8,
                            label = TRUE,
                            repel = TRUE)

print(plot_hd25 | plot_hd26)

#Spatially-informed Clustering with Banksy in sketch data----
#Banksy is a method for clustering spatial omics data
#incorporating information of its spatial neighborhood
DefaultAssay(crc5_filt) <- "Spatial.008um"

crc5_filt <- RunBanksy(crc5_filt,
                       lambda = 0.8,
                       verbose = TRUE,
                       assay = "Spatial.008um",
                       slot = "data",
                       features = "variable",
                       k_geom = 50)

##Perform PCA on BANKSY assay----

DefaultAssay(crc5_filt) <- "BANKSY"

crc5_filt <- RunPCA(crc5_filt,
                    assay = "BANKSY",
                    reduction.name = "pca.banksy",
                    features = rownames(crc5_filt),
                    npcs = 50)

plot_hd27 <- ElbowPlot(crc5_filt, 
                       reduction = "pca.banksy", 
                       ndims = 50)

print(plot_hd27)

##Find neighbors and clusters----

crc5_filt <- FindNeighbors(crc5_filt,
                           reduction = "pca.banksy",
                           dims = 1:30)

crc5_filt <- FindClusters(crc5_filt,
                          cluster.name = "banksy_cluster",
                          resolution = 0.3,
                          algorithm = 4,
                          random.seed = 1234)

##Visualize BANKSY clusters----

crc5_filt <- RunUMAP(crc5_filt,
                     reduction = "pca.banksy",
                     reduction.name = "umap.banksy",
                     dims = 1:30)

Idents(crc5_filt) <- "banksy_cluster"

plot_hd28 <- DimPlot(crc5_filt,
                     reduction = "umap.banksy",
                     label = TRUE,
                     raster = FALSE)

plot_hd29 <- SpatialDimPlot(crc5_filt,
                            group.by = "banksy_cluster",
                            pt.size.factor = 8,
                            label = TRUE,
                            repel = TRUE)

print(plot_hd28 | plot_hd29)
#Compare the clusters with & without Banksy ¿What are the principal differences?

#Spatially-informed Clustering with Banksy in subset sketch data----

DefaultAssay(crc5_sketch) <- "Spatial.008um"

crc5_sketch <- RunBanksy(crc5_sketch,
                         lambda = 0.8,
                         verbose = TRUE,
                         assay = "Spatial.008um",
                         slot = "data",
                         features = "variable",
                         k_geom = 50)

##Perform PCA on BANKSY assay----

DefaultAssay(crc5_sketch) <- "BANKSY"

crc5_sketch <- RunPCA(crc5_sketch,
                      assay = "BANKSY",
                      reduction.name = "pca.banksy.sketch",
                      features = rownames(crc5_sketch),
                      npcs = 50)

plot_hd30 <- ElbowPlot(crc5_sketch, 
                       reduction = "pca.banksy.sketch", 
                       ndims = 50)

print(plot_hd30)

##Find neighbors and clusters----

crc5_sketch <- FindNeighbors(crc5_sketch,
                             reduction = "pca.banksy.sketch",
                             dims = 1:30)

crc5_sketch <- FindClusters(crc5_sketch,
                            cluster.name = "banksy_cluster",
                            resolution = 0.3,
                            algorithm = 4,
                            random.seed = 1234)

##Visualize BANKSY clusters----

crc5_sketch <- RunUMAP(crc5_sketch,
                       reduction = "pca.banksy.sketch",
                       reduction.name = "umap.banksy.sketch",
                       dims = 1:30)

Idents(crc5_sketch) <- "banksy_cluster"

plot_hd31 <- DimPlot(crc5_sketch,
                     reduction = "umap.banksy.sketch",
                     label = TRUE,
                     raster = FALSE)

plot_hd32 <- SpatialDimPlot(crc5_sketch,
                            group.by = "banksy_cluster",
                            pt.size.factor = 8,
                            label = TRUE,
                            repel = TRUE)

print(plot_hd31 | plot_hd32)

#Cell type annotation with RCTD----
#Robust Cell Type Deconvolution (RCTD) permits identify the cell profiles
#present in the scRNA-seq data to label the spatial transcriptomics pixels
#as cell types

##Prepare query data----
#The query data is the one that we would like to annotate (spatial)

DefaultAssay(crc5_filt) <- "sketch"
#We obatain the counts
counts_hd <- crc5_filt[["sketch"]]$counts
#The genes
cells_hd <- colnames(crc5_filt[["sketch"]])
#And also the coordinates of the bins
coords <- GetTissueCoordinates(crc5_filt)[cells_hd, 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

##Prepare reference data----
#The reference data is the one with our annotations of interest (scRNA-seq)
Idents(crc_integrated) <- "seurat_clusters"

counts <- crc_integrated[["RNA"]]$counts

cluster <- as.factor(crc_integrated$seurat_clusters)

nUMI <- crc_integrated$nCount_RNA

cluster <- droplevels(cluster)

reference <- Reference(counts, cluster, nUMI)

##Run RCTD----

RCTD <- create.RCTD(query, reference, max_cores = 4)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

#Append the results of the deconvolution to our query data
crc5_filt <- AddMetaData(crc5_filt, metadata = RCTD@results$results_df)

crc5_filt$first_type <- as.character(crc5_filt$first_type)

crc5_filt$first_type[is.na(crc5_filt$first_type)] <- "Unknown"

crc5_filt <- ProjectData(crc5_filt,
                         assay = "Spatial.008um",
                         full.reduction = "full.pca.sketch",
                         sketched.assay = "sketch",
                         sketched.reduction = "pca.sketch",
                         umap.model = "umap.sketch",
                         dims = 1:30,
                         refdata = list(full_first_type = "first_type"))

Idents(crc5_filt) <- "full_first_type"

SpatialDimPlot(crc5_filt, ncol = 4, pt.size.factor = 8)
#AQUI VA OTRO METODO DE DECONVOLUCION (SPOTLIGHT)----
#Agregarcomo se hacen estos diagramas de torta para los spots 
#recordar probarlo en el visium normal
#Differential expression analysis----
DefaultAssay(crc5_filt) <- "Spatial.008um"

Idents(crc5_filt) <- "seurat_clusters"

#We can change the parameters depending on what's the information
#of interest
hd_markers <- FindAllMarkers(crc5_filt,
                             only.pos = TRUE, 
                             min.pct = 0.25,
                             logfc.threshold = 0.25,
                             assay = "Spatial.008um")

write.csv(hd_markers, file = "hd_markers_crc5.csv")

#Analyze the Top 3 upregulated genes
hd_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

#Bubble o Dot Plot of the Top 3 genes
plot_hd33 <- DotPlot(crc5_filt,
                     features = unique(top3$gene),
                     group.by = "seurat_clusters",
                     dot.min = 0.25) +
  RotatedAxis()

print(plot_hd33)

#Violin Plot
plot_hd34 <- VlnPlot(crc5_filt,
                     features = unique(top3$gene),
                     group.by = "seurat_clusters",
                     pt.size = 0, ncol = 5, )

print(plot_hd34)

#Heatmap construction
crc5_heatmap <- ScaleData(crc5_filt, 
                          assay = "Spatial.008um", 
                          features = unique(top3$gene))

plot_hd35 <- DoHeatmap(crc5_heatmap,
                       features = unique(top3$gene),
                       group.by = "seurat_clusters", combine = TRUE)

print(plot_hd35)
