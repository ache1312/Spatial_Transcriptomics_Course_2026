#!/usr/bin/env Rscript

############################################################
# Bioinformatics Analysis of Spatial Biology Data
# Xenium Workshop - Cell Chat
############################################################

# -------------------- Setup ------------------------------

set.seed(1234)
library(Seurat)
library(ggplot2)
library(CellChat)

plan("multisession", workers = 4)
options(future.globals.maxSize = 8 * 1024^3)


# -------------------- Paths ------------------------------

qc_dir <- "/home/courses/student75/xenium_plots/cellchat"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("Cellchat plots will be saved in:", qc_dir, "\n")

# -------------------- Data loading -----------------------

xen_bc_umap <- readRDS(file = "/home/courses/student75/xen_bc_umap.rds")

DefaultAssay(xen_bc_umap) = "Xenium"
Idents(xen_bc_umap) = "first_type"

# -------------------- Prepare data -----------------------

data_input <- GetAssayData(xen_bc_umap, assay = "Xenium", layer = "data")
spatial_locs <- GetTissueCoordinates(xen_bc_umap, cols = c("imagerow", "imagecol"))
spatial_locs <- spatial_locs[,1:2]

meta = data.frame(labels = Idents(xen_bc_umap), 
                  samples = "Xen_bc", 
                  row.names = names(Idents(xen_bc_umap)))

spot_size <- 10
conversion_factor <- 1
spatial_factors <- data.frame(ratio = conversion_factor,
                              tol = spot_size/2)

# ---------------------- Create CellChat object -------------

cellchat <- createCellChat(object = data_input, 
                           meta = meta, 
                           group.by = "labels",
                           datatype = "spatial", 
                           coordinates = spatial_locs, 
                           spatial.factors = spatial_factors)

CellChatDB <- CellChatDB.human
glimpse(CellChatDB$interaction)
CellChatDB_sub <- subsetDB(CellChatDB)

cellchat@DB <- CellChatDB_sub
cellchat <- subsetData(cellchat)

# ----------------------- Identify genes and interactions -----

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, 
                                              variable.both = F)

# ----------------------- Calculate communication probability -

cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean", 
                              trim = 0.1,
                              distance.use = TRUE, 
                              interaction.range = 250, 
                              scale.distance = 4,
                              contact.dependent = TRUE, 
                              contact.range = 10)
cellchat <- filterCommunication(cellchat, 
                                min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

cat("Check patways:")
cellchat@netP$pathways

cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP")

# --------------------- Plots ----------------------------------

jpeg(file.path(qc_dir, "plot53_net_circle_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netVisual_circle(cellchat@net$count,
                 vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions in CD8 cells",
                 targets.use = "CD8")
dev.off()

jpeg(file.path(qc_dir, "plot54_net_heatmap_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netVisual_heatmap(cellchat, 
                  measure = "count", 
                  color.heatmap = "Blues")
dev.off()

jpeg(file.path(qc_dir, "plot55_net_aggregateCXCL_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netVisual_aggregate(cellchat, 
                    signaling = "CXCL", 
                    layout = "spatial", 
                    edge.width.max = 2, 
                    vertex.size.max = 0.5, 
                    alpha.image = 0.01, 
                    vertex.label.cex = 3.5,sources.use = "Epithelial",	
                    targets.use = "CD8", remove.isolate = TRUE)
dev.off()

jpeg(file.path(qc_dir, "plot56_net_aggregateIL16_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netVisual_aggregate(cellchat, 
                    signaling = "IL16", 
                    layout = "spatial", 
                    edge.width.max = 2, 
                    vertex.size.max = 0.5, 
                    alpha.image = 0.01, 
                    vertex.label.cex = 3.5,
                    targets.use = "Mac",
                    sources.use = "B",
                    remove.isolate = TRUE)
dev.off()

jpeg(file.path(qc_dir, "plot57_net_aggregateCD86_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netVisual_aggregate(cellchat, 
                    signaling = "CD86", 
                    layout = "spatial", 
                    edge.width.max = 2, 
                    vertex.size.max = 0.5, 
                    alpha.image = 0.01, 
                    vertex.label.cex = 3.5,
                    targets.use = "CD4",
                    sources.use = "DC",
                    remove.isolate = TRUE)
dev.off()

jpeg(file.path(qc_dir, "plot58_net_aggregateCD39_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netVisual_aggregate(cellchat, 
                    signaling = "CD39", 
                    layout = "spatial", 
                    edge.width.max = 2, 
                    vertex.size.max = 0.5, 
                    alpha.image = 0.01, 
                    vertex.label.cex = 3.5,
                    targets.use = "Epithelial",
                    sources.use = "TEC",
                    remove.isolate = TRUE)
dev.off()

jpeg(file.path(qc_dir, "plot59_net_signalroleCXCL_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = "CXCL", 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)
dev.off()

jpeg(file.path(qc_dir, "plot60_net_signalroleIL16_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = "IL16", 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)
dev.off()

jpeg(file.path(qc_dir, "plot61_net_signalroleCD86_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = "CD86", 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)
dev.off()

jpeg(file.path(qc_dir, "plot62_net_signalroleCD39_Xenium.jpeg"), width = 6, height = 4,
     res = 300, units = "in")
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = "CD39", 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)
dev.off()

cat("Done. Cellchat plots saved in:", qc_dir, "\n")


                             
