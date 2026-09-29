# ---------------------------------------------------------------------------
# Step 04 - Cluster markers, and what the CCA filter discarded.
#
# ESTABLISHES
#     The marker table behind the step 03 labels, and a characterisation of the
#     cells the variance-ratio filter would have removed.
#
# WHY THIS WAY
#     This step reads TWO objects: the kept object from step 02 and the
#     discarded object from step 02b. That is why the discarding variant is in
#     the repository at all - the markers of the discarded set are what justify
#     not discarding them.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 02, and step 02b for the discarded set
#
# INPUTS
#     <B_DIR>/Clustering_noDiscarding/<project_name>.<res>.<Dim>.<perp>.Seurat_object.rds
#     <B_DIR>/Clustering/discarded/Seurat_object.integrated.discard_<Dim>CC.rds
#
# OUTPUTS
#     marker tables and marker figures under <B_DIR>/MarkerIdentification
#
# USAGE
#     Rscript 04_*.R
# ---------------------------------------------------------------------------

# ===========================================================================
# 0. parameters
# ===========================================================================
# ---------------------------------------------------------------------------
# cluster locations. Point these at your own copies; nothing else changes.
# ---------------------------------------------------------------------------
SEQ_DIR  <- "/project/ctb-jsjoyal/gaelcge/Sequencing/Merging/TimeCourseOIR"
PROJ_DIR <- "/project/def-jsjoyal/gaelcge"
TC_DIR   <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse")
B_DIR    <- file.path(TC_DIR, "Aligned/Cond_Sorting")

# object identity, interpolated into output filenames by the original scripts.
# res=3 / 17 dims is this analysis; the P5-P17 analysis in the sibling folder
# used res=1 / 20 dims on a different object. Do not mix them up.
project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"
res  <- 3
Dim  <- 17
DIM_nb <- c(1:17)
perp <- 30

rm(list = ls())
library(plotly)
library(ggrepel)
library(gtools)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(useful)
library(Matrix)
library(Seurat)


#Set directory of dataset to analyse

setwd(file.path(B_DIR, "Clustering_noDiscarding/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.rds", sep="."))


#Set working directory

setwd(file.path(B_DIR, "MarkerIdentification"))


#Find cluster gene markers

Seurat_object.markers <- FindAllMarkers(object = Seurat_object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

Seurat_object.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file =paste(project_name, res, Dim, perp, "Markers-top10.tsv", sep="."), sep = "\t")

top10 <- read.table(file =paste(project_name, res, Dim, perp, "Markers-top10.tsv", sep="."))

png(filename="ClusterTree.png", width=1000, height=1000, bg = "white", res = 150)
Seurat_object <- BuildClusterTree(object = Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

#
#setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

png(filename=paste(project_name, res, Dim, perp, "HeatMap_cluster.png", sep="."), width=3000, height=3000, bg = "white", res = 100)
DoHeatmap(object =Seurat_object, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
#

markers.retina <- c("RHO", "LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "OTX2", "PRDM1", "RLBP1", "GFAP", 
					"CRABP1", "SLC5A7", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL", "SNHG11", "MEG3")

png(filename=paste(project_name, res, Dim, perp, "VlnPlot_markers.png", sep="."), width=1000, height=4000, bg = "white", res = 100)
VlnPlot(Seurat_object, markers.retina, group.by = "ident", point.size.use=-1, nCol = 1)
dev.off()



png(filename=paste(project_name, res, Dim, perp, "FeaturePlot_VEGFA.png", sep="."), width=1100, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, "VEGFA", min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("lightgrey", "blue"), pch.use = 16, overlay = FALSE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne", no.axes = FALSE, no.legend = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "DotPlot_VEGFA.png", sep="."), width=700, height=1000, bg = "white", res = 150)
DotPlot(Seurat_object, "VEGFA",
  cols.use = c("royalblue1", "red"))
dev.off()


###Subcluster prolif et non prolif VE

Seurat_object.markers <- FindAllMarkers(object = Seurat_object, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.5)

Seurat_object.markers %>% group_by(cluster) %>% top_n(100, avg_logFC) -> top100

write.table(top100, file =paste(project_name, res, Dim, perp, "Markers_VE_prolif-top100.tsv", sep="."), sep = "\t")




#Set directory of the discarded dataset to analyse

setwd(file.path(B_DIR, "Clustering/discarded"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste("Seurat_object.integrated.discard_", Dim, "CC.rds", sep=""))

#Set working directory

setwd(file.path(B_DIR, "MarkerIdentification/Discarded"))



#Find cluster gene markers

Seurat_object <- FindClusters(Seurat_object, reduction.type = "pca",
                                    dims.use = DIM_nb, save.SNN = T, resolution = 1.5, temp.file.location = getwd(), force.recalc=TRUE)

png(filename="UmapPlot-cluster.ident.discard.pca.png", width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename="UmapPlot-cluster.Cond_Sorting.discard.pca.png", width=1500, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = FALSE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename="UmapPlot-cluster.Batch.discard.pca.png", width=1500, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = FALSE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename="TSNEPlot-cluster.pc1-15.ident.discard.pca.png", width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename="TSNEPlot-cluster.pc1-15.Cond_Sorting.discard.pca.png", width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename="TSNEPlot-cluster.pc1-15.Batch.discard.pca.png", width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = FALSE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()


oldidentname <- c("8")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "6")
} 

oldidentname <- c("3", "4", "5", "12")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "0")
}  


png(filename="UmapPlot-cluster.ident.discard.regrouped.pca.png", width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()



Seurat_object.markers <- FindAllMarkers(object = Seurat_object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

Seurat_object.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file =paste(project_name, res, Dim, perp, "Markers-top10.tsv", sep="."), sep = "\t")

top10 <- read.table(file =paste(project_name, res, Dim, perp, "Markers-top10.tsv", sep="."))

png(filename="ClusterTree.png", width=1000, height=1000, bg = "white", res = 150)
Seurat_object <- BuildClusterTree(object = Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

#
#setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

png(filename=paste(project_name, res, Dim, perp, "HeatMap_cluster.png", sep="."), width=3000, height=3000, bg = "white", res = 100)
DoHeatmap(object =Seurat_object, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
#

markers.retina <- c("RHO", "LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "OTX2", "PRDM1", "RLBP1", "GFAP", 
					"CRABP1", "SLC5A7", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL", "SNHG11", "MEG3")

png(filename=paste(project_name, res, Dim, perp, "VlnPlot_markers.png", sep="."), width=1000, height=4000, bg = "white", res = 100)
VlnPlot(Seurat_object, markers.retina, group.by = "ident", point.size.use=-1, nCol = 1)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.initial.png", sep="."), res = 150, width=1000, height=2000)
DotPlot(Seurat_object, genes.plot= markers.retina, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "DotPlot_top10markers.renamed.initial.png", sep="."), res = 150, width=3000, height=2000)
DotPlot(Seurat_object, genes.plot= unique(top10$gene), cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


#



#
#


q("no")
