# ---------------------------------------------------------------------------
# Step 03 - Cell-type labels assigned to the resolution-3 clusters.
#
# ESTABLISHES
#     The annotated object that every remaining step in this folder reads.
#
# WHY THIS WAY
#     Labels are assigned by hand from marker expression over an intentionally
#     over-partitioned clustering, so that several clusters collapse onto one
#     cell type and no cluster carries two. The alternative - clustering at a
#     resolution that yields one cluster per expected type - forces a
#     resolution choice to stand in for a biological judgement.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 02
#
# INPUTS
#     <B_DIR>/Clustering_noDiscarding/<project_name>.<res>.<Dim>.<perp>.Seurat_object.rds
#
# OUTPUTS
#     <B_DIR>/CellTypeAnnotation/<project_name>.<res>.<Dim>.<perp>.Seurat_object.annotated.rds
#
# USAGE
#     Rscript 03_*.R
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

GENE_LIST_DIR      <- file.path(PROJ_DIR, "Gene_list")
GMT_DIR            <- file.path(GENE_LIST_DIR, "Gmt.file")
# MSigDB v6.2 here, NOT the v7.1 used by 1_timecourse_p5_p17. Both releases are
# recorded in reference/PROVENANCE.txt; do not substitute one for the other.
# Three curated senescence .gmx/.gmt files are held by a collaborator and are
# NOT readable by this account - see reference/PROVENANCE.txt.
CUSTOM_GENESET_DIR <- file.path(PROJ_DIR, "..", "jhowa105/projects/Mike")

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
library(tidyr)
library(ggrepel)
library(gtools)
library(data.table)
library(gplots)
library(useful)
library(locfit)
library(Matrix)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Rmagic)
library(readr)
library(phateR)
library(viridis)
library(doubletFinder)

#Set directory of dataset to analyse

setwd(file.path(B_DIR, "Clustering_noDiscarding/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.rds", sep="."))


#Set working directory

setwd(file.path(B_DIR, "CellTypeAnnotation"))

summary(Seurat_object@meta.data)


#rename cluster ID

markers.retina.dotplot <- rev(c("RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "RLBP1", "FBN1", "CLDN5", "LYZ2", "GFAP", "OPTC", "TOP2A", "MKI67"))

png(filename=paste(project_name, res, Dim, perp, "ClusterTree_renamed.initial.png", sep="."), res = 150, width=1000, height=1500)
Seurat_object <- BuildClusterTree(Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.initial.png", sep="."), res = 150, width=1000, height=2000)
DotPlot(Seurat_object, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()



oldidentname <- c("0", "4", "9", "10", "11", "16", "17", "20", "26")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Bipolar cells")
} 



oldidentname <- c("1", "2", "3", "5", "6", "12", "21", "38")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Rods")
} 


oldidentname <- c("8", "32")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Cones")
} 

#
oldidentname <- c("7", "13", "22", "27")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Amacrine cells")
} 


oldidentname <- c("14", "15", "18", "19", "24", "25", "30", "33")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Muller glia")
} 


oldidentname <- c("34", "39", "40")
for (clusterID in oldidentname) {

    Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Immune cells")
} 

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "23", new.ident.name = "Retinal ganglion cells")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "29", new.ident.name = "Endothelial cells")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "31", new.ident.name = "Horizontal cells")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "36", new.ident.name = "Astrocytes")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "28", new.ident.name = "Pericytes")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "35", new.ident.name = "Basal cells 1")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "37", new.ident.name = "Basal cells 2")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = "41", new.ident.name = "Red blood cells")


###

cell_type <- data.frame("cell_type"=Seurat_object@ident, stringsAsFactors=F)

Seurat_object <- AddMetaData(object = Seurat_object, metadata = cell_type)

Seurat_object <- ReorderIdent(Seurat_object, feature = c("CC1", "CC2"), rev = FALSE, aggregate.fxn = mean,
  reorder.numeric = FALSE)


#Cluster proportion

cluster_counts <- Seurat_object@meta.data %>%
  group_by(Sorting) %>%
  count(cell_type) %>%
  mutate(freq = n / sum(n))

cluster_counts$cell_type <- factor(cluster_counts$cell_type,levels=runique(mixedsort(cluster_counts$cell_type)))


png(filename=paste(project_name, res, Dim, perp, "Cluster_proportions_cell_type_renamed.Sorting.png", sep="."), width=1500, height=1000, res=150)
ggplot(cluster_counts,aes(x=reorder(Sorting, cell_type),y=freq,fill=Sorting)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ cell_type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

cluster_counts <- cluster_counts %>% group_by(cell_type) %>% mutate(freq_sum = sum(freq))

png(filename=paste(project_name, res, Dim, perp, "Cluster_proportions_cell_type.full.renamed.Sorting.png", sep="."), width=600, height=600, res=150)
ggplot(cluster_counts,aes(x=reorder(cell_type, -freq_sum),y=freq,fill=Sorting)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cell type ID') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

cluster_counts <- Seurat_object@meta.data %>%
  group_by(Condition) %>%
  count(cell_type) %>%
  mutate(freq = n / sum(n))

cluster_counts$cell_type <- factor(cluster_counts$cell_type,levels=unique(mixedsort(cluster_counts$cell_type)))


png(filename=paste(project_name, res, Dim, perp, "Cluster_proportions_cell_type_renamed.Condition.png", sep="."), width=1500, height=1000, res=150)
ggplot(cluster_counts,aes(x=reorder(Condition, cell_type),y=freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ cell_type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

cluster_counts <- cluster_counts %>% group_by(cell_type) %>% mutate(freq_sum = sum(freq))

png(filename=paste(project_name, res, Dim, perp, "Cluster_proportions_cell_type.full.renamed.Condition.png", sep="."), width=600, height=600, res=150)
ggplot(cluster_counts,aes(x=reorder(cell_type, -freq_sum),y=freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cell type ID') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

cluster_ratio <- Seurat_object@meta.data %>%
  group_by(cell_type) %>%
  count(Sorting) %>%
  mutate(ratio = n / sum(n))

png(filename=paste(project_name, res, Dim, perp, "Cluster_ratio_cell_type.full.renamed.Sorting.png", sep="."), width=1200, height=700, res = 150)
ggplot(cluster_ratio,aes(reorder(x=reorder(cell_type, -ratio), ratio),ratio,fill=Sorting)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cell type ID') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

cluster_ratio <- Seurat_object@meta.data %>%
  group_by(cell_type) %>%
  count(Condition) %>%
  mutate(ratio = n / sum(n))

png(filename=paste(project_name, res, Dim, perp, "Cluster_ratio_cell_type.full.renamed.Condition.png", sep="."), width=1200, height=700, res = 150)
ggplot(cluster_ratio,aes(reorder(x=reorder(cell_type, -ratio), ratio),ratio,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cell type ID') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



bp <- ggplot(cluster_ratio,aes(reorder(cell_type, ratio),ratio,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cell type ID')


png(filename=paste(project_name, res, Dim, perp, "RelativeProportion_ID_condition_retina.png", sep="."), width=1500, height=1500, res = 150)
bp + coord_polar("x", start=0) + theme(
   panel.background = element_rect(fill = "transparent", colour = NA),
   plot.background = element_rect(fill = "transparent", colour = NA)
 )
dev.off()


png(filename=paste(project_name, res, Dim, perp, "Cell_count.png", sep="."))
par(mar=c(5, 15 ,4 ,2))
barplot(sort(summary(Seurat_object@ident)), horiz=TRUE, las=2, col=c("limegreen"), log = "x", xlab = "Cell count", ann=FALSE)
dev.off()

###add cell cycle

cc.genes <- readLines(con=CC_GENES)

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]

Seurat_object <- CellCycleScoring(Seurat_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)

phase_counts <- Seurat_object@meta.data %>%
  group_by(cell_type) %>%
  count(Phase) %>%
  mutate(freq = n / sum(n))  

phase_counts$Phase  = factor(phase_counts$Phase, levels=c("G1", "S", "G2M"))

png(filename=paste(project_name, res, Dim, perp, "CellCycle_proportions_cell_type.png", sep="."), width=1200, height=1000, res = 150) 
ggplot(phase_counts,aes(x=Phase,freq,fill=Phase)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ cell_type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


##Making final plots

set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = 1)(length(identities)))

#hue_pal()
#Arguments
#h = range of hues to use, in [0, 360]
#c =chroma (intensity of colour), maximum value varies depending on combination of hue and luminance.
#l =luminance (lightness), in [0, 100]
#h.start = hue to start at
#direction = direction to travel around the colour wheel, 1 = clockwise, -1 = counter-clockwise


png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.noLabel.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 8, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.noLabel.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 8, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.Label.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.Label.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.TimePoint_legend.png", sep="."), width=1000, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 5, no.legend = FALSE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.TimePoint_legend.png", sep="."), width=1000, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 5, no.legend = FALSE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()


png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.Condition_legend.png", sep="."), width=1000, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Condition", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 5, no.legend = FALSE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.Condition_legend.png", sep="."), width=1000, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Condition", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 5, no.legend = FALSE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()


png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.Sorting_legend.png", sep="."), width=1000, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 5, no.legend = FALSE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "UMAPPlot.renamed.final.Sorting_legend.png", sep="."), width=1000, height=800, bg = "white", res = 100)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 5, no.legend = FALSE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()


markers.retina <- c("RHO", "OPN1MW", "TRPM1", "SNHG11", "MEG3", "CLDN5", "PECAM1", "RLBP1", "KCNJ8", "LYZ2", "GFAP", "GRIA2", "OPTC")

png(filename=paste(project_name, res, Dim, perp, "VlnPlot_markers.renamed.final.png", sep="."), width=800, height=2000, bg = "white", res = 50)
VlnPlot(Seurat_object, markers.retina, ident.include = NULL, nCol = 1,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = "ident", y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE,
  return.plotlist = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "VlnPlot_markers.renamed.imputed.final.png", sep="."), width=800, height=7000, bg = "white", res = 50)
VlnPlot(Seurat_object, markers.retina, ident.include = NULL, nCol = 1,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = "ident", y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE,
  return.plotlist = FALSE, use.imputed=TRUE)
dev.off()


markers.retina.dotplot <- c("RHO", "OPN1SW", "TRPM1", "RLBP1", "OPTC", "KCNJ8", "LYZ2", "CLDN5")

png(filename=paste(project_name, res, Dim, perp, "ClusterTree_renamed.final.png", sep="."), res = 150, width=1000, height=1500)
Seurat_object <- BuildClusterTree(Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.final.png", sep="."), res = 150, width=1000, height=800)
DotPlot(Seurat_object, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 10,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

saveRDS(Seurat_object, paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))



lfile <- as.loom(Seurat_object, assay = NULL, filename = "Seurat_object.new.1.loom", verbose = TRUE, overwrite = TRUE, max.size = "4000mb")


lfile <- as.loom(Seurat_object, assay = NULL, filename = "Seurat_object.new.loom", max.size = "400mb",
  chunk.dims = NULL, chunk.size = NULL, overwrite = FALSE,
  verbose = TRUE)


Convert(Seurat_object, "loom", "Seurat_object.OIR_TimeCourse1.loom", chunk.dims = "auto",
  chunk.size = 1000, overwrite = TRUE, display.progress = TRUE,
  anndata.raw = "raw.data", anndata.X = "data")


###Continue working with the seurat object for more annotation

setwd(file.path(B_DIR, "CellTypeAnnotation"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))


####
# Seurat_object <- UpdateSeuratObject(object = Seurat_object)

# Seurat_object <- NormalizeData(object = Seurat_object)

# Seurat_object <- FindVariableFeatures(Seurat_object)

# Seurat_object <- ScaleData(object = Seurat_object)
# Seurat_object <- RunPCA(object = Seurat_object)
# Seurat_object <- FindNeighbors(object = Seurat_object)
# Seurat_object <- FindClusters(object = Seurat_object)
# Seurat_object <- RunTSNE(object = Seurat_object)

# Subset time points
setwd(file.path(B_DIR, "CellTypeAnnotation/P14_P17"))

Subseted_cells <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14", "P17")))

Seurat_object.Subseted <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)


###Make processed matrix for GEO submission without Activated Muller glia

Seurat_object.Subseted <- SubsetData(Seurat_object.Subseted, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Activated Muller Glia"), accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

df <- as.matrix(Seurat_object.Subseted@data)

head(colnames(df))

head(rownames(df))

write.table(df, "retina_NORM_OIR_P14_P17_C57_WR_CD73FT_noamg_normalizedUMI_Count_DGEmatrix.txt", sep="\t",
                                      quote=FALSE,
                                      row.names=TRUE,
                                      col.names=TRUE)


###Write meta.data

df <- as.matrix(Seurat_object.Subseted@meta.data)

head(colnames(df))

head(rownames(df))

write.table(df, "retina_NORM_OIR_P14_P17_C57_WR_CD73FT_noamg_meta.data.txt", sep="\t",
                                      quote=FALSE,
                                      row.names=TRUE,
                                      col.names=TRUE)
###Renamed and removed unwanted IDs

Seurat_object.Subseted <- RenameIdent(Seurat_object.Subseted, old.ident.name = c("Quiescent Muller Glia"), new.ident.name = "Muller glia")
Seurat_object.Subseted <- RenameIdent(Seurat_object.Subseted, old.ident.name = c("Activated Muller Glia"), new.ident.name = "Muller glia")

Seurat_object.Subseted <- SubsetData(Seurat_object.Subseted, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Basal cells 1", "Basal cells 2", "Red blood cells", "Mki67 Muller Glia", "Early Muller Glia"), accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


summary(Seurat_object.Subseted@ident)

### Add new meatdata

df <- as.data.frame(paste(Seurat_object.Subseted@meta.data$Cond_Sorting, Seurat_object.Subseted@meta.data$TimePoint, sep="_"))
rownames(df) <- rownames(Seurat_object.Subseted@meta.data)
colnames(df) <- "Cond_Sorting_TimePoint"

Seurat_object.Subseted <- AddMetaData(object =Seurat_object.Subseted, metadata = df, col.name = "Cond_Sorting_TimePoint")

df <- as.data.frame(paste(Seurat_object.Subseted@meta.data$Phase, Seurat_object.Subseted@meta.data$Condition, sep="_"))
rownames(df) <- rownames(Seurat_object.Subseted@meta.data)
colnames(df) <- "Phase_Condition"

Seurat_object.Subseted <- AddMetaData(object =Seurat_object.Subseted, metadata = df, col.name = "Phase_Condition")

df <- as.data.frame(paste(Seurat_object.Subseted@meta.data$cell_type, Seurat_object.Subseted@meta.data$Condition, sep="_"))
rownames(df) <- rownames(Seurat_object.Subseted@meta.data)
colnames(df) <- "cell_type_Condition"

Seurat_object.Subseted <- AddMetaData(object =Seurat_object.Subseted, metadata = df, col.name = "cell_type_Condition")

##3Re-cluster
Seurat_object.Subseted <- RunTSNE(Seurat_object.Subseted,
                               reduction.use = "cca.aligned",
                               dims.use = DIM_nb, do.fast = T, dim_embed=2, perplexity=perp)



set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = 1)(length(identities)))

#hue_pal()
#Arguments
#h = range of hues to use, in [0, 360]
#c =chroma (intensity of colour), maximum value varies depending on combination of hue and luminance.
#l =luminance (lightness), in [0, 100]
#h.start = hue to start at
#direction = direction to travel around the colour wheel, 1 = clockwise, -1 = counter-clockwise


png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.noLabel.P17.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object.Subseted, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 8, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

svg(paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.noLabel.P17.svg", sep="."), width=8, height=8, bg = "white")
DimPlot(Seurat_object.Subseted, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 8, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.noLabel.P14_P17.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object.Subseted, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 8, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.Label.P17.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object.Subseted, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

svg(paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.Label.P17.svg", sep="."), width=8, height=8, bg = "white")
DimPlot(Seurat_object.Subseted, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 8, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.Label.P14_P17.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object.Subseted, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()


Seurat_object.Subseted_immunecells <- SubsetData(Seurat_object.Subseted, cells.use = NULL, subset.name =NULL, ident.use = "Immune cells",
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


png(filename=paste(project_name, res, Dim, perp,  "TSNEPlot.renamed.final.Label.P14_P17_immunecells.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object.Subseted_immunecells, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()

png(filename=paste(project_name, res, Dim, perp,  "UMAPPlot.renamed.final.Label.P14_P17_immunecells.png", sep="."), width=800, height=800, bg = "white", res = 100)
DimPlot(Seurat_object.Subseted_immunecells, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = TRUE,
  dark.theme = FALSE)
dev.off()


### Cell cycle annotation

phase_counts <- Seurat_object.Subseted@meta.data %>%
  group_by(cell_type_Condition) %>%
  count(Phase_Condition) %>%
  mutate(freq = n / sum(n))  

phase_counts$Phase_Condition  = factor(phase_counts$Phase_Condition, levels=c("G1_NORM", "G1_OIR", "S_NORM", "S_OIR", "G2M_NORM", "G2M_OIR"))

phase_counts <- phase_counts %>%
  separate(cell_type_Condition,c("cell_type", "Condition"),sep="_")


png(filename=paste(project_name, res, Dim, perp, "CellCycle_proportions_cell_type_Phase_Condition_P17.png", sep="."), width=1200, height=1000, res = 150) 
ggplot(phase_counts,aes(x=Phase_Condition,freq,fill=Phase_Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ cell_type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette="YlGnBu", direction=1)
dev.off()


phase_counts <- Seurat_object.Subseted@meta.data %>%
  group_by(Condition) %>%
  count(Phase_Condition) %>%
  mutate(freq = n / sum(n))  

phase_counts$Phase_Condition  = factor(phase_counts$Phase_Condition, levels=c("G1_NORM", "G1_OIR", "S_NORM", "S_OIR", "G2M_NORM", "G2M_OIR"))

png(filename=paste(project_name, res, Dim, perp, "CellCycle_proportions_Retina_Phase_Condition_P17.png", sep="."), width=1200, height=1000, res = 150) 
ggplot(phase_counts,aes(x=Phase_Condition,freq,fill=Phase_Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette="YlGnBu", direction=1)
dev.off()


q("no")
