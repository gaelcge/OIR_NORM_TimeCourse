# ---------------------------------------------------------------------------
# 04_pericytes - 01_subcluster.R
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# pericytes, run on cells taken from the annotated
# object of that folder's step 03.
#
# WHICH OBJECT
#     Cells are subset from
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30.Seurat_object.annotated.rds
#     - P14/P17 only, Seurat v2 lineage. NOT the P5-P17 object of the sibling
#     folder.
#
# WHY SUBCLUSTER PER COMPARTMENT
#     A single clustering resolution cannot resolve microglial activation states
#     and rod subtypes at once: the resolution that splits the former shatters
#     the latter. Each compartment is therefore re-embedded and re-clustered on
#     its own cells, and ../11_subcluster_convergence.R is where the independent
#     label sets are reconciled onto one object.
#
# RUN ORDER
#     After ../03_annotate_celltypes.R, before ../11_subcluster_convergence.R.
#     Within this folder, in numbered order.
#
# RELIC CODE
#     Blocks commented out with `#~` were carried in from other projects and are
#     disabled, not deleted - they operate on objects this script never loads,
#     or write into directories that no longer exist. See ../../EDIT_POLICY.md.
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
library(magrittr)
library(harmony)


setwd(file.path(B_DIR, "Clustering/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

SeuratObject <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.rds", sep="."))

summary(as.factor(SeuratObject@meta.data$Cond_Sorting))


####Subclustering

setwd(file.path(B_DIR, "Subclustering/Pericytes/Mapping"))


markers.retina.dotplot <- rev(c("RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "RLBP1", "FBN1", "CLDN5", "LYZ2", "GFAP", "OPTC", "TOP2A", "MKI67"))

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.initial.png", sep="."), res = 150, width=1000, height=2000)
DotPlot(SeuratObject, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


SeuratObject <- SubsetData(SeuratObject, cells.use = NULL, subset.name =NULL, ident.use = c(33,38),
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


set.seed(001)

require(scales)

identities <- levels(SeuratObject@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="UMAPPlot-Subclustered_pericytes.png", width=500, height=500, bg = "white", res = 150)
DimPlot(SeuratObject, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

pericytes <- SeuratObject

pericytes <- ScaleData(object = pericytes, vars.to.regress = c("nUMI", "percent.mito", "Batch", "percent.crystal"))

png(filename="MeanVarPlot.png")
pericytes <- FindVariableGenes(object = pericytes, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.20, x.high.cutoff = 4, y.cutoff = 0.5)
dev.off()

length(x=pericytes@var.genes)


##DImension reduction

## PCA and ICA

pericytes <- RunPCA(object = pericytes, pc.genes = pericytes@var.genes, do.print = TRUE, pcs.print = 1:10, 
    genes.print = 10)
    
PrintPCA(object = pericytes, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

pericytes <- ProjectPCA(object = pericytes, do.print = FALSE)

pericytes <- RunICA(object = pericytes, ic.genes = pericytes@var.genes)

pericytes <- ProjectDim(pericytes, reduction.type = "ica", dims.print = 1:5,
  dims.store = 30, genes.print = 30, replace.dim = FALSE,
  do.center = FALSE, do.print = TRUE, assay.type = "RNA")


png(filename="VizPCA.png", width=800, height=2000, bg = "white", res = 150)
VizPCA(object =pericytes, pcs.use = 1:15)
dev.off()

png(filename="VizICA.png", width=800, height=2000, bg = "white", res = 150)
VizICA(object =pericytes, ics.use = 1:15)
dev.off()

png(filename="PCAPlot.batch.pc1vs2.pca.png")
DimPlot(object =pericytes, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()

png(filename="ICAPlot.batch.pc1vs2.pca.png")
DimPlot(object =pericytes, reduction.use = "ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()

png(filename="PCAPlot.batch.pc1vs2.pca.Cond_Sorting.png")
DimPlot(object =pericytes, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()

png(filename="ICAPlot.batch.pc1vs2.pca.Cond_Sorting.png")
DimPlot(object =pericytes, reduction.use = "ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()


png(filename="PCHeatmap_mult.pca.png", width=2000, height=4000, bg = "white", res = 150)
DimHeatmap(pericytes, assay.use = "RNA", reduction.type = "pca", dim.use = 1:12,
  cells.use = NULL, num.genes = 30, use.full = FALSE, disp.min = -2.5,
  disp.max = 2.5, do.return = FALSE, col.use = PurpleAndYellow(),
  use.scale = TRUE, do.balanced = FALSE, remove.key = FALSE,
  label.columns = NULL, check.plot = TRUE)
dev.off()

png(filename="ICHeatmap_mult.ica.png", width=2000, height=4000, bg = "white", res = 150)
DimHeatmap(pericytes, assay.use = "RNA", reduction.type = "ica", dim.use = 1:12,
  cells.use = 100, num.genes = 30, use.full = FALSE, disp.min = -2.5,
  disp.max = 2.5, do.return = FALSE, col.use = PurpleAndYellow(),
  use.scale = TRUE, do.balanced = FALSE, remove.key = FALSE,
  label.columns = NULL, check.plot = TRUE)
dev.off()

png(filename="PCElbowPlot.png")
PCElbowPlot(object = pericytes)
dev.off()

png(filename="ICAElbowPlot.png")
DimElbowPlot(pericytes, reduction.type = "ica", dims.plot = 20)
dev.off()


PrintDim(object = pericytes, reduction.type = "pca", dims.print = 1:2, genes.print = 10)

p1 <- DimPlot(object = pericytes, reduction.use = "pca", group.by = "Batch", pt.size = 1, 
    do.return = TRUE)
p2 <- VlnPlot(object = pericytes, features.plot = "PC1", group.by = "Batch", do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_pca.Batch.png", width=1500, height=1000, bg = "white", res = 100)
plot_grid(p1, p2)
dev.off()

PrintDim(object = pericytes, reduction.type = "ica", dims.print = 1:2, genes.print = 10)

p1 <- DimPlot(object = pericytes, reduction.use = "ica", group.by = "Batch", pt.size = 1, 
    do.return = TRUE)
p2 <- VlnPlot(object = pericytes, features.plot = "IC1", group.by = "Batch", do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_ica.Batch.png", width=1500, height=1000, bg = "white", res = 100)
plot_grid(p1, p2)
dev.off()


###### RUn harmony

pericytes <- RunHarmony(pericytes, "Sorting", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100)

pericytes <- ProjectDim(pericytes, reduction.type = "harmony", dims.print = 1:5,
  dims.store = 30, genes.print = 30, replace.dim = FALSE,
  do.center = FALSE, do.print = TRUE, assay.type = "RNA")


png(filename="VizHarmony.png", width=800, height=2000, bg = "white", res = 150)
VizDimReduction(pericytes, reduction.type = "harmony", dims.use = 1:15,
  num.genes = 30, use.full = FALSE, font.size = 0.5, nCol = NULL,
  do.balanced = FALSE)
dev.off()


png(filename="PCAPlot.batch.pc1vs2.harmony.png")
DimPlot(object =pericytes, reduction.use = "harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()


png(filename="PCHeatmap_mult.png", width=2000, height=4000, bg = "white", res = 150)
PCHeatmap(object = pericytes, pc.use = 1:12, cells.use = 100, do.balanced = TRUE, 
    label.columns = FALSE, use.full = FALSE)
dev.off()

png(filename="HarmonyHeatmap_mult.png", width=2000, height=4000, bg = "white", res = 150)
DimHeatmap(pericytes, assay.use = "RNA", reduction.type = "harmony", dim.use = 1:12,
  cells.use = 100, num.genes = 30, use.full = FALSE, disp.min = -2.5,
  disp.max = 2.5, do.return = FALSE, col.use = PurpleAndYellow(),
  use.scale = TRUE, do.balanced = FALSE, remove.key = FALSE,
  label.columns = NULL, check.plot = TRUE)
dev.off()

# No standard deviation info stored for harmony, so elbow on pca
png(filename="PCAElbowPlot.pca_harmony.png")
DimElbowPlot(pericytes, reduction.type = "pca", dims.plot = 20)
dev.off()


PrintDim(object = pericytes, reduction.type = "harmony", dims.print = 1:2, genes.print = 10)

p1 <- DimPlot(object = pericytes, reduction.use = "harmony", group.by = "Batch", pt.size = 1, 
    do.return = TRUE)
p2 <- VlnPlot(object = pericytes, features.plot = "Harmony1", group.by = "Batch", do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_Harmony.Batch.png", width=1500, height=1000, bg = "white", res = 100)
plot_grid(p1, p2)
dev.off()


##Rn TSNE

pericytes <- RunTSNE(object = pericytes, dims.use = 1:4, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.type = "pca", reduction.name = "tsne_pca")

pericytes <- RunTSNE(object = pericytes, dims.use = 1:4, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.type = "ica", reduction.name = "tsne_ica")

pericytes <- RunTSNE(object = pericytes, dims.use = 1:4, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.type = "harmony", reduction.name = "tsne_harmony")


##Run UMAP

pericytes <- RunUMAP(pericytes, cells.use = NULL, dims.use = 1:4, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_pca", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)

pericytes <- RunUMAP(pericytes, cells.use = NULL, dims.use = 1:4, reduction.use = "ica",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_ica", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)

pericytes <- RunUMAP(pericytes, cells.use = NULL, dims.use = 1:4, reduction.use = "harmony",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_harmony", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)


#######Find cluster

##For PCA
pericytes <- FindClusters(object= pericytes , reduction.type = "pca", dims.use = c(1:4), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN=1/8)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(pericytes@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="TSNEPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="UMAPPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


##For ICA

pericytes <- FindClusters(object= pericytes , reduction.type = "ica", dims.use = c(1:4), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN=1/8)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(pericytes@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))


png(filename="TSNEPlot-cluster.ica.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.ica.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="UMAPPlot-cluster.ica.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "umap_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.ica.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "umap_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



##For harmony

pericytes <- FindClusters(object= pericytes , reduction.type = "harmony", dims.use = c(1:4), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN=1/8)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(pericytes@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))


png(filename="TSNEPlot-cluster.harmony.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.harmony.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="UMAPPlot-cluster.harmony.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "umap_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.harmony.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "umap_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

### Choose dimensionality reduction for cluster finding

pericytes <- FindClusters(object= pericytes , reduction.type = "pca", dims.use = c(1:4), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN = 1/8)


##### QC on clusters

png(filename="VlnPlotQC.pericytes.pca.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(pericytes, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

png(filename="Tree_initial.png", res = 150, width=1000, height=1500)
BuildClusterTree(pericytes, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()


##Find cell types

setwd(file.path(B_DIR, "Subclustering/Pericytes/CellIdentification"))


## Regrouping

oldidentname <- c("2", "3", "4")
for (clusterID in oldidentname) {
  pericytes <- RenameIdent(pericytes, old.ident.name = clusterID, new.ident.name = "1")
} 

pericytes <-ValidateSpecificClusters(pericytes, cluster1 = 1, cluster2 = 8,
  pc.use = 1:5, top.genes = 30, acc.cutoff = 0.9)

png(filename="TSNEPlot-cluster.pca.regrouping.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="Tree_regrouping.png", res = 150, width=1000, height=1500)
BuildClusterTree(pericytes, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()

png(filename="VlnPlotQC.pericytes.harmony.regrouping.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(pericytes, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

#####

pericytes.markers <- FindAllMarkers(object = pericytes, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

pericytes.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file ="pericytes-subclustered-Markers-top10.tsv", sep = "\t")

top10.markers <- unique(top10$gene)

png(filename="DotPlot_pericytes_subclustered_top10_markers.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(pericytes, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="ClusterTree.png", width=1000, height=1000, bg = "white", res = 150)
pericytes <- BuildClusterTree(object = pericytes, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

pericytes.markersNode <- FindAllMarkersNode(object =pericytes, node = NULL, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)

pericytes.markersNode %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10.node

write.table(top10.node, file ="pericytes-subclustered-MarkersNode-top10.avg_logFC.tsv", sep = "\t")



###Gene set for cluster identification

neutro_gene <- read.csv("Neutrophil_gene_set.txt", header=F)

neutro_gene <-apply(neutro_gene, 1, toupper)

neutro_gene <- intersect(neutro_gene, rownames(pericytes@data))


png(filename="DotPlot_pericytes_subclustered_neutro_genes.png", width=1000, height=1000, bg = "white", res = 150)
DotPlot(pericytes, neutro_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="FeaturePlot_pericytes_subclustered_TLR2_CD14.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(pericytes, c("TLR2", "CD14"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 2,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_pericytes_subclustered_ICAM1_CD63.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(pericytes, c("ICAM1", "CD63"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 2,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()


Macrophage_gene <- read.csv("Macrophage_gene_set.txt", header=F)

Macrophage_gene <-apply(Macrophage_gene, 1, toupper)

Macrophage_gene <- intersect(Macrophage_gene, rownames(pericytes@data))


png(filename="DotPlot_pericytes_subclustered_Macrophage_genes.png", width=1000, height=1000, bg = "white", res = 150)
DotPlot(pericytes, Macrophage_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

pericytes_gene <- read.csv("pericytes_gene_set.txt", header=F)

pericytes_gene <-apply(pericytes_gene, 1, toupper)

pericytes_gene <- intersect(pericytes_gene, rownames(pericytes@data))


png(filename="DotPlot_pericytes_subclustered_pericytes_genes.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(pericytes, pericytes_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

Myeloid_gene <- read.csv("Myeloid_gene_set.txt", header=F)

Myeloid_gene <-apply(Myeloid_gene, 1, toupper)

Myeloid_gene <- intersect(Myeloid_gene, rownames(pericytes@data))


png(filename="DotPlot_pericytes_subclustered_Myeloid_genes.png", width=1500, height=1000, bg = "white", res = 150)
DotPlot(pericytes, Myeloid_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


neutroUp_gene <- read.csv("NeuroVsMono_UpNeutro.txt", header=F)

neutroUp_gene <-apply(neutroUp_gene, 1, toupper)

neutroUp_gene <- intersect(neutroUp_gene, rownames(pericytes@data))

png(filename="DotPlot_pericytes_subclustered_neutroUp_gene.png", width=1500, height=1000, bg = "white", res = 150)
DotPlot(pericytes, neutroUp_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


MonoUp_gene <- read.csv("NeuroVsMono_UpMono.txt", header=F)

MonoUp_gene <-apply(MonoUp_gene, 1, toupper)

MonoUp_gene <- intersect(MonoUp_gene, rownames(pericytes@data))

png(filename="DotPlot_pericytes_subclustered_MonoUp_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes, MonoUp_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


neutro_CellMarker_gene <- read.csv("CellMarker.Neutrophils.csv", header=T)

neutro_CellMarker_gene <- intersect(neutro_CellMarker_gene$Cell.Marker, rownames(pericytes.NeutroVsMono@data))

png(filename="DotPlot_pericytes_subclustered_neutro_CellMarker_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes, neutro_CellMarker_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

neutro_humanKidney_gene <- c("ADGRG3", "ANXA3", "AQP9", "ARG1", "BASP1", "BCL2A1", "BTNL8", "CDA", "CLC", "CLEC4E", "CMTM2", "CSF3R", "FCAR", "CXCL8", "CXCR2", 
"FCGR3B", "FPR1", "G0S2", "GCA", "IFITM2", "IL1R2", "LINC00694", "LST1", "MNDA", "NAMPT", "NCF1", "NFE4", "ORM1", "PGLYRP1", "PHOSPHO1", 
"PROK2", "PTGS2", "QPCT", "RGS2", "RP6-159A1.4", "S100A11", "S100A12", "S100A8", "S100A9", "S100P", "SOD2", "TNFAIP6", "USP10")

neutro_humanKidney_gene <- intersect(neutro_humanKidney_gene, rownames(pericytes@data))

png(filename="DotPlot_pericytes_subclustered_neutro_humanKidney_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes, neutro_humanKidney_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


###Renaming clusters

pericytes <- RenameIdent(pericytes, old.ident.name = 0, new.ident.name = "Acta2 cells")

pericytes <- RenameIdent(pericytes, old.ident.name = 1, new.ident.name = "Kcnj8 cells")

pericytes <- RenameIdent(pericytes, old.ident.name = 5, new.ident.name = "Mki67 cells")

pericytes <- RenameIdent(pericytes, old.ident.name = 6, new.ident.name = "Contamination")


png(filename="TSNEPlot-cluster_renamed.png", width=700, height=700, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Batch_renamed.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="ClusterTree.renamed.png", width=1000, height=1000, bg = "white", res = 150)
pericytes <- BuildClusterTree(object = pericytes, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  



pericytes <- StashIdent(pericytes, save.name = "Cell_type")

pericytes <- SubsetData(pericytes, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Contamination"), accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


batch_long <- pericytes@meta.data
batch_assigned <- batch_long 
batch_assigned <- batch_assigned %>%
   mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
batch_assigned <- as.data.frame(batch_assigned[,"Cond_TimePoint"])
 
rownames(batch_assigned) <- rownames(batch_long)
 
colnames(batch_assigned) <- "Cond_TimePoint"
 
pericytes <- AddMetaData(object = pericytes, metadata = batch_assigned, col.name = "Cond_TimePoint")

saveRDS(pericytes, "pericytes_subclustered_renamed.rds")

### Make distribution plot for P12 to P17

setwd(file.path(B_DIR, "Subclustering/Pericytes/CellIdentification"))

pericytes <- readRDS("pericytes_subclustered_renamed.rds")

png(filename="TSNEPlot-cluster_renamed_final.png", width=700, height=700, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="TSNEPlot-cluster_condition_renamed_final.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(pericytes, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Condition", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()




#Subset from P12 to P17

Subseted_cells_P5_P7_P10 <- rownames(subset(pericytes@meta.data, TimePoint %in% c("P5", "P7", "P10")))

pericytes.P5_P7_P10 <- SubsetData(pericytes, cells.use = Subseted_cells_P5_P7_P10, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


cluster_counts <- pericytes.P12_P14_P17@meta.data %>%
  group_by(Condition) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(Condition,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.P12_P14_P17.Condition.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.P12_P14_P17.Condition.png")


#Subset P17

Subseted_cells_P17 <- rownames(subset(pericytes@meta.data, TimePoint %in% c("P17")))

pericytes.P17 <- SubsetData(pericytes, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


cluster_counts <- pericytes.P17@meta.data %>%
  group_by(Condition) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(Condition,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.P17.Condition.png", width = 5, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.P17.Condition.png")



png("Pericytes_ACTA2_P17.png", width=2000, height=800, bg = "white", res = 150)
SplitDotPlotGG(pericytes.P17, grouping.var="Condition", genes.plot=c("ACTA2", "MYH11", "CSPG4", "PDGFRB", "CLDN5", "CDH5", "PECAM1"),
  cols.use = c("red", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 6, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()


png("Pericytes_ROBO_P17.png", width=1000, height=800, bg = "white", res = 150)
SplitDotPlotGG(pericytes.P17, grouping.var="Condition", genes.plot=c("ROBO1", "ROBO2"),
  cols.use = c("red", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 6, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png("pericytes.P5_P7_P10.proliferation.png", width=1000, height=800, bg = "white", res = 150)
SplitDotPlotGG(pericytes.P5_P7_P10, grouping.var="TimePoint", genes.plot=c("TBX18", "LAMA1", "COL1A1"),
  cols.use = c("red", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 6, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

, "MYO9B", "LAMA1", "COL1A1", "VHL"

###Subset OIR and NORM

Subseted_cells <- rownames(subset(pericytes@meta.data, Condition %in% c("OIR")))


pericytes.OIR <- SubsetData(pericytes, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

Subseted_cells <- rownames(subset(pericytes@meta.data, Condition %in% c("NORM")))


pericytes.NORM <- SubsetData(pericytes, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

cluster_counts <- pericytes.NORM@meta.data %>%
  group_by(TimePoint) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts <- cluster_counts %>% arrange(factor(cluster_counts$TimePoint, levels = c("P5", "P7", "P10", "P12", "P14", "P17")))

cluster_portions <-  ggplot(cluster_counts,aes(TimePoint,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()+
  scale_x_discrete(limits=cluster_counts$TimePoint)

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.NORM.TimePoint.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.NORM.TimePoint.png")

###Subset OIR and NORM for P12 P14 P17

Subseted_cells.OIR <- rownames(subset(pericytes.P12_P14_P17@meta.data, Condition %in% c("OIR")))

pericytes.OIR.P12_P14_P17 <- SubsetData(pericytes, cells.use = Subseted_cells.OIR, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

Subseted_cells.NORM <- rownames(subset(pericytes.P12_P14_P17@meta.data, Condition %in% c("NORM")))

pericytes.NORM.P12_P14_P17 <- SubsetData(pericytes, cells.use = Subseted_cells.NORM, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

cluster_counts <- pericytes.NORM.P12_P14_P17@meta.data %>%
  group_by(TimePoint) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(TimePoint,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.NORM.P12_P14_P17.TimePoint.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.NORM.P12_P14_P17.TimePoint.png")

###Find cluster with specific gene set for monocyte/ neutrophile

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/Mapping/Neutrophils"))

neutroUp_gene <- read.csv("NeuroVsMono_UpNeutro.txt", header=F)

neutroUp_gene <-apply(neutroUp_gene, 1, toupper)

neutroUp_gene <- intersect(neutroUp_gene, rownames(pericytes@data))

MonoUp_gene <- read.csv("NeuroVsMono_UpMono.txt", header=F)

MonoUp_gene <-apply(MonoUp_gene, 1, toupper)

MonoUp_gene <- intersect(MonoUp_gene, rownames(pericytes@data))



pericytes.NeutroVsMono <- RunPCA(object = pericytes, pc.genes = c(neutroUp_gene, MonoUp_gene) , do.print = TRUE, pcs.print = 1:10, 
    genes.print = 10)
    
PrintPCA(object = pericytes.NeutroVsMono, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

pericytes.NeutroVsMono <- ProjectPCA(object = pericytes.NeutroVsMono, do.print = FALSE)

png(filename="VizPCA.png")
VizPCA(object =pericytes.NeutroVsMono, pcs.use = 1:2)
dev.off()

png(filename="PCHeatmap_mult.png", width=2000, height=4000, bg = "white", res = 150)
PCHeatmap(object = pericytes.NeutroVsMono, pc.use = 1:12, cells.use = 100, do.balanced = TRUE, 
    label.columns = FALSE, use.full = FALSE)
dev.off()

png(filename="PCElbowPlot.png")
PCElbowPlot(object = pericytes.NeutroVsMono)
dev.off()

pericytes.NeutroVsMono <- RunTSNE(object = pericytes.NeutroVsMono, dims.use = 1:7, do.fast = T, dim_embed=2, perplexity=20, 
  reduction.type = "pca", reduction.name = "tsne_pca")


pericytes.NeutroVsMono <- FindClusters(object= pericytes.NeutroVsMono , reduction.type = "pca", dims.use = c(1:7), resolution = 0.5, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN = 1/15)


#FindClusters(object, genes.use = NULL, reduction.type = "pca",
#  dims.use = NULL, k.param = 30, plot.SNN = FALSE, prune.SNN = 1/15,
#  print.output = TRUE, distance.matrix = NULL, save.SNN = FALSE,
#  reuse.SNN = FALSE, force.recalc = FALSE, nn.eps = 0,
#  modularity.fxn = 1, resolution = 0.8, algorithm = 1, n.start = 100,
#  n.iter = 10, random.seed = 0, temp.file.location = NULL,
#  edge.file.name = NULL)
#Make plots

set.seed(001)

require(scales)

identities <- levels(pericytes.NeutroVsMono@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="TSNEPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes.NeutroVsMono, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(pericytes.NeutroVsMono, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="VlnPlotQC.pericytes.NeutroVsMono.ident.png", width=3000, height=1000, bg = "white", res = 150)
VlnPlot(object = pericytes.NeutroVsMono, c("nGene", "nUMI", "percent.mito"), group.by = "ident", nCol = 3, x.lab.rot = TRUE)
dev.off()


#####Find markers

oldidentname <- c("1")
for (clusterID in oldidentname) {
  pericytes.NeutroVsMono <- RenameIdent(pericytes.NeutroVsMono, old.ident.name = clusterID, new.ident.name = "0")
} 

png(filename="TSNEPlot-cluster.pca.regrouped.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes.NeutroVsMono, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="VlnPlotQC.pericytes.NeutroVsMono.ident.regrouped.png", width=3000, height=1000, bg = "white", res = 150)
VlnPlot(object = pericytes.NeutroVsMono, c("nGene", "nUMI", "percent.mito"), group.by = "ident", nCol = 3, x.lab.rot = TRUE)
dev.off()


pericytes.NeutroVsMono.markers <- FindAllMarkers(object = pericytes.NeutroVsMono, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

pericytes.NeutroVsMono.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file ="pericytes.NeutroVsMono-subclustered-Markers-top10.tsv", sep = "\t")

top10.markers <- unique(top10$gene)

png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_top10_markers.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="ClusterTree.png", width=1000, height=1000, bg = "white", res = 150)
pericytes.NeutroVsMono <- BuildClusterTree(object = pericytes.NeutroVsMono, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

pericytes.NeutroVsMono.markersNode <- FindAllMarkersNode(object =pericytes.NeutroVsMono, node = NULL, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)

pericytes.NeutroVsMono.markersNode %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10.node

write.table(top10.node, file ="pericytes.NeutroVsMono-subclustered-MarkersNode-top10.avg_logFC.tsv", sep = "\t")

###Gene set for cluster identification

neutro_gene <- read.csv("Neutrophil_gene_set.txt", header=F)

neutro_gene <-apply(neutro_gene, 1, toupper)

neutro_gene <- intersect(neutro_gene, rownames(pericytes.NeutroVsMono@data))


png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_neutro_genes.png", width=1000, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, neutro_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="FeaturePlot_pericytes.NeutroVsMono_subclustered_TLR2_CD14.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(pericytes.NeutroVsMono, c("TLR2", "CD14"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 2,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_pericytes.NeutroVsMono_subclustered_ICAM1_CD63.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(pericytes.NeutroVsMono, c("ICAM1", "CD63"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 2,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()


Macrophage_gene <- read.csv("Macrophage_gene_set.txt", header=F)

Macrophage_gene <-apply(Macrophage_gene, 1, toupper)

Macrophage_gene <- intersect(Macrophage_gene, rownames(pericytes.NeutroVsMono@data))


png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_Macrophage_genes.png", width=1000, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, Macrophage_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

pericytes_gene <- read.csv("pericytes_gene_set.txt", header=F)

pericytes_gene <-apply(pericytes_gene, 1, toupper)

pericytes_gene <- intersect(pericytes_gene, rownames(pericytes.NeutroVsMono@data))


png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_pericytes_genes.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, pericytes_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

Myeloid_gene <- read.csv("Myeloid_gene_set.txt", header=F)

Myeloid_gene <-apply(Myeloid_gene, 1, toupper)

Myeloid_gene <- intersect(Myeloid_gene, rownames(pericytes.NeutroVsMono@data))


png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_Myeloid_genes.png", width=1500, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, Myeloid_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


neutroUp_gene <- read.csv("NeuroVsMono_UpNeutro.txt", header=F)

neutroUp_gene <-apply(neutroUp_gene, 1, toupper)

neutroUp_gene <- intersect(neutroUp_gene, rownames(pericytes.NeutroVsMono@data))

png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_neutroUp_gene.png", width=1500, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, neutroUp_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


MonoUp_gene <- read.csv("NeuroVsMono_UpMono.txt", header=F)

MonoUp_gene <-apply(MonoUp_gene, 1, toupper)

MonoUp_gene <- intersect(MonoUp_gene, rownames(pericytes.NeutroVsMono@data))

png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_MonoUp_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, MonoUp_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

neutro_binet_gene <- read.csv("Neutrophil_Binet_gene_set.txt", header=F)

neutro_binet_gene <-apply(neutro_binet_gene, 1, toupper)

neutro_binet_gene <- intersect(neutro_binet_gene, rownames(pericytes.NeutroVsMono@data))

png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_neutro_Binet_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, neutro_binet_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


neutro_CellMarker_gene <- read.csv("CellMarker.Neutrophils.csv", header=T)

neutro_CellMarker_gene <- intersect(neutro_CellMarker_gene$Cell.Marker, rownames(pericytes.NeutroVsMono@data))

png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_neutro_CellMarker_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, neutro_CellMarker_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


neutro_humanKidney_gene <- c("ADGRG3", "ANXA3", "AQP9", "ARG1", "BASP1", "BCL2A1", "BTNL8", "CDA", "CLC", "CLEC4E", "CMTM2", "CSF3R", "FCAR", "CXCL8", "CXCR2", 
"FCGR3B", "FPR1", "G0S2", "GCA", "IFITM2", "IL1R2", "LINC00694", "LST1", "MNDA", "NAMPT", "NCF1", "NFE4", "ORM1", "PGLYRP1", "PHOSPHO1", 
"PROK2", "PTGS2", "QPCT", "RGS2", "RP6-159A1.4", "S100A11", "S100A12", "S100A8", "S100A9", "S100P", "SOD2", "TNFAIP6", "USP10")

neutro_humanKidney_gene <- intersect(neutro_humanKidney_gene, rownames(pericytes.NeutroVsMono@data))

png(filename="DotPlot_pericytes.NeutroVsMono_subclustered_neutro_humanKidney_gene.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(pericytes.NeutroVsMono, neutro_humanKidney_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


pericytes.NeutroVsMono <- StashIdent(pericytes.NeutroVsMono, save.name = "Cell_type")


saveRDS(pericytes.NeutroVsMono, "pericytes.NeutroVsMono_subclustered_renamed.rds")

### Make distribution plot for P12 to P17
setwd(file.path(B_DIR, "Subclustering/ImmuneCells/Mapping/Neutrophils"))

pericytes.NeutroVsMono <- readRDS("pericytes.NeutroVsMono_subclustered_renamed.rds")

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/CellIdentification/Neutrophils"))


#Subset from P12 to P17

Subseted_cells_P12_P14_P17 <- rownames(subset(pericytes.NeutroVsMono@meta.data, TimePoint %in% c("P12", "P14", "P17")))

pericytes.NeutroVsMono.P12_P14_P17 <- SubsetData(pericytes.NeutroVsMono, cells.use = Subseted_cells_P12_P14_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


cluster_counts <- pericytes.NeutroVsMono.P12_P14_P17@meta.data %>%
  group_by(Condition) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(Condition,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.P12_P14_P17.Condition.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.P12_P14_P17.Condition.png")


#Subset P17

Subseted_cells_P17 <- rownames(subset(pericytes.NeutroVsMono@meta.data, TimePoint %in% c("P17")))

pericytes.NeutroVsMono.P17 <- SubsetData(pericytes.NeutroVsMono, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


cluster_counts <- pericytes.NeutroVsMono.P17@meta.data %>%
  group_by(Condition) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(Condition,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.P17.Condition.png", width = 5, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.P17.Condition.png")


###Subset OIR and NORM

Subseted_cells <- rownames(subset(pericytes.NeutroVsMono@meta.data, Condition %in% c("OIR")))


pericytes.NeutroVsMono.OIR <- SubsetData(pericytes.NeutroVsMono, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

Subseted_cells <- rownames(subset(pericytes.NeutroVsMono@meta.data, Condition %in% c("NORM")))


pericytes.NeutroVsMono.NORM <- SubsetData(pericytes.NeutroVsMono, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

cluster_counts <- pericytes.NeutroVsMono.NORM@meta.data %>%
  group_by(TimePoint) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts <- cluster_counts %>% arrange(factor(cluster_counts$TimePoint, levels = c("P5", "P7", "P10", "P12", "P14", "P17")))

cluster_portions <-  ggplot(cluster_counts,aes(TimePoint,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()+
  scale_x_discrete(limits=cluster_counts$TimePoint)

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.NORM.TimePoint.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.NORM.TimePoint.png")

###Subset OIR and NORM for P12 P14 P17

Subseted_cells.OIR <- rownames(subset(pericytes.NeutroVsMono.P12_P14_P17@meta.data, Condition %in% c("OIR")))

pericytes.NeutroVsMono.OIR.P12_P14_P17 <- SubsetData(pericytes.NeutroVsMono, cells.use = Subseted_cells.OIR, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

Subseted_cells.NORM <- rownames(subset(pericytes.NeutroVsMono.P12_P14_P17@meta.data, Condition %in% c("NORM")))

pericytes.NeutroVsMono.NORM.P12_P14_P17 <- SubsetData(pericytes.NeutroVsMono, cells.use = Subseted_cells.NORM, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

cluster_counts <- pericytes.NeutroVsMono.OIR.P12_P14_P17@meta.data %>%
  group_by(TimePoint) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(TimePoint,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.ImmuneCells.OIR.P12_P14_P17.TimePoint.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.ImmuneCells.OIR.P12_P14_P17.TimePoint.png")


##OVerley clusters from different clusterisation type

cells_cluster_2 <- FastWhichCells(object = pericytes.NeutroVsMono, group.by = "ident", subset.value = 2)

#cells_cluster_0 <- FastWhichCells(object = pericytes.NeutroVsMono, group.by = "ident", subset.value = 0)

pericytes_test <- SetIdent(pericytes, cells.use = cells_cluster_2, ident.use = "Cluster 2")

#pericytes_test <- SetIdent(pericytes_test, cells.use = cells_cluster_0, ident.use = "Cluster 0")

png(filename="TSNEPlot-cluster.pericytes_cluster2_test.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(pericytes_test, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


#END
q("no")
