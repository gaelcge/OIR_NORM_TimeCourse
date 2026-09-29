# ---------------------------------------------------------------------------
# 05_glial_ecs - 01_subcluster.R
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# glia and endothelium together, run on cells taken from the annotated
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


setwd(file.path(B_DIR, "CellTypeAnnotation/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

NeuroGlia <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))

Seurat_object <- NeuroGlia 
summary(as.factor(Seurat_object@meta.data$Cond_Sorting))

summary(Seurat_object@meta.data)

####Subclustering

setwd(file.path(B_DIR, "Subclustering/Glial_ECs/Mapping"))

Subseted_cells_P12_P14_P17 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P12", "P14", "P17")))

Seurat_object <- SubsetData(Seurat_object, cells.use = Subseted_cells_P12_P14_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Quiescent Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Activated Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Mki67 Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Early Muller Glia"), new.ident.name = "Muller glia")


markers.retina.dotplot <- rev(c("RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "RLBP1", "FBN1", "CLDN5", "LYZ2", "GFAP", "OPTC", "TOP2A", "MKI67"))

Seurat_object <- ReorderIdent(Seurat_object, feature = c("VEGFA", "GFAP", "LYZ2", "PECAM1"), rev = FALSE, aggregate.fxn = mean,
  reorder.numeric = FALSE)

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.VEGFFamilly.png", sep="."), res = 150, width=800, height=800)
DotPlot(Seurat_object, genes.plot= rev(c("VEGFA", "RLBP1", "GFAP", "LYZ2", "PECAM1")), cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name =NULL, ident.use = c("Astrocytes", "Muller glia", "Immune cells", "Endothelial cells"),
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)


set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="UMAPPlot-Subclustered_Seurat_object.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Subclustered_Seurat_object.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

Seurat_object <- ScaleData(object = Seurat_object, vars.to.regress = c("nUMI", "Batch", "percent.mito", "percent.crystal"))

png(filename="MeanVarPlot.final.png")
Seurat_object <- FindVariableGenes(object = Seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.15, x.high.cutoff = 4, y.cutoff = 0.5)
dev.off()

length(x=Seurat_object@var.genes)


##DImension reduction

## PCA and ICA

Seurat_object <- RunPCA(object = Seurat_object, pc.genes = Seurat_object@var.genes, do.print = TRUE, pcs.print = 1:10, 
    genes.print = 10)
    
PrintPCA(object = Seurat_object, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

Seurat_object <- ProjectPCA(object = Seurat_object, do.print = FALSE)

Seurat_object <- RunICA(object = Seurat_object, ic.genes = Seurat_object@var.genes)

Seurat_object <- ProjectDim(Seurat_object, reduction.type = "ica", dims.print = 1:5,
  dims.store = 30, genes.print = 30, replace.dim = FALSE,
  do.center = FALSE, do.print = TRUE, assay.type = "RNA")


png(filename="VizPCA.png", width=800, height=2000, bg = "white", res = 150)
VizPCA(object =Seurat_object, pcs.use = 1:15)
dev.off()

png(filename="VizICA.png", width=800, height=2000, bg = "white", res = 150)
VizICA(object =Seurat_object, ics.use = 1:15)
dev.off()

png(filename="PCAPlot.batch.pc1vs2.pca.png")
DimPlot(object =Seurat_object, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()

png(filename="ICAPlot.batch.pc1vs2.pca.png")
DimPlot(object =Seurat_object, reduction.use = "ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()

png(filename="PCAPlot.batch.pc1vs2.pca.Cond_Sorting.png")
DimPlot(object =Seurat_object, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()

png(filename="ICAPlot.batch.pc1vs2.pca.Cond_Sorting.png")
DimPlot(object =Seurat_object, reduction.use = "ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()


png(filename="PCHeatmap_mult.pca.png", width=2000, height=4000, bg = "white", res = 150)
DimHeatmap(Seurat_object, assay.use = "RNA", reduction.type = "pca", dim.use = 1:12,
  cells.use = 100, num.genes = 30, use.full = FALSE, disp.min = -2.5,
  disp.max = 2.5, do.return = FALSE, col.use = PurpleAndYellow(),
  use.scale = TRUE, do.balanced = TRUE, remove.key = FALSE,
  label.columns = NULL, check.plot = TRUE)
dev.off()

png(filename="ICHeatmap_mult.ica.png", width=2000, height=4000, bg = "white", res = 150)
DimHeatmap(Seurat_object, assay.use = "RNA", reduction.type = "ica", dim.use = 1:12,
  cells.use = 100, num.genes = 30, use.full = FALSE, disp.min = -2.5,
  disp.max = 2.5, do.return = FALSE, col.use = PurpleAndYellow(),
  use.scale = TRUE, do.balanced = TRUE, remove.key = FALSE,
  label.columns = NULL, check.plot = TRUE)
dev.off()

png(filename="PCElbowPlot.png")
PCElbowPlot(object = Seurat_object)
dev.off()

png(filename="ICAElbowPlot.png")
DimElbowPlot(Seurat_object, reduction.type = "ica", dims.plot = 20)
dev.off()


PrintDim(object = Seurat_object, reduction.type = "pca", dims.print = 1:2, genes.print = 10)

p1 <- DimPlot(object = Seurat_object, reduction.use = "pca", group.by = "Batch", pt.size = 1, 
    do.return = TRUE)
p2 <- VlnPlot(object = Seurat_object, features.plot = "PC1", group.by = "Batch", do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_pca.Batch.png", width=1500, height=1000, bg = "white", res = 100)
plot_grid(p1, p2)
dev.off()

PrintDim(object = Seurat_object, reduction.type = "ica", dims.print = 1:2, genes.print = 10)

p1 <- DimPlot(object = Seurat_object, reduction.use = "ica", group.by = "Batch", pt.size = 1, 
    do.return = TRUE)
p2 <- VlnPlot(object = Seurat_object, features.plot = "IC1", group.by = "Batch", do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_ica.Batch.png", width=1500, height=1000, bg = "white", res = 100)
plot_grid(p1, p2)
dev.off()


###### RUn harmony

Seurat_object <- RunHarmony(Seurat_object, "Cond_Sorting", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100)

Seurat_object <- ProjectDim(Seurat_object, reduction.type = "harmony", dims.print = 1:5,
  dims.store = 30, genes.print = 30, replace.dim = FALSE,
  do.center = FALSE, do.print = TRUE, assay.type = "RNA")


png(filename="VizHarmony.png", width=800, height=2000, bg = "white", res = 150)
VizDimReduction(Seurat_object, reduction.type = "harmony", dims.use = 1:15,
  num.genes = 30, use.full = FALSE, font.size = 0.5, nCol = NULL,
  do.balanced = FALSE)
dev.off()


png(filename="PCAPlot.batch.pc1vs2.harmony.png")
DimPlot(object =Seurat_object, reduction.use = "harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 3, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.label = FALSE, label.size = 1, no.legend = FALSE)
dev.off()


png(filename="PCHeatmap_mult.final.png", width=2000, height=4000, bg = "white", res = 150)
PCHeatmap(object = Seurat_object, pc.use = 1:15, cells.use = 100, do.balanced = TRUE, 
    label.columns = FALSE, use.full = FALSE)
dev.off()


png(filename="HarmonyHeatmap_mult.png", width=2000, height=4000, bg = "white", res = 150)
DimHeatmap(Seurat_object, assay.use = "RNA", reduction.type = "harmony", dim.use = 1:12,
  cells.use = 100, num.genes = 30, use.full = FALSE, disp.min = -2.5,
  disp.max = 2.5, do.return = FALSE, col.use = PurpleAndYellow(),
  use.scale = TRUE, do.balanced = TRUE, remove.key = FALSE,
  label.columns = NULL, check.plot = TRUE)
dev.off()

# No standard deviation info stored for harmony, so elbow on pca
png(filename="PCAElbowPlot.pca_harmony.png")
DimElbowPlot(Seurat_object, reduction.type = "pca", dims.plot = 20)
dev.off()


PrintDim(object = Seurat_object, reduction.type = "harmony", dims.print = 1:2, genes.print = 10)

p1 <- DimPlot(object = Seurat_object, reduction.use = "harmony", group.by = "Batch", pt.size = 1, 
    do.return = TRUE)
p2 <- VlnPlot(object = Seurat_object, features.plot = "Harmony1", group.by = "Batch", do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_Harmony.Batch.png", width=1500, height=1000, bg = "white", res = 100)
plot_grid(p1, p2)
dev.off()


##Rn TSNE

Seurat_object <- RunTSNE(object = Seurat_object, dims.use = 1:12, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.type = "pca", reduction.name = "tsne_pca")

Seurat_object <- RunTSNE(object = Seurat_object, dims.use = 1:12, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.type = "ica", reduction.name = "tsne_ica")

Seurat_object <- RunTSNE(object = Seurat_object, dims.use = 1:12, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.type = "harmony", reduction.name = "tsne_harmony")


##Run UMAP

Seurat_object <- RunUMAP(Seurat_object, cells.use = NULL, dims.use = 1:12, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_pca", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)

Seurat_object <- RunUMAP(Seurat_object, cells.use = NULL, dims.use = 1:4, reduction.use = "ica",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_ica", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)

Seurat_object <- RunUMAP(Seurat_object, cells.use = NULL, dims.use = 1:4, reduction.use = "harmony",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_harmony", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)


#######Find cluster

##For PCA
Seurat_object <- FindClusters(object= Seurat_object , reduction.type = "pca", dims.use = c(1:12), resolution = 2, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 30,  algorithm = 1, prune.SNN=1/15)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="TSNEPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-TimePoint.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="DotPlot_markers.ident.VEGFFamilly.pca.png", res = 150, width=800, height=800)
DotPlot(Seurat_object, genes.plot= rev(c("VEGFA", "RLBP1", "GFAP", "LYZ2", "PECAM1")), cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()



png(filename="UMAPPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


##For ICA

Seurat_object <- FindClusters(object= Seurat_object , reduction.type = "ica", dims.use = c(1:12), resolution = 2, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN=1/8)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))


png(filename="TSNEPlot-cluster.ica.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.ica.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-TimePoint.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="UMAPPlot-cluster.ica.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.ica.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_ica", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



##For harmony

Seurat_object <- FindClusters(object= Seurat_object , reduction.type = "harmony", dims.use = c(1:12), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN=1/8)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))


png(filename="TSNEPlot-cluster.harmony.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.harmony.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-TimePoint.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="UMAPPlot-cluster.harmony.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.harmony.png", width=800, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_harmony", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

### Choose dimensionality reduction for cluster finding

Seurat_object <- FindClusters(object= Seurat_object , reduction.type = "pca", dims.use = c(1:12), resolution = 2, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 30,  algorithm = 1, prune.SNN=1/15)


##### QC on clusters

png(filename="VlnPlotQC.Seurat_object.pca.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(Seurat_object, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

png(filename="Tree_initial.png", res = 150, width=1000, height=1500)
BuildClusterTree(Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()


##Find cell types

setwd(file.path(B_DIR, "Subclustering/Glial_ECs/CellIdentification"))


## Regrouping

oldidentname <- c("2", "3", "5", "6", "8", "9", "10", "12")
for (clusterID in oldidentname) {
  Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "0")
} 

oldidentname <- c("15", "16")
for (clusterID in oldidentname) {
  Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "11")
} 

oldidentname <- c("14")
for (clusterID in oldidentname) {
  Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "1")
} 


Seurat_object <-ValidateSpecificClusters(Seurat_object, cluster1 = 1, cluster2 = 8,
  pc.use = 1:5, top.genes = 30, acc.cutoff = 0.9)

png(filename="TSNEPlot-cluster.pca.regrouping.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster.pca.regrouping.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="UmapPlot-pANNPredictions.pca.regrouping.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "pANNPredictions", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="Tree_regrouping.png", res = 150, width=1000, height=1500)
BuildClusterTree(Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()

png(filename="VlnPlotQC.Seurat_object.harmony.regrouping.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(Seurat_object, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

#####

Seurat_object.markers <- FindAllMarkers(object = Seurat_object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

Seurat_object.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file ="Seurat_object-subclustered-Markers-top10.tsv", sep = "\t")

top10.markers <- unique(top10$gene)

png(filename="DotPlot_Seurat_object_subclustered_top10_markers.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(Seurat_object, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="ClusterTree.png", width=1000, height=1000, bg = "white", res = 150)
Seurat_object <- BuildClusterTree(object = Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

Seurat_object.markersNode <- FindAllMarkersNode(object =Seurat_object, node = NULL, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)

Seurat_object.markersNode %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10.node

write.table(top10.node, file ="Seurat_object-subclustered-MarkersNode-top10.avg_logFC.tsv", sep = "\t")



###Gene set for cluster identification

neutro_gene <- read.csv("Neutrophil_gene_set.txt", header=F)

neutro_gene <-apply(neutro_gene, 1, toupper)

neutro_gene <- intersect(neutro_gene, rownames(Seurat_object@data))


png(filename="DotPlot_Seurat_object_subclustered_neutro_genes.png", width=1000, height=1000, bg = "white", res = 150)
DotPlot(Seurat_object, neutro_gene, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="FeaturePlot_Seurat_object_subclustered_NDUFA4L2_VEGFA.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("NDUFA4L2", "VEGFA"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 2,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_Seurat_object_subclustered_ICAM1_CD63.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("ICAM1", "CD63"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 2,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

markers.retina.dotplot <- rev(c("RHO", "LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "OTX2", "PRDM1", "RLBP1", "GFAP", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL", "TOP2A", "MKI67"))

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.regrouped.png", sep="."), res = 150, width=1500, height=1000)
DotPlot(Seurat_object, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


###Renaming clusters

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 0, new.ident.name = "Muller Glia")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 1, new.ident.name = "Rho Contamination")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 4, new.ident.name = "Activated Muller Glia")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 7, new.ident.name = "Endothelial cells")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 11, new.ident.name = "Microglia")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 13, new.ident.name = "Astrocytes")

Seurat_object <- StashIdent(Seurat_object, save.name = "Cell_type")


png(filename="TSNEPlot-cluster_renamed.png", width=700, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Batch_renamed.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Batch", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="DotPlot_Seurat_object_subclustered_Markers.renamed.png", width=3000, height=800, bg = "white", res = 150)
DotPlot(Seurat_object, unique(top10$gene), cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 5,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


png(filename="ClusterTree.renamed.png", width=1000, height=1000, bg = "white", res = 150)
Seurat_object <- BuildClusterTree(object = Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

#### Remove unwanted clusters to re-cluster and make final plot

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Rho Contamination"), accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

png(filename="UmapPlot-cluster.pca.regrouping.final.png", width=500, height=500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


Seurat_object <- ScaleData(object = Seurat_object, vars.to.regress = c("nUMI", "percent.mito", "Batch", "percent.crystal"),
  do.par=TRUE, num.cores=10)

png(filename="MeanVarPlot.final.png")
Seurat_object <- FindVariableGenes(object = Seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.15, x.high.cutoff = 4, y.cutoff = 0.5)
dev.off()

length(x=Seurat_object@var.genes)

Seurat_object <- RunPCA(object = Seurat_object, pc.genes = Seurat_object@var.genes, do.print = TRUE, pcs.print = 1:10, 
    genes.print = 10, weight.by.var = FALSE)

png(filename="PCHeatmap_mult.final.png", width=2000, height=4000, bg = "white", res = 150)
PCHeatmap(object = Seurat_object, pc.use = 1:15, cells.use = 100, do.balanced = TRUE, 
    label.columns = FALSE, use.full = FALSE)
dev.off()

# No standard deviation info stored for harmony, so elbow on pca
#png(filename="PCAElbowPlot.pca_final.png")
#DimElbowPlot(Seurat_object, reduction.type = "pca", dims.plot = 20)
#dev.off()

Seurat_object <- RunTSNE(object = Seurat_object, dims.use = 1:10, do.fast = T, dim_embed=2, perplexity=50, 
  reduction.use = "pca", reduction.name = "tsne_pca")

Seurat_object <- RunUMAP(Seurat_object, cells.use = NULL, dims.use = 1:12, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_pca", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)


Seurat_object <- FindClusters(object= Seurat_object , reduction.type = "pca", dims.use = c(1:10), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 50,  algorithm = 1, prune.SNN=1/15)


png(filename="TsnePlot-cluster.ident.final.png", width=900, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TsnePlot-cluster.Cell_type.final.png", width=900, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cell_type", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="UmapPlot-cluster.ident.final.png", width=900, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()




oldidentname <- c("0", "1", "2", "3", "4")
for (clusterID in oldidentname) {
  Seurat_object <- RenameIdent(Seurat_object, old.ident.name = clusterID, new.ident.name = "Muller Glia")
} 

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 5, new.ident.name = "Activated Muller Glia")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 7, new.ident.name = "Endothelial cells")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 8, new.ident.name = "Astrocytes")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 9, new.ident.name = "Microglia")

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = 6, new.ident.name = "Microglia")


Seurat_object <- StashIdent(Seurat_object, save.name = "Cell_type_final")


set.seed(001)

require(scales)

identities <- levels(Seurat_object@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(55, 360), c = 150, l = 60, h.start = 0, direction = 1)(length(identities)))

png(filename="FeaturePlot_Seurat_object_subclustered_nGene.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("nGene"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("blue","red"), pch.use = 16, overlay = FALSE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_Seurat_object_subclustered_nUMI.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("nUMI"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("blue","red"), pch.use = 16, overlay = FALSE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_Seurat_object_subclustered_VEGFA.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("VEGFA"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("blue","red"), pch.use = 16, overlay = FALSE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_Seurat_object_subclustered_NDUFA4L2.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("NDUFA4L2"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("blue","red"), pch.use = 16, overlay = FALSE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="FeaturePlot_Seurat_object_subclustered_NDUFA4L2_GFAP.png", width=900, height=500, bg = "white", res = 150)
FeaturePlot(Seurat_object, c("NDUFA4L2", "GFAP"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("lightblue", "orange", "yellow", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_pca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()


png(filename="VlnPlotQC.Seurat_object.Subtype.renamed.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(Seurat_object, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = NULL, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

png(filename="TsnePlot-Condition.renamed.png", width=1100, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = c("mediumseagreen", "mediumslateblue"), group.by = "Condition", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="Tree_final.png", res = 150, width=1000, height=1500)
BuildClusterTree(Seurat_object, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()


png(filename="TsnePlot-cluster.renamed.final.png", width=900, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-cluster_TimePoint_renamed_final.png", width=1100, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster.renamed.final.png", width=900, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 6, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster_TimePoint_renamed_final.png", width=1100, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


saveRDS(Seurat_object, "Seurat_object_subclustered_renamed.rds")

### Make distribution plot for P12 to P17

setwd(file.path(B_DIR, "Subclustering/Glial_ECs/CellIdentification"))

Seurat_object <- readRDS("Seurat_object_subclustered_renamed.rds")

batch_long <- Seurat_object@meta.data
batch_assigned <- batch_long 
batch_assigned <- batch_assigned %>%
   mutate("Cond_cellType"=paste(Condition,Cell_type_final,sep="_"))
 
batch_assigned <- as.data.frame(batch_assigned[,"Cond_cellType"])
 
rownames(batch_assigned) <- rownames(batch_long)
 
colnames(batch_assigned) <- "Cond_cellType"
 
Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned, col.name = "Cond_cellType")

saveRDS(Seurat_object, "Seurat_object_subclustered_renamed.rds")



png(filename="TSNEPlot-cluster_renamed_final.png", width=700, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="TSNEPlot-cluster_condition_renamed_final.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Condition", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-cluster_TimePoint_renamed_final.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-cluster_Cond_TimePoint_renamed_final.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(Seurat_object, reduction.use = "tsne_pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


#Subset OIR and NORM

Subseted_cells_NORM <- rownames(subset(Seurat_object@meta.data, Condition %in% c("NORM")))

Seurat_object.NORM <- SubsetData(Seurat_object, cells.use = Subseted_cells_NORM, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.OIR, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png(filename="UmapPlot-cluster.renamed.final.OIR.png", width=900, height=700, bg = "white", res = 150)
#~ DimPlot(Seurat_object.OIR, reduction.use = "umap_pca", dim.1 = 1, dim.2 = 2,
#~   cells.use = NULL, pt.size = 2, do.return = FALSE, do.bare = FALSE,
#~   cols.use = colors, group.by = "ident", pt.shape = NULL,
#~   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#~   do.label = FALSE, label.size = 6, no.legend = TRUE, no.axes = FALSE,
#~   dark.theme = FALSE)
#~ dev.off()



Seurat_object_Cond_cellType <- SetAllIdent(Seurat_object, id="Cond_cellType")

Seurat_object.markers <- FindMarkers(Seurat_object_Cond_cellType, ident.1="OIR_Activated Muller Glia", ident.2 ="NORM_Muller Glia")

top100 <- head(Seurat_object.markers, 100)

write.table(top100, file ="Seurat_object-ActivatedMullerglia.DEG.wilcox.tsv", sep = "\t")

top100 <- read.table("Seurat_object-ActivatedMullerglia.DEG.wilcox.tsv")

top10.markers <- unique(rownames(head(top100,10)))

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png(filename="DotPlot_Seurat_object_ActivatedMullerglia.DEG.wilcox.top10.png", width=1200, height=800, bg = "white", res = 200)
#~ DotPlot(Seurat_object, rev(top10.markers), cols.use = c("blue", "red"),
#~   col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
#~   scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
#~   plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()
#~
#~
#~
#~ cluster_counts <- Seurat_object.P12_P14_P17@meta.data %>%
#~   group_by(Cond_TimePoint) %>%
#~   count(Cell_type_final) %>%
#~   mutate(freq = n / sum(n))
#~
#~ cluster_counts$Cell_type_final <- factor(cluster_counts$Cell_type_final,levels=unique(mixedsort(cluster_counts$Cell_type_final)))
#~
#~ cluster_portions <-  ggplot(cluster_counts,aes(Cond_TimePoint,freq,fill=Cond_TimePoint)) +
#~   geom_bar(stat="identity",col="black") +
#~   facet_wrap(~ Cell_type_final, scale="free") +
#~   theme_light()
#~
#~ ggsave(cluster_portions,file="Cluster_proportions_CellType.Seurat_object.P12_P14_P17.Cond_TimePoint.png", width = 12, height = 5)
#~
#~ cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type_final,freq,fill=Cond_TimePoint)) +
#~   geom_bar(stat="identity",col="black") +
#~   theme_light()
#~
#~ ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.Seurat_object.P12_P14_P17.Cond_TimePoint.png")
#~
#~
#~ #Subset P17
#~
#~ Subseted_cells_P14 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14")))
#~
#~ Seurat_object.P14 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P14, subset.name =NULL, ident.use = NULL,
#~   ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
#~   do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
#~   random.seed = 1)
#~
#~
#~ cluster_counts <- Seurat_object.P14@meta.data %>%
#~   group_by(Condition) %>%
#~   count(Cell_type_final) %>%
#~   mutate(freq = n / sum(n))
#~
#~ cluster_counts$Cell_type_final <- factor(cluster_counts$Cell_type_final,levels=unique(mixedsort(cluster_counts$Cell_type_final)))
#~
#~ cluster_portions <-  ggplot(cluster_counts,aes(Condition,freq,fill=Condition)) +
#~   geom_bar(stat="identity",col="black") +
#~   facet_wrap(~ Cell_type_final, scale="free") +
#~   theme_light()
#~
#~ ggsave(cluster_portions,file="Cluster_proportions_CellType.Seurat_object.P14.Condition.png", width = 6, height = 3)
#~
#~ cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type_final,freq,fill=Condition)) +
#~   geom_bar(stat="identity",col="black") +
#~   theme_light()
#~
#~ ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.Seurat_object.P14.Condition.png")
#~
#~
#~ png("Seurat_object_VEGFA_P14.png", width=900, height=600, bg = "white", res = 150)
#~ SplitDotPlotGG(Seurat_object.P14, grouping.var="Condition", genes.plot=c("VEGFA"),
#~   cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 6, group.by="Ident", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()


#END
q("no")
