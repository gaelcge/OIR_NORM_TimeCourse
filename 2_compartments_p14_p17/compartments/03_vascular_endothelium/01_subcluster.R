# ---------------------------------------------------------------------------
# 03_vascular_endothelium - 01_subcluster.R
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# vascular endothelium, run on cells taken from the annotated
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
library(SingleR)

setwd(file.path(B_DIR, "Clustering_noDiscarding/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

retina <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.rds", sep="."))

summary(as.factor(retina@meta.data$Cond_Sorting))


####Subclustering

setwd(file.path(B_DIR, "Subclustering/VascularEndothelium/Mapping"))


markers.retina.dotplot <- rev(c("RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "RLBP1", "FBN1", "CLDN5", "LYZ2", "GFAP", "OPTC", "TOP2A", "MKI67"))

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.initial.png", sep="."), res = 150, width=1000, height=2000)
DotPlot(retina, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


vascendothelium <- SubsetData(retina, cells.use = NULL, subset.name =NULL, ident.use = c(29),
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

summary(as.factor(vascendothelium@meta.data$Cond_Sorting))

vascendothelium <- StashIdent(vascendothelium, save.name = "IniClustNBIdent")

set.seed(001)

require(scales)

identities <- levels(vascendothelium@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="UMAPPlot-Subclustered_vascendothelium.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TsnePlot-Subclustered_vascendothelium.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


vascendothelium <- ScaleData(object = vascendothelium, vars.to.regress = c("nUMI", "percent.mito", "Batch", "percent.crystal"))

png("Variable_gene_vascendothelium.png")
vascendothelium <- FindVariableGenes(object = vascendothelium, do.plot = T, display.progress = F, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.20, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()


vascendothelium <- RunPCA(object = vascendothelium, features = vascendothelium@var.genes)

png("ElbowPlot_pca.png")
PCElbowPlot(object = vascendothelium)
dev.off()

png(filename="DimHeatmap.pca.png", width=1200, height=1700, res=150)
DimHeatmap(object = vascendothelium, reduction.type = "pca", cells.use = 100, 
    dim.use = 1:12, do.balanced = TRUE)
dev.off()


vascendothelium <- RunTSNE(object = vascendothelium, dims.use = 1:6, do.fast = T, dim_embed=2, perplexity=20, reduction.use = "pca",
  reduction.name = "tsne", tsne.method = "Rtsne")

vascendothelium <- RunUMAP(vascendothelium, cells.use = NULL, dims.use = 1:6, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap", reduction.key = "UMAP", n_neighbors = 30L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)



vascendothelium <- FindClusters(object= vascendothelium , reduction.type = "pca", dims.use = c(1:10), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1, prune.SNN=1/15)

png(filename="TSNEPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-TimePoint.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "TimePoint", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


png(filename="UMAPPlot-cluster.pca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.pca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



#png(filename="MeanVarPlot.png")
#vascendothelium <- FindVariableGenes(object = vascendothelium, mean.function = ExpMean, dispersion.function = LogVMR, 
#    x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.5)
#dev.off()

#length(x=vascendothelium@var.genes)


###Subset the different dataset

test <- unique(vascendothelium@meta.data$Cond_Sorting)

for (i in test) {
  subset_cells <- rownames(subset(vascendothelium@meta.data, Sorting %in% i))
  subset_retina <- SubsetData(vascendothelium, cells.use = subset_cells, subset.name = NULL, ident.use = NULL,
    ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
    do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
    random.seed = 1)
  subset_retina <- FindVariableGenes(object = subset_retina, do.plot = T, display.progress = F, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = 0.5)
  print(length(x=subset_retina@var.genes))
  assign(paste("subset_retina",i,sep="_"), subset_retina)
}

# Determine genes to use for CCA, must be highly variable in at least 2 datasets


ob.list <- list()

for (i in test) {
  vecname<-paste("subset_retina",i,sep="_")
  ob.list[[i]]<-get(vecname)
  }

genes.use <- c()

for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 3000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

length(genes.use)

#Perform linear dimensional reduction
#Perform PCA on the scaled data. By default, the genes in object@var.genes are used as input, but can be defined using pc.genes. We have typically found that running dimensionality reduction on genes with high-dispersion can improve performance. However, with UMI data - particularly after using RegressOut, we often see that PCA returns similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome.


Seurat_object.intregrated <- RunCCA(object=ob.list[1]$WR, object2=ob.list[2]$Cd73ft, group.by=NULL, num.cc = 20, genes.use=genes.use)

    
png(filename="CCA_1v2_before_alignment.Cond_Sorting.png", width=1500, height=700, res=150)
VlnPlot(object = Seurat_object.intregrated, features.plot = "CC1", group.by = "Sorting", 
        do.return = TRUE, x.lab.rot = TRUE)
dev.off()


png(filename="MetageneBicorPlot.Seurat_object.intregrated.png", width=3500, height=1000, res=150)
MetageneBicorPlot(Seurat_object.intregrated, grouping.var = "Sorting", dims.eval = 1:12)
dev.off()

png(filename="DimHeatmap.cca.Seurat_object.intregrated.png", width=1200, height=1700, res=150)
DimHeatmap(object = Seurat_object.intregrated, reduction.type = "cca", cells.use = 100, 
    dim.use = 1:12, do.balanced = TRUE)
dev.off()

#Alignement

Seurat_object.intregrated <- AlignSubspace(object = Seurat_object.intregrated, 
  reduction.type = "cca", grouping.var = "Sorting", 
    dims.align = 1:5)

#Visualize the aligned CCA and perform integrated analysis

p1 <- VlnPlot(object = Seurat_object.intregrated, features.plot = "ACC1", group.by = "Sorting", 
    do.return = TRUE, x.lab.rot = TRUE)
p2 <- VlnPlot(object = Seurat_object.intregrated, features.plot = "ACC2", group.by = "Sorting", 
    do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_after_alignment.Sorting.png", width=2000, height=1000, bg = "white", res = 150)
plot_grid(p1, p2)
dev.off()


## Write a table with the genes that correlate most with each CCA component to find meaningful signals that
## contribute to each CCA
dim_top_genes <- DimTopGenes(Seurat_object.intregrated,
                             reduction.type = "cca", 
                             dim.use = DIM_nb,
                             do.balanced=TRUE)

write.table(file="retina_cca_top_genes.txt",
            dim_top_genes,
            sep="\t",
            col.names=F,
            row.names=F,
            quote=FALSE)



##Rn TSNE

Seurat_object.intregrated <- RunTSNE(object = Seurat_object.intregrated, dims.use = 1:5, do.fast = T, dim_embed=2, perplexity=20, reduction.use = "cca.aligned",
  reduction.name = "tsne_cca", tsne.method = "Rtsne")


##Run UMAP

Seurat_object.intregrated <- RunUMAP(Seurat_object.intregrated, cells.use = NULL, dims.use = 1:5, reduction.use = "cca.aligned",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap_cca", reduction.key = "UMAP", n_neighbors = 30L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)


#######Find cluster

##For CCA
vascendothelium <- FindClusters(object= Seurat_object.intregrated , reduction.type = "cca.aligned", dims.use = c(1:5), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(vascendothelium@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="TSNEPlot-cluster.cca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.cca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="UMAPPlot-cluster.cca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.cca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


##### QC on clusters


png(filename="VlnPlotQC.vascendothelium.pca.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(vascendothelium, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = colors, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

png(filename="Tree_initial.png", res = 150, width=1000, height=1500)
BuildClusterTree(vascendothelium, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()


##Find cell types

## Regrouping

oldidentname <- c("0", "2", "3")
for (clusterID in oldidentname) {
  vascendothelium <- RenameIdent(vascendothelium, old.ident.name = clusterID, new.ident.name = "1")
} 

vascendothelium <-ValidateSpecificClusters(vascendothelium, cluster1 = 5, cluster2 = 1,
  pc.use = 1:6, top.genes = 30, acc.cutoff = 0.9)

png(filename="TSNEPlot-cluster.pca.regrouping.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster.pca.regrouping.png", width=500, height=500, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="Tree_regrouping.png", res = 150, width=1000, height=1500)
BuildClusterTree(vascendothelium, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
  do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
  show.progress = TRUE)
dev.off()

png(filename="VlnPlotQC.vascendothelium.regrouping.png", width=2000, height=2000, bg = "white", res = 150)
VlnPlot(vascendothelium, c("nGene", "nUMI", "percent.mito"), ident.include = NULL, nCol = 2,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
  size.y.use = 16, size.title.use = 20, adjust.use = 1,
  point.size.use = -1, cols.use = NULL, group.by = NULL, y.log = FALSE,
  x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
  single.legend = TRUE, remove.legend = FALSE, do.return = FALSE)
dev.off()

#####

vascendothelium.markers <- FindAllMarkers(object = vascendothelium, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

vascendothelium.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file ="vascendothelium-subclustered-Markers-top10.tsv", sep = "\t")

top10.markers <- unique(top10$gene)

png(filename="DotPlot_vascendothelium_subclustered_top10_markers.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(vascendothelium, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="DotPlot_vascendothelium_subclustered_top10_markers_RetinaMarkers.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(vascendothelium, unique(c(top10.markers, "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "RLBP1", "FBN1", "CLDN5", "EGR1", "RHOB", "LPCAT2", "DOCK2", "CSF3R", "LYZ2", "C1QA", "CX3CR1", "GFAP", "OPTC", "TOP2A", "MKI67")), cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


png(filename="ClusterTree.png", width=1000, height=1000, bg = "white", res = 150)
vascendothelium <- BuildClusterTree(object = vascendothelium, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  

#vascendothelium.markersNode <- FindAllMarkersNode(object =vascendothelium, node = NULL, only.pos = FALSE, min.pct = 0.15, thresh.use = 0.15)

#vascendothelium.markersNode %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10.node

#write.table(top10.node, file ="vascendothelium-subclustered-MarkersNode-top10.avg_logFC.tsv", sep = "\t")

png(filename="FeaturePLot.vascendothelium.TOP2A.MKI67.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(vascendothelium, c("TOP2A", "MKI67"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("aquamarine2", "orange", "pink", "red"), pch.use = 16, overlay = TRUE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, coord.fixed = FALSE,
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()


###Gene set for cluster identification



###Renaming clusters

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "1", new.ident.name = "Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "5", new.ident.name = "Proliferative Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "4", new.ident.name = "Rod/Bip conta")


set.seed(001)

require(scales)

identities <- levels(vascendothelium@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))


png(filename="TSNEPlot-cluster_renamed.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster_renamed.png", width=1000, height=700, bg = "white", res = 150)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="ClusterTree.renamed.png", width=1000, height=1000, bg = "white", res = 150)
vascendothelium <- BuildClusterTree(object = vascendothelium, genes.use = NULL, pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE)
dev.off()  


vascendothelium <- StashIdent(vascendothelium, save.name = "Cell_type")

vascendothelium <- SubsetData(vascendothelium, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Rod/Bip conta"), accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

png("Variable_gene_vascendothelium.renamed.png")
vascendothelium <- FindVariableGenes(object = vascendothelium, do.plot = T, display.progress = F, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.20, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()


vascendothelium <- RunPCA(object = vascendothelium, features = vascendothelium@var.genes)

png("ElbowPlot_pca.renamed.png")
PCElbowPlot(object = vascendothelium)
dev.off()

png(filename="DimHeatmap.pca.renamed.png", width=1200, height=1700, res=150)
DimHeatmap(object = vascendothelium, reduction.type = "pca", cells.use = 100, 
    dim.use = 1:12, do.balanced = TRUE)
dev.off()


vascendothelium <- FindClusters(object= vascendothelium , reduction.type = "pca", dims.use = c(1:6), resolution = 0.8, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 15,  algorithm = 1, prune.SNN=1/15)


vascendothelium <- RunTSNE(object = vascendothelium, dims.use = 1:6, do.fast = T, dim_embed=2, perplexity=15, reduction.use = "pca",
  reduction.name = "tsne", tsne.method = "Rtsne")


vascendothelium <- RunUMAP(vascendothelium, cells.use = NULL, dims.use = 1:6, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap", reduction.key = "UMAP", n_neighbors = 15L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)


png(filename="TSNEPlot-cluster_renamed.final.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster_renamed.final.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


vascendothelium.markers <- FindAllMarkers(object = vascendothelium, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

vascendothelium.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file ="vascendothelium-subclustered-Markers-top10.final.tsv", sep = "\t")

top10.markers <- unique(top10$gene)

png(filename="DotPlot_vascendothelium_subclustered_top10_markers.final.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(vascendothelium, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

vascendothelium <- SubsetData(vascendothelium, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c(5), accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)


vascendothelium<- RunTSNE(object = vascendothelium, dims.use = 1:6, do.fast = T, dim_embed=2, perplexity=15, reduction.use = "pca",
  reduction.name = "tsne", tsne.method = "Rtsne")


vascendothelium <- RunUMAP(vascendothelium, cells.use = NULL, dims.use = 1:6, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap", reduction.key = "UMAP", n_neighbors = 15L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)



vascendothelium.markers <- FindAllMarkers(object = vascendothelium, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

vascendothelium.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file ="vascendothelium-subclustered-Markers-top10.final.filtered.tsv", sep = "\t")

top10.markers <- unique(top10$gene)

png(filename="DotPlot_vascendothelium_subclustered_top10_markers.final.filtered.png", width=2000, height=1000, bg = "white", res = 150)
DotPlot(vascendothelium, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png(filename="DotPlot_VascularEndothelium_subclustered_TipStalkMarkers.png", width=1000, height=600, bg = "white", res = 150)
DotPlot(vascendothelium, c("APOLD1", "BMX", "PDGFB", "ID1", "CXCR4", "APLNR", "FLT1", "LFNG", "MFNG", "RFNG", "CLDN5", "UNC5B", "EFNB2", "NOTCH1", "HEY1", "HES1", "JAG1", "JAG2", "KDR", "CTHRC1", "FLT4", "DLL3", "DLL4", "ESM1", "ANGPT2", "APLN", "GUF1", "TOP2A"), cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()



vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "0", new.ident.name = "Quiescent Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "6", new.ident.name = "Proliferative Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "1", new.ident.name = "TM4SF1 Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "2", new.ident.name = "APLN Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "3", new.ident.name = "SKIL Endothelial cells")

vascendothelium <- RenameIdent(vascendothelium, old.ident.name = "4", new.ident.name = "LY6C1 Endothelial cells")


vascendothelium <- StashIdent(vascendothelium, save.name = "Cell_type")


png(filename="DotPlot_vascendothelium_subclustered_top10_markers.final.filtered.renamed.png", width=2500, height=1000, bg = "white", res = 150)
DotPlot(vascendothelium, top10.markers, cols.use = c("blue", "red"),
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()



png(filename="TSNEPlot-cluster_renamed.final.filtered.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 2, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster_renamed.final.filtered.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 2, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-cluster_renamed.final.filtered.legend.png", width=1000, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster_renamed.final.filtered.legend.png", width=1000, height=700, bg = "white", res = 200)
DimPlot(vascendothelium, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()


saveRDS(vascendothelium, "vascendothelium_subclustered_renamed.rds")

### Make distribution plot for P12 to P17


setwd(file.path(B_DIR, "Subclustering/VascularEndothelium/Mapping"))

vascendothelium <- readRDS("vascendothelium_subclustered_renamed.rds")

batch_long <- vascendothelium@meta.data
batch_assigned <- batch_long 
batch_assigned <- batch_assigned %>%
   mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
batch_assigned <- as.data.frame(batch_assigned[,"Cond_TimePoint"])
 
rownames(batch_assigned) <- rownames(batch_long)
 
colnames(batch_assigned) <- "Cond_TimePoint"
 
vascendothelium <- AddMetaData(object = vascendothelium, metadata = batch_assigned, col.name = "Cond_TimePoint")

saveRDS(vascendothelium, "vascendothelium_subclustered_renamed.rds")




cluster_counts <- vascendothelium@meta.data %>%
  group_by(Cond_TimePoint) %>%
  count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(Cond_TimePoint,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(cluster_portions,file="Cluster_proportions_CellType.VascEndoCells.Cond_TimePoint.png", width = 8, height = 5)

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.VascEndoCells.Cond_TimePoint.png")


cluster_counts <- vascendothelium@meta.data %>%
  group_by(Cell_type) %>%
  count(Cond_TimePoint) %>%
  mutate(freq = n / sum(n))

cluster_counts$Cond_TimePoint <- factor(cluster_counts$Cond_TimePoint,levels=unique(mixedsort(cluster_counts$Cond_TimePoint)))


cluster_portions <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=Cell_type)) +
  geom_bar(mapping = aes(reorder(Cell_type, freq),freq,fill=Cell_type), stat="identity",col="black") +
  facet_wrap(~ Cond_TimePoint, scale="free") +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(cluster_portions,file="Cluster_proportions_CellType.VascEndoCells.Cell_type.png", width = 10, height = 7)


##Dimentionality reduction for P17 time point


Subseted_cells_P17 <- rownames(subset(vascendothelium@meta.data, TimePoint %in% c("P17")))

vascendothelium_P17 <- SubsetData(vascendothelium, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


vascendothelium_P17 <- RunTSNE(object = vascendothelium_P17, dims.use = 1:3, do.fast = T, dim_embed=2, perplexity=15, reduction.use = "pca",
  reduction.name = "tsne", tsne.method = "Rtsne")


vascendothelium_P17 <- RunUMAP(vascendothelium_P17, cells.use = NULL, dims.use = 1:6, reduction.use = "pca",
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap", reduction.key = "UMAP", n_neighbors = 15L,
  min_dist = 0.5, metric = "correlation", seed.use = 42)


png(filename="TSNEPlot-cluster_renamed.final.P17.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium_P17, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UmapPlot-cluster_renamed.final.P17.png", width=700, height=700, bg = "white", res = 200)
DimPlot(vascendothelium_P17, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()




###Make Single R Object

singler = CreateSinglerObject(microglia@data, annot = NULL, "Microglia", min.genes = 100,
  technology = "Dropseq", species = "Mouse", citation = "unpublished",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T, 
  reduce.file.size = T, numCores = SingleR.numCores)

singler$seurat = microglia # (optional)
singler$meta.data$orig.ident = microglia@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = microglia@dr$tsne_cca@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = microglia@ident # the Seurat clusters (if 'clusters' not provided)

save(singler,file='singler_microglia.RData')

load(file='singler_microglia.RData')


png("SingleR.PlotTsne.png", width=700)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
      singler$meta.data$xy, do.label = FALSE, do.letters = F,
      labels=singler$meta.data$clusters,label.size = 4, 
      dot.size = 3)
out$p
dev.off()


png("SingleR.DrawHeatmap.main.png", width=3000, height=1000, res=150)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$clusters)
dev.off()

png("SingleR.DrawHeatmap.png", width=3000, height=1500, res=150)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,
                    clusters = singler$meta.data$clusters)
dev.off()

png("SingleR.DrawHeatmap.2.main.png", width=3000, height=1000, res=150)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$clusters)
dev.off()

png("SingleR.DrawHeatmap.2.png", width=3000, height=1500, res=150)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50,
                    clusters = singler$meta.data$clusters)
dev.off()


png("SingleR.PlotTsne.labels.png", width=2500, height=2200, res=300)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
        singler$meta.data$xy, do.label=FALSE,
        do.letters = F,labels = singler$singler[[1]]$SingleR.single$labels,
        label.size = 4, dot.size = 3)
out$p
dev.off()

png("SingleR.PlotTsne.labels.main.png", width=2500, height=2200, res=300)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,
        singler$meta.data$xy, do.label=FALSE,
        do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
        label.size = 4, dot.size = 3)
out$p
dev.off()


png("SingleR.PlotTsne.2_labels.png", width=1800, height=1700, res=300)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
        singler$meta.data$xy,do.label=FALSE,
        do.letters = F,labels = singler$singler[[2]]$SingleR.single$labels,
        label.size = 4, dot.size = 3)
out$p
dev.off()

png("SingleR.PlotTsne.2_labels.main.png", width=1800, height=1700, res=300)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,
        singler$meta.data$xy,do.label=FALSE,
        do.letters = F,labels = singler$singler[[2]]$SingleR.single.main$labels,
        label.size = 4, dot.size = 3)
out$p
dev.off()


microglia <- AddMetaData(object = microglia, metadata = singler$singler[[2]]$SingleR.single$labels, col.name = "SingleR.single.labels")


Subseted_cells_P14_P17 <- rownames(subset(microglia@meta.data, TimePoint %in% c("P14", "P17")))

microglia.P14_P17 <- SubsetData(microglia, cells.use = Subseted_cells_P14_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


cluster_counts <- microglia.P14_P17@meta.data %>%
  group_by(Condition) %>%
  count(SingleR.single.labels) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$SingleR.single.labels <- factor(cluster_counts$SingleR.single.labels,levels=unique(mixedsort(cluster_counts$SingleR.single.labels)))

cluster_portions <-  ggplot(cluster_counts,aes(Condition,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ SingleR.single.labels, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.SingleR.single.labels.2_ImmuneCells.P14_P17.Condition.png", width = 10, height = 10)

cluster_portions_full <-  ggplot(cluster_counts,aes(SingleR.single.labels,freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.SingleR.single.labels.2.ImmuneCells.P14_P17.Condition.png")



cluster_counts <- microglia.P14_P17@meta.data %>%
  group_by(Cond_TimePoint) %>%
  count(SingleR.single.labels) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$SingleR.single.labels <- factor(cluster_counts$SingleR.single.labels,levels=unique(mixedsort(cluster_counts$SingleR.single.labels)))

cluster_portions <-  ggplot(cluster_counts,aes(Cond_TimePoint,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ SingleR.single.labels, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.microglia.SingleR.single.labels.2.P14_P17.Cond_TimePoint.png", width = 12, height = 10)

write.table(cluster_counts, "cluster_counts.CellType.microglia.SingleR.single.labels.2.P14_P17.Cond_TimePoint.txt")

cluster_portions_full <-  ggplot(cluster_counts,aes(SingleR.single.labels,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity", col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.microglia.SingleR.single.labels.2.P14_P17.Cond_TimePoint.png")



###Some bash script
filenames <- list.files(path=getwd()) 

#filenames <- "GSM3351760_PV_Precl4d_2_4.coutt.csv"

dataframe_Stich <- fread("GSM3351760_PV_Precl4d_2_4.coutt.csv",sep="\t",header=T)
setkey(dataframe_Stich, by="GENEID")
colnames(dataframe_Stich) <- paste(i, colnames(dataframe_Stich), sep = "_")
colnames(dataframe_Stich)[1] <- "gene"

dataframe_Stich <- dataframe_Stich[,1]

for (i in filenames){
  i <- "GSM3351760_PV_Precl4d_2_4.coutt.csv"
  expression_matrix <- fread(i,sep="\t",header=T)
  setkey(expression_matrix, by="GENEID") 
  colnames(expression_matrix) <- paste(i, colnames(expression_matrix), sep = "_")
  colnames(expression_matrix)[1] <- "gene"
  dataframe_Stich <- full_join(dataframe_Stich, expression_matrix, by="gene")
}


head(dataframe_Stich[1:5,1:5])
dim(dataframe_Stich)

dataframe_Stich[is.na(dataframe_Stich)] = 0

head(dataframe_Stich[1:5,1:5]) 

write.table(dataframe_Stich, file="Neuroinflammation.txt", quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)


#END
q("no")
