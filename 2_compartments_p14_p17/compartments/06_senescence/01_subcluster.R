# ---------------------------------------------------------------------------
# 06_senescence - 01_subcluster.R
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# the senescence-merged compartment, run on cells taken from the annotated
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
library(ktplots)

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/Mapping"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

retina <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.subseted.rds", sep="."))

summary(as.factor(retina@meta.data$Cond_Sorting))


###Remove timepoint

retina <- StashIdent(retina, save.name = "cell_subtype")

retina <- SetAllIdent(retina, id = "TimePoint")

retina <- SubsetData(retina, cells.use = NULL, subset.name = NULL, ident.use = NULL,
  ident.remove = c("P5", "P7", "P10", "P12"))


####Subclustering

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/Mapping"))

retina <- SetAllIdent(retina, id = "cell_type")


markers.retina.dotplot <- rev(c("RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "RLBP1", "FBN1", "CLDN5", "LYZ2", "GFAP", "OPTC", "TOP2A", "MKI67"))

png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.initial.png", sep="."), res = 150, width=1000, height=1000)
DotPlot(retina, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()

exp_mat <- tryCatch(retina@data, error = function(e) {
            tryCatch(GetAssayData(object = retina), error = function(e) {
                stop(paste0("are you sure that your data is normalized?"))
return(NULL)
            })
        })


png(filename=paste(project_name, res, Dim, perp, "geneDotPlot_markers.test.png", sep="."), res = 300, width=2000, height=1000)
exp_mat <- retina@data
metadata <- retina@meta.data
geneDotPlot(retina, # object 
  idents = "cell_type", # column name in meta data that holds the cell-cluster ID/assignment
  genes = c("GFAP", "KCNJ8", "RLBP1", "LYZ2"), # genes to plot
  split.by = "Condition", # column name in the meta data that you want to split the plotting by. If not provided, it will just plot according to idents
  save.plot = FALSE,
  heat_cols = rev(RColorBrewer::brewer.pal(9, "RdBu"))) # If TRUE, it will save to a location that you can specify via filepath and filename
dev.off()


set.seed(001)

require(scales)

identities <- levels(retina@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="UMAPPlot-Subclustered_retina.png", width=500, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TsnePlot-Subclustered_retina.png", width=500, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "tsne", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()




#Find markers before subclusterisation:

retina <- ScaleData(object = retina, vars.to.regress = c("nUMI", "percent.mito", "Batch", "percent.crystal"))

#png(filename="MeanVarPlot.png")
#retina <- FindVariableGenes(object = retina, mean.function = ExpMean, dispersion.function = LogVMR, 
#    x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.5)
#dev.off()

#length(x=retina@var.genes)


###Subset the different dataset

test <- unique(retina@meta.data$Sorting)

for (i in test) {
  subset_cells <- rownames(subset(retina@meta.data, Sorting %in% i))
  subset_retina <- SubsetData(retina, cells.use = subset_cells, subset.name = NULL, ident.use = NULL,
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

svg(filename="DimHeatmap.cca.Seurat_object.intregrated.svg", width=12, height=17)
DimHeatmap(object = Seurat_object.intregrated, reduction.type = "cca", cells.use = 50, 
    dim.use = 1:9, do.balanced = TRUE)
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
retina <- FindClusters(object= Seurat_object.intregrated , reduction.type = "umap_cca", dims.use = c(1:2), resolution = 1, 
  print.output = 0, save.SNN = T, temp.file.location = getwd(), force.recalc=TRUE, k.param = 20,  algorithm = 1)

#, prune.SNN=0.1)


#Make plots

set.seed(001)

require(scales)

identities <- levels(retina@ident)

colors <- rainbow(length(identities))

#colors <- terrain.colors(length(unique(Seurat_object@ident)))

colors <- viridis(length(identities))

require(scales)

# Create vector of default ggplot2 colors
colors <- sample(hue_pal(h = c(0, 360), c = 150, l = 60, h.start = 0, direction = -1)(length(identities)))



png(filename="TSNEPlot-cluster.cca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "tsne_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="TSNEPlot-Cond_Sorting.cca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "tsne_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



png(filename="TSNEPlot-cell_subtype.cca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "tsne_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "cell_subtype", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()



graphplot <- DimPlot(retina, reduction.use = "tsne_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = TRUE, do.bare = FALSE,
  cols.use = colors, group.by = "cell_subtype", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)

ggsave(graphplot, filename="TSNEPlot-cell_subtype.cca.eps", width=7, height=5, bg = "white", dpi = 300)


png(filename="UMAPPlot-cluster.cca.png", width=500, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "umap_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

png(filename="UMAPPlot-Cond_Sorting.cca.png", width=800, height=500, bg = "white", res = 150)
DimPlot(retina, reduction.use = "umap_cca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = colors, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE)
dev.off()

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png(filename="UMAPPlot-cell_subtype.cca.png", width=800, height=500, bg = "white", res = 150)
#~ DimPlot(retina, reduction.use = "umap_cca", dim.1 = 1, dim.2 = 2,
#~   cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
#~   cols.use = colors, group.by = "cell_subtype", pt.shape = NULL,
#~   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#~   do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
#~   dark.theme = FALSE)
#~ dev.off()
#~
#~ ##############################
#~
#~
#~ ### Add GSVA to seurat object assay
#~
#~ #GO_score=data.frame(fread("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina//OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/GSVA/retina_senescence_merged_data_gsva_out_April2020.csv", sep=",", header=TRUE), row.names=1)
#~
#~ GO_score= read.table("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina//OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/GSVA/retina_senescence_merged_data_gsva_out_April2020.csv", sep = ",", header = T, row.names=1, stringsAsFactors=F)
#~
#~ Seurat_object_subset <- SubsetData(Seurat_object.P17, cells.use = colnames(GO_score), subset.name =NULL, ident.use = NULL,
#~   ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
#~   do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
#~   random.seed = 1)
#~
#~ GO_score <- GO_score[,colnames(retina@data)]
#~
#~ retina <- SetAssayData(retina, assay.type = "GO", slot = "raw.data", new.data = GO_score)
#~
#~ retina <- NormalizeData(retina, assay.type = "GO", normalization.method = "genesCLR")
#~
#~ retina <- ScaleData(retina, assay.type = "GO", display.progress = FALSE, model.use = "negbinom")
#~
#~ retina@assay$GO@data <- retina@assay$GO@scale.data
#~
#~ retina <- SetAllIdent(retina, id="cell_subtype")
#~
#~
#~ ######### 
#~
#~ ####Gating on  enrichment score for specific function
#~
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/Mapping/Gating")
#~
#~
#~ ##Ridgeplot 
#~ grep.genes <- grep("NEUTROPHIL_BINET", rownames(GO_score))
#~
#~ grep.genes <- grep("HALLMARK_KRAS_SIGNALING_UP", rownames(GO_score))
#~
#~ rownames(GO_score[grep.genes,])
#~
#~ #pathway <- "WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1"
#~
#~ pathway <- rownames(GO_score[grep.genes,])
#~
#~ retina <- SetAllIdent(retina, id="cell_subtype")
#~
#~ retina <- ReorderIdent(retina, feature = pathway, rev = FALSE, aggregate.fxn = mean,
#~     reorder.numeric = FALSE)
#~
#~ png(paste(pathway,"Ridgeplot.png",sep="_"))
#~ RidgePlot(retina, pathway, do.return = TRUE, , size.title.use = 10, do.sort = FALSE)
#~ dev.off()


png(paste(pathway,"Ridgeplot.subset_Senescent.png",sep="_"))
RidgePlot(retina, pathway, do.return = TRUE, size.title.use = 10, do.sort = FALSE, ident.include = c("Astrocytes", "Endothelial cells", "Muller glia", "Pericytes"))
dev.off()

graphplot <- RidgePlot(retina, pathway, do.return = TRUE, size.title.use = 10, do.sort = FALSE, ident.include = c("Astrocytes", "Endothelial cells", "Muller glia", "Pericytes"))
ggsave(graphplot, filename=paste(pathway,"Ridgeplot.subset_Senescent.eps",sep="_"), width=7, height=5, bg = "white", dpi = 300)

png(paste(pathway,"Ridgeplot.subset_Endothelial_Senescent.png", sep="_"))
RidgePlot(retina, pathway, do.return = TRUE, size.title.use = 10, do.sort = FALSE, ident.include = c("Endothelial cells"), cols=viridis(1, dir=-1))
dev.off()

graphplot <- RidgePlot(retina, pathway, do.return = TRUE, size.title.use = 10, do.sort = FALSE, ident.include = c("Endothelial cells"), cols=viridis(1, dir=-1))
ggsave(graphplot, filename=paste(pathway,"Ridgeplot.subset_Endothelial_Senescent.eps",sep="_"), width=7, height=5, bg = "white", dpi = 300)


png(paste(pathway,"Ridgeplot.subset_immune.png",sep="_"))
RidgePlot(retina, pathway, do.return = TRUE, , size.title.use = 10, do.sort = FALSE, 
  ident.include = c("Immune cells 1", "Immune cells 2", "Immune cells 3", "Immune cells 4", "Immune cells 5"))
dev.off()

#gate

retina <- SetAllIdent(retina, id="cell_subtype")

Seurat_object_gated_cells <- WhichCells(retina, ident = c("Astrocytes", "Endothelial cells", "Muller glia", "Pericytes"), ident.remove = NULL, cells.use = NULL,
  subset.name = "HALLMARK_KRAS_SIGNALING_UP", accept.low = 0.5, accept.high = Inf,
  accept.value = NULL, max.cells.per.ident = Inf, random.seed = 1)

retina <- SetIdent(retina, cells.use = Seurat_object_gated_cells, ident.use = "KRAS positive ")

Seurat_object_gated_cells <- WhichCells(retina, ident = c("Astrocytes", "Endothelial cells", "Muller glia", "Pericytes"), ident.remove = NULL, cells.use = NULL,
  subset.name = "HALLMARK_KRAS_SIGNALING_UP", accept.low = -Inf, accept.high = 0.5,
  accept.value = NULL, max.cells.per.ident = Inf, random.seed = 1)

retina <- SetIdent(retina, cells.use = Seurat_object_gated_cells, ident.use = "KRAS negative ")

Seurat_object_gated_cells <- WhichCells(retina, ident = c("Immune cells 1", "Immune cells 2", "Immune cells 3", "Immune cells 4", "Immune cells 5"), ident.remove = NULL, cells.use = NULL,
  subset.name = "HALLMARK_KRAS_SIGNALING_UP", accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, max.cells.per.ident = Inf, random.seed = 1)

retina <- SetIdent(retina, cells.use = Seurat_object_gated_cells, ident.use = "")

retina <- StashIdent(retina, save.name = "KRAS_cells")

df <- as.data.frame(paste(retina@meta.data$KRAS_cells,retina@meta.data$cell_subtype, sep=""))

rownames(df) <- rownames(retina@meta.data)

colnames(df) <- "KRAS_cell_subtype"

retina <- AddMetaData(retina, df, col.name = "KRAS_cell_subtype")

retina <- SetAllIdent(retina, id = "KRAS_cell_subtype")

###
grep.genes <- grep("HALLMARK_KRAS_SIGNALING_UP", rownames(GO_score))

rownames(GO_score[grep.genes,])

pathway <- rownames(GO_score[grep.genes,])

retina <- ReorderIdent(retina, feature = pathway, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)

png(paste(pathway,"Ridgeplot.KRAS_cells.png",sep="_"))
RidgePlot(retina, pathway, do.return = TRUE, , size.title.use = 10, do.sort = FALSE)
dev.off()

png(paste(pathway,"Ridgeplot.subset_Senescent.KRAS_cells.png",sep="_"))
RidgePlot(retina, pathway, do.return = TRUE, , size.title.use = 10, do.sort = FALSE, 
                              ident.include = c("KRAS positive Astrocytes", "KRAS positive Endothelial cells", 
                                                "KRAS positive Muller glia", "KRAS positive Pericytes",
                                                "KRAS negative Astrocytes", "KRAS negative Endothelial cells", 
                                                "KRAS negative Muller glia", "KRAS negative Pericytes"))
dev.off()

graphplot <- RidgePlot(retina, pathway, do.return = TRUE, , size.title.use = 10, do.sort = FALSE, 
                              ident.include = c("KRAS positive Astrocytes", "KRAS positive Endothelial cells", 
                                                "KRAS positive Muller glia", "KRAS positive Pericytes",
                                                "KRAS negative Astrocytes", "KRAS negative Endothelial cells", 
                                                "KRAS negative Muller glia", "KRAS negative Pericytes"))

ggsave(graphplot, filename=paste(pathway,"Ridgeplot.subset_Senescent.KRAS_cells.eps",sep="_"), width=7, height=5, bg = "white", dpi = 300)


png(paste(pathway,"Ridgeplot.subset_immune.KRAS_cells.png",sep="_"))
RidgePlot(retina, pathway, do.return = TRUE, , size.title.use = 10, do.sort = FALSE, 
  ident.include = c("Immune cells 1", "Immune cells 2", "Immune cells 3", "Immune cells 4", "Immune cells 5"))
dev.off()




saveRDS(retina, "../retina_subclustered_renamed.rds")

retina <- readRDS("retina_subclustered_renamed.rds")





#END
q("no")
