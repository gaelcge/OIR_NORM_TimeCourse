# ---------------------------------------------------------------------------
# Step 02b - Clustering of the retained and the discarded sets.
#
# ESTABLISHES
#     Two objects: the CCA-filtered set, and the discarded set with its own
#     PCA, t-SNE and UMAP so that what was removed can be seen rather than
#     assumed. Step 04 reads the discarded one.
#
# WHY THIS WAY
#     The discarded cells are embedded and plotted for PECAM1, CLDN5 and MKI67
#     - endothelial and proliferation markers. The question being asked is
#     whether the variance-ratio filter preferentially removes a real
#     population rather than low-quality droplets, which is the failure mode
#     that would make the filter unsafe.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 01b
#
# INPUTS
#     <B_DIR>/Clustering/<project_name>.Seurat_object.rds
#
# OUTPUTS
#     <B_DIR>/Clustering/discarded/Seurat_object.integrated.discard_<Dim>CC.rds
#
# USAGE
#     Rscript 02b_*.R
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
library(harmony)

# Load the dataset

setwd(file.path(B_DIR, "Clustering/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30


Seurat_object.intregrated <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))


#Before we align the subspaces, we first search for cells whose expression profile cannot be well-explained by low-dimensional CCA, compared to low-dimensional PCA.


Seurat_object.intregrated <- CalcVarExpRatio(object = Seurat_object.intregrated, reduction.type = "pca", grouping.var = "Cond_Sorting", 
    dims.use = DIM_nb)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
Seurat_object.integrated.all.save <- Seurat_object.intregrated
Seurat_object.integrated <- SubsetData(object = Seurat_object.integrated.all.save, subset.name = "var.ratio.pca", 
  accept.low = 0.5)

#For illustrative purposes, we can also look at the discarded cells. Note that discarded cells tend to have lower gene counts, but a subset also express high levels of PF4. This is because megakaryocytes are only present in the 10X dataset, and thus are correctly discarded as ‘dataset-specific’]

setwd(file.path(B_DIR, "Clustering/discarded"))

Seurat_object.integrated.discard <- SubsetData(object = Seurat_object.integrated.all.save, subset.name = "var.ratio.pca", 
    accept.high = 0.5)
median(x = Seurat_object.integrated@meta.data[, "nGene"])

median(x = Seurat_object.integrated.discard@meta.data[, "nGene"])


png(paste("Discarded.cells.",Dim,".png",sep=""), width=2000, height=500, bg = "white", res = 150)
VlnPlot(object = Seurat_object.integrated.discard, features.plot = c("PECAM1", "CLDN5", "MKI67"), group.by = "Cond_Sorting", x.lab.rot = TRUE)
dev.off()

png(filename=paste("VlnPlotQC", project_name, res, Dim, perp, "Cond_Sorting.Discarded.cells.png", sep="."), width=2500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object.integrated.discard, c("nGene", "nUMI", "percent.mito"), group.by = "Cond_Sorting", nCol = 2, x.lab.rot = TRUE)
dev.off()


table(Seurat_object.integrated.discard@meta.data$Cond_Sorting)


Seurat_object.integrated.discard <- RunPCA(object = Seurat_object.integrated.discard, pc.genes = Seurat_object.integrated.discard@var.genes, do.print = TRUE, pcs.print = 1:10, 
    genes.print = 10)

Seurat_object.integrated.discard <- ProjectPCA(object = Seurat_object.integrated.discard, do.print = FALSE)


Seurat_object.integrated.discard <- RunTSNE(Seurat_object.integrated.discard,
                               reduction.use = "pca",
                               dims.use = DIM_nb, do.fast = T, dim_embed=2, perplexity=perp)

#
png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "Cond_Sorting.discard.pca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
TSNEPlot(object =Seurat_object.integrated.discard, do.return = T, no.legend = F, do.label = F, group.by = "Cond_Sorting")
dev.off()


Seurat_object.integrated.discard <- RunUMAP(Seurat_object.integrated.discard, reduction.use = "pca", dims.use = DIM_nb)


png(filename=paste("UmapPlot-cluster.",Dim,".Cond_Sorting.discard.pca.png",sep=""), width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object.integrated.discard, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()


saveRDS(Seurat_object.integrated.discard, paste("Seurat_object.integrated.discard_",Dim,"CC.rds",sep=""))


#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

setwd(file.path(B_DIR, "Clustering/"))


png(filename="MetageneBicorPlot.retina.integrated.overlapped.png", width=2000, height=1000, bg = "white", res = 150)
MetageneBicorPlot(Seurat_object.integrated, grouping.var = "Cond_Sorting", dims.eval = 1:30, 
                  display.progress = FALSE)
dev.off()                  


#Alignement

Seurat_object.integrated <- AlignSubspace(object = Seurat_object.integrated, 
  reduction.type = "cca", grouping.var = "Cond_Sorting", 
    dims.align = DIM_nb)

#Visualize the aligned CCA and perform integrated analysis

p1 <- VlnPlot(object = Seurat_object.integrated, features.plot = "ACC1", group.by = "Cond_Sorting", 
    do.return = TRUE, x.lab.rot = TRUE)
p2 <- VlnPlot(object = Seurat_object.integrated, features.plot = "ACC2", group.by = "Cond_Sorting", 
    do.return = TRUE, x.lab.rot = TRUE)

png("CC_1v2_after_alignment.Cond_Sorting.png", width=2000, height=1000, bg = "white", res = 150)
plot_grid(p1, p2)
dev.off()

## Write a table with the genes that correlate most with each CCA component to find meaningful signals that
## contribute to each CCA
dim_top_genes <- DimTopGenes(Seurat_object.integrated,
                             reduction.type = "cca", 
                             dim.use = DIM_nb,
                             do.balanced=TRUE)

write.table(file="retina_cca_top_genes.txt",
            dim_top_genes,
            sep="\t",
            col.names=F,
            row.names=F,
            quote=FALSE)

##Find cluster


Seurat_object.integrated <- FindClusters(Seurat_object.integrated, reduction.type = "cca.aligned",
                                    dims.use = DIM_nb, save.SNN = T, resolution = res, temp.file.location = getwd(), force.recalc=TRUE)


PrintFindClustersParams(object = Seurat_object.integrated)

#Run Rtsne 
#t-SNE method does not require the removal of duplicates. The fact that it is a default feature in Rtsne does not imply its requirement. It is useful for some short-term event monitoring. For characterising long-term trends and/or patterns with big data sets, I see little use. The Rtsne default setup can be more inclined for characterising events in the time-domain, without any studies in Fourier domain.
#Assume you have points in the time-domain. The duplicate algorithm causes significant amount of false positives because the duplicate checking is mostly designed on the time-domain signal. Fourier space can show that those events which are considered by the algorithm duplicate are not necessary so.
#So my observation is that the algorithm is greedy about duplicate points in the time-domain, which is not useful for me when considering long-term signals, long-terms trends and long-term patterns. The fact that the point is duplicate in the time-domain does not actually mean that it is duplicate also in Fourier domain. I think it will be more a coincidence if is a duplicate in a time domain in the real-life applications. So turning off the feature, should be ok. To estimate how much of the points are really duplicates in both domains is specific on the case study. I get significantly better descriptors of events and/or phenomena by considering long-term data sets without the duplicate check in many real-life applications.
#I think the Rtsne documentation is not clear about the case in saying [turn off check_duplicates and] don't wast processing power. There are really other reasons as described above why the check_duplicates can be turned off as realised also by some other implementations of the method. The check_duplicates=TRUE is a personal selection of the Rtsne developer by default at the moment. I would love to hear if there is any implementation reasons for the decision. 
#dim_embed {The dimensional space of the resulting tSNE embedding (default is 2). For example, set to 3 for a 3d tSNE}
#item{\dots}{Additional arguments to the tSNE call. Most commonly used is perplexity (expected number of neighbors default is 30)}
#



###Run mapping algo


Seurat_object.integrated <- RunTSNE(Seurat_object.integrated,
                               reduction.use = "cca.aligned",
                               dims.use = DIM_nb, do.fast = T, dim_embed=2, perplexity=perp)


#QC on clustering
#
png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "ident.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
TSNEPlot(object = Seurat_object.integrated, do.return = T, no.legend = T, do.label = T)
dev.off()
#
png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "Cond_Sorting.cca.png", sep="."), width=1500, height=1000, bg = "white", res = 150)
TSNEPlot(object =Seurat_object.integrated, do.return = T, no.legend = F, do.label = F, group.by = "Cond_Sorting")
dev.off()

#
png(filename=paste("VlnPlotQC", project_name, res, Dim, perp, "Cond_Sorting.regressed.png", sep="."), width=1500, height=2000, bg = "white", res = 50)
VlnPlot(object = Seurat_object.integrated, c("nGene", "nUMI", "percent.mito"), group.by = "Cond_Sorting", nCol = 2, x.lab.rot = TRUE)
dev.off()

png(filename=paste("VlnPlotQC", project_name, res, Dim, perp, "Ident.regressed.png", sep="."), width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object.integrated, c("nGene", "nUMI", "percent.mito"), group.by = "ident", nCol = 2, x.lab.rot = TRUE)
dev.off()


Seurat_object.integrated <- RunUMAP(Seurat_object.integrated, reduction.use = "cca.aligned", dims.use = DIM_nb)

png(filename=paste("UmapPlot", project_name, res, Dim, perp, "ident.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename=paste("UmapPlot", project_name, res, Dim, perp, "Cond_Sorting.cca.png", sep="."), width=1500, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

markers.retina.dotplot <- rev(c("RHO", "LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "OTX2", "PRDM1", "RLBP1", "GFAP", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL"))


png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.cca.png", sep="."), res = 150, width=1000, height=2000)
DotPlot(Seurat_object.integrated, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
  plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
dev.off()


# #Run harmony (optional)

# Seurat_object.integrated@dr$pca <- Seurat_object.integrated@dr$cca.aligned

# Seurat_object.integrated <- RunHarmony(Seurat_object.integrated, "Cond_Sorting", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100)



# Seurat_object.integrated <- FindClusters(Seurat_object.integrated, reduction.type = "harmony",
#                                     dims.use = DIM_nb, save.SNN = T, resolution = res, temp.file.location = getwd(), force.recalc=TRUE)



# #Run mapping algo

# Seurat_object.integrated <- RunUMAP(Seurat_object.integrated, reduction.use = "harmony", dims.use = DIM_nb)

# png(filename=paste("UmapPlot", project_name, res, Dim, perp, "ident.cca.harmony.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
# DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
#   cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
#   cols.use = NULL, group.by = "ident", pt.shape = NULL,
#   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#   do.label = TRUE, label.size = 4, no.legend = TRUE,
#   coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
#   plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
#   sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
# dev.off()

# png(filename=paste("UmapPlot", project_name, res, Dim, perp, "Cond_Sorting.cca.harmony.png", sep="."), width=1500, height=1000, bg = "white", res = 150)
# DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
#   cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
#   cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
#   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#   do.label = FALSE, label.size = 4, no.legend = FALSE,
#   coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
#   plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
#   sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
# dev.off()

# png(filename=paste("UmapPlot", project_name, res, Dim, perp, "Sorting.cca.harmony.png", sep="."), width=1500, height=1000, bg = "white", res = 150)
# DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
#   cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
#   cols.use = NULL, group.by = "Sorting", pt.shape = NULL,
#   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#   do.label = FALSE, label.size = 4, no.legend = FALSE,
#   coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
#   plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
#   sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
# dev.off()

# markers.retina.dotplot <- rev(c("LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "OTX2", "PRDM1", "RLBP1", "GFAP", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL"))


# png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.cca.harmony.png", sep="."), res = 150, width=1000, height=2000)
# DotPlot(Seurat_object.integrated, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
#   col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
#   plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
# dev.off()

#
#
#doublet detection


Seurat_object <- FindVariableGenes(object = Seurat_object.integrated, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = 0.5)

Seurat_object <- RunPCA(object = Seurat_object, pc.genes = Seurat_object@var.genes, do.print = TRUE, pcs.print = 1:10, 
    genes.print = 10)

Seurat_object <- ProjectPCA(object = Seurat_object, do.print = FALSE)

Seurat_object <- FindClusters(object= Seurat_object , reduction.type = "pca", dims.use = DIM_nb, 
  resolution = res, print.output = 0, save.SNN = T, temp.file.location = getwd())

cell_count_Seurat_object <- length(rownames(Seurat_object@meta.data))

Seurat_object <- doubletFinder(Seurat_object, expected.doublets = 0.05*cell_count_Seurat_object, proportion.artificial = 0.25, proportion.NN = 0.01)

doublet <- select(Seurat_object@meta.data, pANNPredictions)

Seurat_object.integrated <- AddMetaData(object = Seurat_object.integrated, metadata = doublet, col.name = "pANNPredictions")


png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "pANNPredictions.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
TSNEPlot(object =Seurat_object.integrated, do.return = T, no.legend = F, do.label = F, group.by = "pANNPredictions", pt.size = 0.1)
dev.off()

png(filename=paste("UmapPlot", project_name, res, Dim, perp, "pANNPredictions.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = c("red", "blue"), group.by = "pANNPredictions", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()


#Remove doublets


#Singlet_cells_names <- rownames(subset(Seurat_object.integrated@meta.data, pANNPredictions %in% c("Singlet")))

#Seurat_object.integrated <- SubsetData(Seurat_object.integrated, cells.use = Singlet_cells_names)


##Check Marker expression


# markers.retina.dotplot <- rev(c("LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "OTX2", "PRDM1", "RLBP1", "GFAP", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL"))

# png(filename=paste(project_name, res, Dim, perp, "ClusterTree_renamed.final.cca.png", sep="."), res = 150, width=1000, height=1500)
# Seurat_object.integrated <- BuildClusterTree(Seurat_object.integrated, genes.use = NULL, pcs.use = NULL, SNN.use = NULL,
#   do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE,
#   show.progress = TRUE)
# dev.off()

# png(filename=paste(project_name, res, Dim, perp, "DotPlot_markers.renamed.final.cca.png", sep="."), res = 150, width=1000, height=2000)
# DotPlot(Seurat_object.integrated, genes.plot= markers.retina.dotplot, cols.use = c("blue", "red"),
#   col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8,
#   plot.legend = TRUE, do.return = FALSE, x.lab.rot = TRUE)
# dev.off()

#Make final plot

#Seurat_object <- RunTSNE(object = Seurat_object, dims.use = DIM_nb, do.fast = T, dim_embed=2, perplexity= perp)

# png(filename=paste(project_name, res, Dim, perp, "TSNEPlot.ident.singlet.cca.png", sep="."), width=600, height=600, bg = "white", res = 100)
# TSNEPlot(object =Seurat_object.integrated, do.return = T, no.legend = T, do.label = T, group.by = "ident")
# dev.off()


#Seurat_object <- RunUMAP(Seurat_object, reduction.use = "pca", dims.use = DIM_nb)


# png(filename=paste(project_name, res, Dim, perp, "UMAPlot.ident.singlet.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
# DimPlot(Seurat_object.integrated, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
#   cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
#   cols.use = NULL, group.by = "ident", pt.shape = NULL,
#   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#   do.label = TRUE, label.size = 4, no.legend = TRUE,
#   coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
#   plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
#   sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
# dev.off()




saveRDS(Seurat_object.integrated, paste(project_name, res, Dim, perp, "Seurat_object.rds", sep="."))



q("no")
