# ---------------------------------------------------------------------------
# 01_immune_cells - 03_de_by_condition.R
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# immune subtypes, run on cells taken from the annotated
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
library(magrittr)
library(harmony)
library(RColorBrewer)
library(limma)

#Set directory of dataset to analyse


setwd(file.path(B_DIR, "Subclustering/ImmuneCells/CellIdentification_5"))

microglia <- readRDS("Microglia_subclustered_renamed.rds")

summary(microglia@ident)
summary(microglia@meta.data)
summary(microglia@meta.data$Cond_TimePoint)


# for differential biological function between condition

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/Conditions/Functions"))

Seurat_object <- microglia 

GO_score=data.frame(fread(file.path(B_DIR, "Subclustering/ImmuneCells/GSVA/NotImputed/gsva.exprs.Full_eset.NonImputed.scaled.poisson.sup_h.c2.c5.c6.c7.csv"), sep=",", header=TRUE), row.names=1)

GO_score_t <- as.data.frame(t(GO_score))

colnames(GO_score_t)[grep("INNATE_IMMUNE", colnames(GO_score_t))]

png("Histogram_GO_INNATE_IMMUNE_RESPONSE.scaled.png")
ggplot(GO_score_t, aes(x=GO_score_t$GO_INNATE_IMMUNE_RESPONSE)) + geom_histogram(binwidth=.005)
dev.off()


expr.csv =data.frame(fread(file.path(B_DIR, "Subclustering/ImmuneCells/GSVA/NotImputed/expr.csv"), sep=",", header=TRUE), row.names=1)

expr.csv_t <- as.data.frame(t(expr.csv))

colnames(expr.csv_t)[grep("CX", colnames(expr.csv_t))]

png("Histogram_ACTB.png")
ggplot(expr.csv_t, aes(x=expr.csv_t$ACTB)) + geom_histogram(binwidth=.005)
dev.off()

GO_score <- GO_score[,colnames(Seurat_object@data)]

Seurat_object <- SetAssayData(Seurat_object, assay.type = "non_imputed_GO", slot = "raw.data", new.data = GO_score)

Seurat_object <- NormalizeData(Seurat_object, assay.type = "non_imputed_GO", normalization.method = "genesCLR")

Seurat_object <- ScaleData(Seurat_object, assay.type = "non_imputed_GO", display.progress = FALSE, model.use = "negbinom")


Seurat_object_scaled_score <- GetAssayData(Seurat_object, assay.type = "non_imputed_GO", slot = "scale.data")

Seurat_object_scaled_score_t <- as.data.frame(t(Seurat_object_scaled_score))

png("Histogram_GO_INNATE_IMMUNE_RESPONSE.afterNorm_Scaled.png")
ggplot(Seurat_object_scaled_score_t, aes(x=Seurat_object_scaled_score_t$GO_INNATE_IMMUNE_RESPONSE)) + geom_histogram(binwidth=.005)
dev.off()



###Function of interest

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/Functions/Selection")

#i <- "Cones"

grep.genes <- grep("PHOTO", rownames(GO_score))

rownames(GO_score[grep.genes,])

pathway <- "WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1"

#pathway <- rownames(GO_score)

for (i in unique(Seurat_object@ident)) {
  Seurat_object_subset <- SubsetData(Seurat_object, ident.use =i)
  Seurat_object_subset <- SetAllIdent(Seurat_object_subset, id = "Dataset")
  for (p in pathway) {
    Seurat_object_subset <- ReorderIdent(Seurat_object_subset, feature = p, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
    RidgePlot <- RidgePlot(Seurat_object_subset, p, do.return = TRUE, , size.title.use = 10)
    ggsave(RidgePlot, filename=paste(i, p, "RidgePlot_Markers.png", sep="."), width = 15, height = 7, dpi=150, units = "cm")
    VlnPlot <- VlnPlot(Seurat_object_subset, p, ident.include = NULL, nCol = NULL,
    do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
    size.y.use = 16, size.title.use = 10, adjust.use = 1,
    point.size.use = -1, cols.use = NULL, group.by = NULL, y.log = FALSE,
    x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
    single.legend = TRUE, remove.legend = FALSE, do.return = TRUE)
    ggsave(VlnPlot, filename=paste(i, p, "VlnPlot_Markers.png", sep="."), width = 15, height = 15, dpi=150, units = "cm")
  }
}

Seurat_object.GO.markers <- FindAllMarkers(object =Seurat_object, only.pos = FALSE, assay.type = "non_imputed_GO", test.use = "wilcox")

Seurat_object.GO.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file =paste("CellType_GO.Markers.c5-top10.tsv", sep="."), sep = "\t")


png(filename=paste("HeatMap_GO_marker.png", sep="."), width=3000, height=2500, bg = "white", res = 150)
DoHeatmap(Seurat_object, genes.use = unique(top10$gene), assay.type = "non_imputed_GO", 
    slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE)
dev.off()


#i <- "Cones"

#p <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"

for (i in unique(Seurat_object@ident)) {
  Seurat_object_subset <- SubsetData(Seurat_object, ident.use =i)
  Seurat_object_subset <- SetAllIdent(Seurat_object_subset, id = "orig.ident")
  for (p in rownames(GO_score)) {
    RidgePlot <- RidgePlot(Seurat_object_subset, p, do.return = TRUE, , size.title.use = 10)
    ggsave(RidgePlot, filename=paste(project_name, res, Dim, perp, i, p, "RidgePlot_Markers.png", sep="."), width = 15, height = 7, dpi=150, units = "cm")
    VlnPlot <- VlnPlot(Seurat_object_subset, p, ident.include = NULL, nCol = NULL,
    do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
    size.y.use = 16, size.title.use = 10, adjust.use = 1,
    point.size.use = -1, cols.use = NULL, group.by = NULL, y.log = FALSE,
    x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
    single.legend = TRUE, remove.legend = FALSE, do.return = TRUE)
    ggsave(VlnPlot, filename=paste(project_name, res, Dim, perp, i, p, "VlnPlot_Markers.png", sep="."), width = 15, height = 15, dpi=150, units = "cm")
  }
}



for (i in unique(top10$gene)){
  Seurat_object <- ReorderIdent(Seurat_object, feature = i, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
	png(filename=paste(i, "FeaturePlot.png", sep="."), height=1000, width=1000, res=150)
	FeaturePlot(Seurat_object, i, min.cutoff = NA, max.cutoff = NA,
  	dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1.5,
  	cols.use = c("blue", "green", "yellow", "orange", "red"), pch.use = 16, overlay = FALSE,
  	do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  	reduction.use = "umap_harmony", use.imputed = FALSE, nCol = NULL,
  	no.axes = FALSE, no.legend = TRUE, coord.fixed = FALSE,
  	dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
	dev.off()
  RidgePlot <- RidgePlot(Seurat_object, i, do.return = TRUE, , size.title.use = 10)
  ggsave(RidgePlot, filename=paste(i,"RidgePlot_Markers.png", sep="."), width = 15, height = 7, dpi=150, units = "cm")
}

for (i in rownames(GO_score)){

	png(filename=paste(project_name, res, Dim, perp, i, "FeaturePlot.png", sep="."), height=1000, width=1000, res=150)
	FeaturePlot(Seurat_object, i, min.cutoff = NA, max.cutoff = NA,
  	dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1.5,
  	cols.use = c("blue", "purple", "green", "yellow", "orange", "red"), pch.use = 16, overlay = FALSE,
  	do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  	reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
  	no.axes = FALSE, no.legend = TRUE, coord.fixed = FALSE,
  	dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
	dev.off()

}


#### For differential expressed functions between conditions globally

##Analysis global

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/Functions/DEFs/Global")

#cell_type_list <- sort(unique(Seurat_object@ident))

cell_type_list <- c("Pericytes")                 

cell_type_list <- unique(Seurat_object@ident)

for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Dataset")
  
  #Perform DGE analysis using one of the model above
  DGE_test <- FindAllMarkers(cells_in_this_cluster, genes.use = NULL,
          logfc.threshold = 0, test.use = "wilcox", min.pct = 0.001,
          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
          min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
          assay.type = "non_imputed_GO")

  ## Write table for all differentially expressed genes containing testing results
  write.table(DGE_test,file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t",
                                     quote=FALSE,
                                     row.names=TRUE,
                                     col.names=TRUE)
  DGE_test <- read.table(file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t")

  DGE_test %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
  color <- colorRampPalette(brewer.pal(11,"Spectral"))(3)
  RidgePlot <- RidgePlot(cells_in_this_cluster, top10$gene, do.sort = TRUE, do.return = TRUE, , size.title.use = 7, cols.use =rev(color))
  ggsave(RidgePlot, filename=paste("cluster_",this_cluster,"RidgePlot_Top10_functions_Markers.png", sep="."), width = 50, height = 50, dpi=100, units = "cm")
}

Venn_Cluster1 <- subset(DGE_test, cluster == "OIR.P17.CD31.WT.S129")

Venn_Cluster2 <- subset(DGE_test, cluster == "NORM.P17.CD31.WT.S129")

Venn_Cluster3 <- subset(DGE_test, cluster == "OIR.P17.CD31.Sirt3KO.S129")

my_list <- list(Venn_Cluster1$gene, Venn_Cluster2$gene, Venn_Cluster3$gene)

library(VennDiagram)

venn.diagram(x = list(OIR.WT=Venn_Cluster1$gene, NORM.WT=Venn_Cluster2$gene, OIR.Sirt3KO=Venn_Cluster3$gene),
        filename = paste("cluster_",this_cluster,"_VennDiagram_significant_Markers_genes.Wilcox.png")
        )

}



###Optional: Subset cells:

Subseted_cells <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14", "P17")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = FALSE, do.scale = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE)

#For all geneset

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/Functions/DEFs/OIRvsNORM")

#cell_type_list <- "Pericytes"

cell_type_list <- unique(Seurat_object_Subset@ident)

for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Dataset")
  
  #Perform DGE analysis using one of the model above
  DGE_test <- FindAllMarkers(cells_in_this_cluster, genes.use = NULL,
          logfc.threshold = 0, test.use = "wilcox", min.pct = 0.001,
          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
          min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
          assay.type = "non_imputed_GO")

  ## Write table for all differentially expressed genes containing testing results
  write.table(DGE_test,file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t",
                                     quote=FALSE,
                                     row.names=TRUE,
                                     col.names=TRUE)
  DGE_test <- read.table(file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t")

  DGE_test %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
  top10 <- as.data.frame(top10)
  color <- colorRampPalette(brewer.pal(11,"Spectral"))(2)
  #RidgePlot <- RidgePlot(cells_in_this_cluster, top10$gene)
  #ggsave(RidgePlot, filename=paste("cluster_",this_cluster,"RidgePlot_Top10_functions_Markers.png", sep="."), width = 50, height = 50, dpi=100, units = "cm")
  cells_in_this_cluster_Rep <- SetAllIdent(cells_in_this_cluster, id = "Batch")
  cells_in_this_cluster_Rep@data <- GetAssayData(cells_in_this_cluster_Rep, assay.type = "non_imputed_GO", slot = "scale.data")
  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster_Rep, genes.use = top10$gene, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")
  my_palette <- colorRampPalette(c("blue", "white", "red"))
  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)
  png(filename=paste("cluster_",this_cluster,"DoHeatmap_Top10_functions_Markers.png", sep="."), width=1000, height=1000, res = 100)  
  heatmap.2(average_cells_in_this_cluster, col=my_palette, symbreak=FALSE, trace='none', cexRow=1, cexCol= 1,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(1,10), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,25), key.ylab=NA, srtCol=45)
  dev.off()
}

#For GO all geneset

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/Functions/GO/OIRvsNORM")

cell_type_list <- unique(Seurat_object_Subset@ident)

GO.genes <- grep("GO_", rownames(GO_score))
GO.genes <- rownames(GO_score[GO.genes,])

for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Dataset")
  
  #Perform DGE analysis using one of the model above
  DGE_test <- FindAllMarkers(cells_in_this_cluster, genes.use = GO.genes,
          logfc.threshold = 0, test.use = "wilcox", min.pct = 0.001,
          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
          min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
          assay.type = "non_imputed_GO")

  ## Write table for all differentially expressed genes containing testing results
  write.table(DGE_test,file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t",
                                     quote=FALSE,
                                     row.names=TRUE,
                                     col.names=TRUE)
  DGE_test <- read.table(file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t")
  DGE_test <- 
  

  DGE_test %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
  top10 <- as.data.frame(top10)
  color <- colorRampPalette(brewer.pal(11,"Spectral"))(2)
  #RidgePlot <- RidgePlot(cells_in_this_cluster, top10$gene)
  #ggsave(RidgePlot, filename=paste("cluster_",this_cluster,"RidgePlot_Top10_functions_Markers.png", sep="."), width = 50, height = 50, dpi=100, units = "cm")
  cells_in_this_cluster_Rep <- SetAllIdent(cells_in_this_cluster, id = "Batch")
  cells_in_this_cluster_Rep@data <- GetAssayData(cells_in_this_cluster_Rep, assay.type = "non_imputed_GO", slot = "data")
  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster_Rep, genes.use = top10$gene, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")
  my_palette <- colorRampPalette(c("blue", "white", "red"))
  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)
  png(filename=paste("cluster_",this_cluster,"DoHeatmap_Top10_GO_functions_Markers.png", sep="."), width=1000, height=1000, res = 100)  
  heatmap.2(average_cells_in_this_cluster, col=my_palette, symbreak=FALSE, trace='none', cexRow=1, cexCol= 1,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(1,10), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,25), key.ylab=NA, srtCol=45)
  dev.off()
}

#For Selected geneset

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/Conditions/Functions/Selection"))

cell_type_list <- unique(Seurat_object_Subset@ident)

Selected.genes_1 <- grep(c("FATTY"), rownames(GO_score))

Selected.genes_2 <- grep(c("MITOCHONDRIA"), rownames(GO_score))

Selected.genes_1 <- rownames(GO_score[Selected.genes_1,])

Selected.genes_2 <- rownames(GO_score[Selected.genes_2,])

Selected.genes <- unique(c(Selected.genes_1,Selected.genes_2))

for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Dataset")
  
  #Perform DGE analysis using one of the model above
  DGE_test <- FindAllMarkers(cells_in_this_cluster, genes.use = Selected.genes,
          logfc.threshold = 0, test.use = "wilcox", min.pct = 0.001,
          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
          min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
          assay.type = "non_imputed_GO")

  ## Write table for all differentially expressed genes containing testing results
  write.table(DGE_test,file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t",
                                     quote=FALSE,
                                     row.names=TRUE,
                                     col.names=TRUE)
  DGE_test <- read.table(file=paste("cluster_",this_cluster,"_significant_Markers_function.Wilcox.txt"), sep="\t")
  DGE_test %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
  top10 <- as.data.frame(top10)
  color <- colorRampPalette(brewer.pal(11,"Spectral"))(2)
  #RidgePlot <- RidgePlot(cells_in_this_cluster, top10$gene)
  #ggsave(RidgePlot, filename=paste("cluster_",this_cluster,"RidgePlot_Top10_functions_Markers.png", sep="."), width = 50, height = 50, dpi=100, units = "cm")
  cells_in_this_cluster_Rep <- SetAllIdent(cells_in_this_cluster, id = "Batch")
  cells_in_this_cluster_Rep@data <- GetAssayData(cells_in_this_cluster_Rep, assay.type = "non_imputed_GO", slot = "data")
  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster_Rep, genes.use = unique(top10$gene), return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")
  my_palette <- colorRampPalette(c("blue", "white", "red"))
  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)
  png(filename=paste("cluster_",this_cluster,"DoHeatmap_FATTY_ACID_MITOCHONDRIA_functions_Markers.png", sep="."), width=1000, height=1000, res = 100)  
  heatmap.2(average_cells_in_this_cluster, col=my_palette, symbreak=FALSE, trace='none', cexRow=1, cexCol= 1,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(1,10), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,25), key.ylab=NA, srtCol=45)
  dev.off()
}


for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)
  cells_in_this_cluster_Rep <- SetAllIdent(cells_in_this_cluster, id = "Batch")
  cells_in_this_cluster_Rep@data <- GetAssayData(cells_in_this_cluster_Rep, assay.type = "non_imputed_GO", slot = "data")
  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster_Rep, genes.use = Selected.genes, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")
  my_palette <- colorRampPalette(c("blue", "white", "red"))
  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)
  png(filename=paste("cluster_",this_cluster,"DoHeatmap_FATTY_ACID_MITOCHONDRIA_functions_Markers.png", sep="."), width=1500, height=2000, res = 100)  
  heatmap.2(average_cells_in_this_cluster, col=my_palette, symbreak=FALSE, trace='none', cexRow=1, cexCol= 1,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(1,10), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,45), key.ylab=NA, srtCol=45)
  dev.off()
}

Selected.genes <- c("NEUTROPHIL_BINET")


for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)
  png(filename=paste("cluster",this_cluster,Selected.genes,"RidgePlot_Markers.png", sep="."), width = 700, height = 500, res=150)
  p <- RidgePlot(cells_in_this_cluster, Selected.genes, ident.include = NULL, nCol = NULL,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 8,
  size.y.use = 8, size.title.use = 8, cols.use = NULL,
  group.by = "Cell_type", y.log = FALSE, x.lab.rot = FALSE, y.lab.rot = FALSE,
  legend.position = "right", single.legend = TRUE, remove.legend = FALSE,
  do.return = FALSE, return.plotlist = FALSE)
  print(p)
  dev.off()
}


## For all cell type 

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/Conditions/Functions/Selection"))

Selected.genes <- c(grep(c("NEUTROPHIL"), rownames(GO_score))
  #,grep(c("GRANULO"), rownames(GO_score))
  )

Selected.genes <- rownames(GO_score[Selected.genes,])

Selected.genes <- unique(c(Selected.genes))

unwanted_genes <- Selected.genes[grep("GSE", Selected.genes)]

Selected.genes <- Selected.genes[!Selected.genes %in% unwanted_genes]

#cell_type_list <- "Pericytes"

cell_type_list <- unique(Seurat_object_Subset@ident)

dataframe <- data.frame()

dataframe <- dataframe[1:length(Selected.genes),]

rownames(dataframe) <- Selected.genes

for(this_cluster in cell_type_list){

  #####this_cluster <- "Cluster 1"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Condition")
  
  #Perform DGE analysis using one of the model above
  cells_in_this_cluster@data <- GetAssayData(cells_in_this_cluster, assay.type = "non_imputed_GO", slot = "scale.data")


  average_cells_in_this_cluster_seurat <- AverageExpression(cells_in_this_cluster, genes.use = NULL, return.seurat = TRUE,
  add.ident = NULL, use.scale = TRUE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")


  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster, genes.use = NULL, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")


  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)


  design <- model.matrix(~0+average_cells_in_this_cluster_seurat@ident)

  colnames(design) <- c("NORM","OIR")

  contrast_dir <- paste(colnames(design)[1], "vs", colnames(design)[2], sep="")

  contrast <- makeContrasts(OIR - NORM, levels = design)

  fit <- lmFit(average_cells_in_this_cluster, design)

  fit <- contrasts.fit(fit, contrast)

  sel_FC <- as.data.frame(fit$coefficients[Selected.genes,])

  colnames(sel_FC) <- paste(this_cluster)

  dataframe <- cbind(dataframe,sel_FC)
}

my_palette <- colorRampPalette(c("blue", "white", "red"))

dataframe <- as.matrix(dataframe)

#dataframe2 <- dataframe[-9,]

png(filename=paste("DoHeatmap_Selection_NEUTROPHILS_functions_Markers_updated2.png", sep="."), width=1500, height=1000, res = 120)  
heatmap.2(dataframe2, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, Rowv=TRUE,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(3,12), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(20,45), key.ylab=NA, srtCol=45)
dev.off()

svg(filename=paste("DoHeatmap_Selection_NEUTROPHILS_functions_Markers_updated2.svg", sep="."), width=12, height=8)  
heatmap.2(dataframe2, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, Rowv=TRUE,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(3,12), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(20,45), key.ylab=NA, srtCol=45)
dev.off()


###FOR OIR WT vs SIRT3KO
#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/Functions/DEFs/OIRSirt3KOvsWT")

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("OIR.P17.CD31.Sirt3KO.S129","OIR.P17.CD31.WT.S129")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = FALSE, do.scale = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE)

Selected.genes <- c(grep(c("HALLMARK_GLYCOLYSIS"), rownames(GO_score)),
  grep(c("HALLMARK_FATTY_ACID_METABOLISM"), rownames(GO_score)),
  grep(c("REACTOME_FATTY_ACID"), rownames(GO_score)),
  grep(c("GO_REGULATION_OF_FATTY_ACID_TRANSPORT"), rownames(GO_score)),
  grep(c("GO_POSITIVE_REGULATION_OF_FATTY_ACID_BIOSYNTHETIC"), rownames(GO_score)),
  grep(c("MOOTHA_GLYCOLYSIS"), rownames(GO_score))
  )

Selected.genes <- rownames(GO_score[Selected.genes,])

Selected.genes <- unique(c(Selected.genes))

cell_type_list <- "Pericytes"

cell_type_list <- unique(Seurat_object_Subset@ident)

dataframe <- data.frame()

dataframe <- dataframe[1:length(Selected.genes),]

rownames(dataframe) <- Selected.genes

for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Dataset")
  
  #Perform DGE analysis using one of the model above
  cells_in_this_cluster@data <- GetAssayData(cells_in_this_cluster, assay.type = "non_imputed_GO", slot = "scale.data")
  average_cells_in_this_cluster_seurat <- AverageExpression(cells_in_this_cluster, genes.use = NULL, return.seurat = TRUE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")
  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster, genes.use = NULL, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "non_imputed_GO")
  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)


  design <- model.matrix(~0+average_cells_in_this_cluster_seurat@ident)

  colnames(design) <- c("OIR.P17.CD31.Sirt3KO.S129","OIR.P17.CD31.WT.S129")

  contrast_dir <- paste(colnames(design)[1], "vs", colnames(design)[2], sep="")

  contrast <- makeContrasts(OIR.P17.CD31.Sirt3KO.S129 - OIR.P17.CD31.WT.S129, levels = design)

  fit <- lmFit(average_cells_in_this_cluster, design)

  fit <- contrasts.fit(fit, contrast)

  sel_FC <- as.data.frame(fit$coefficients[Selected.genes,])

  colnames(sel_FC) <- paste(this_cluster)

  dataframe <- cbind(dataframe,sel_FC)
}

my_palette <- colorRampPalette(c("blue", "white", "red"))

dataframe <- as.matrix(dataframe)

png(filename=paste("DoHeatmap_Selection_GLYCOLYSIS_FATTY_ACID_functions_Markers.png", sep="."), width=1500, height=1000, res = 120)  
heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, Rowv=TRUE,
  Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(3,10), lwid=c(2,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(20,45), key.ylab=NA, srtCol=45)
dev.off()


Selected.genes <- c("HALLMARK_FATTY_ACID_METABOLISM")


for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)
  png(filename=paste("cluster",this_cluster,Selected.genes,"RidgePlot_Markers.png", sep="."), width = 700, height = 500, res=150)
  p <- RidgePlot(cells_in_this_cluster, Selected.genes, ident.include = NULL, nCol = NULL,
  do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 8,
  size.y.use = 8, size.title.use = 8, cols.use = NULL,
  group.by = "Dataset", y.log = FALSE, x.lab.rot = FALSE, y.lab.rot = FALSE,
  legend.position = "right", single.legend = TRUE, remove.legend = FALSE,
  do.return = FALSE, return.plotlist = FALSE)
  print(p)
  dev.off()
}


#### For DGE between conditions globally

##Analysis global

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/DEGs/Global")

#cell_type_list <- sort(unique(Seurat_object@ident))

cell_type_list <- c("Pericytes")                 

#cell_type_list <- unique(Seurat_object@ident)

for(this_cluster in cell_type_list){

  #####this_cluster <- "Amacrine cells 1"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
        cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Dataset")
  
        #Perform DGE analysis using one of the model above
        DGE_test <- FindAllMarkers(cells_in_this_cluster, genes.use = NULL,
          logfc.threshold = 1, test.use = "wilcox", min.pct = 0.1,
          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
          min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
          assay.type = "RNA")

    ## Write table for all differentially expressed genes containing testing results
        write.table(DGE_test,file=paste("cluster_",this_cluster,"_significant_Markers_genes.Wilcox.txt"), sep="\t",
                                     quote=FALSE,
                                     row.names=TRUE,
                                     col.names=TRUE)


        DGE_test %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

        dotplot <- DotPlot(cells_in_this_cluster, unique(top10$gene), cols.use = c("blue", "red"),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 7,
        scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
        plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE)

        ggsave(dotplot, filename=paste("cluster_",this_cluster,"_top10_significant_DE_genes.Wilcox.png"), dpi=150, width = 15, height=5)


Venn_Cluster1 <- subset(DGE_test, cluster == "OIR.P17.CD31.WT.S129")

Venn_Cluster2 <- subset(DGE_test, cluster == "NORM.P17.CD31.WT.S129")

Venn_Cluster3 <- subset(DGE_test, cluster == "OIR.P17.CD31.Sirt3KO.S129")

my_list <- list(Venn_Cluster1$gene, Venn_Cluster2$gene, Venn_Cluster3$gene)

library(VennDiagram)

venn.diagram(x = list(OIR.WT=Venn_Cluster1$gene, NORM.WT=Venn_Cluster2$gene, OIR.Sirt3KO=Venn_Cluster3$gene),
        filename = paste("cluster_",this_cluster,"_VennDiagram_significant_Markers_genes.Wilcox.png")
        )

}



q("no")
