# ---------------------------------------------------------------------------
# 01_immune_cells - 02_de_by_subtype.R
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
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(GSVA) 
data(c2BroadSets)
library(gplots)
library(heatmap3)
library(devtools)
library(fgsea)
library(tidyverse)
library(R.utils)






#Set directory of dataset to analyse


setwd(file.path(B_DIR, "Subclustering/ImmuneCells/CellIdentification_5"))

microglia <- readRDS("Microglia_subclustered_renamed.rds")

summary(microglia@ident)
table(microglia@meta.data$Cond_TimePoint)

# for gene expression between cell types 

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/CellType/PathwaysGenes"))

png(filename="FeaturePlot_microglia_subclustered_IL1R1.png", width=700, height=500, bg = "white", res = 150)
FeaturePlot(microglia, c("IL1R1"), min.cutoff = NA, max.cutoff = NA,
  dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
  cols.use = c("aquamarine2", "red"), pch.use = 16, overlay = FALSE,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  reduction.use = "umap_cca", use.imputed = FALSE, nCol = NULL,
  no.axes = FALSE, no.legend = FALSE, 
  dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
dev.off()

png(filename="DotPlot_microglia_subclustered_CX3CR1_Cond_TimePoint.png", width=700, height=800, bg = "white", res = 150)
DotPlot(microglia, "CX3CR1", cols.use = c("blue", "red"),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 7,
        scale.by = "radius", scale.min = NA, scale.max = NA, group.by="Cond_TimePoint",
        plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE)
dev.off()        


markers <- FindMarkers(microglia, ident.1="Angioblast", ident.2 = NULL, genes.use = NULL,
  logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1,
  min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
  min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
  assay.type = "RNA")


head(markers)

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on VascularEndothelium, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png(filename="FeaturePlot_VascularEndothelium_subclustered_HBA-A1_HBA-A2.png", width=800, height=500, bg = "white", res = 150)
#~ FeaturePlot(VascularEndothelium, c("HBA-A1", "HBA-A2"), min.cutoff = NA, max.cutoff = NA,
#~   dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
#~   cols.use = c("aquamarine2", "orange", "plum", "red"), pch.use = 16, overlay = TRUE,
#~   do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
#~   reduction.use = "umap_harmony", use.imputed = FALSE, nCol = NULL,
#~   no.axes = FALSE, no.legend = FALSE, 
#~   dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
#~ dev.off()

# for differential pathways between cell types 

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/CellType/Functions"))

Seurat_object <- microglia 

GO_score=data.frame(fread(file.path(B_DIR, "Subclustering/ImmuneCells/GSVA/NotImputed/gsva.exprs.Full_eset.NonImputed.scaled.poisson.sup_h.c2.c5.c6.c7.csv"), sep=",", header=TRUE), row.names=1)

GO_score_t <- as.data.frame(t(GO_score))

colnames(GO_score_t)[grep("BINET", colnames(GO_score_t))]

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


## Differential expressed function between cell types
setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/CellType/Functions/DEFs"))

Seurat_object.GO.markers <- FindAllMarkers(object =Seurat_object, only.pos = FALSE, assay.type = "non_imputed_GO", test.use = "wilcox")

Seurat_object.GO.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file =paste("CellType_GO.Markers.c5-top10.tsv", sep="."), sep = "\t")


png(filename=paste("HeatMap_GO_marker.png", sep="."), width=3000, height=2500, bg = "white", res = 150)
DoHeatmap(Seurat_object, genes.use = unique(top10$gene), assay.type = "non_imputed_GO", 
    slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE)
dev.off()


#i <- "Cones"

p <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"

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


###Function of interest accross cell types

##In NORM

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/CellType/Functions/Selection/NORM"))

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Condition %in% c("NORM")))

Seurat_object_subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)


Seurat_object.GO.markers <- FindAllMarkers(object =Seurat_object_subset, only.pos = FALSE, assay.type = "non_imputed_GO", test.use = "wilcox")

Seurat_object.GO.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file =paste("CellType_GO.Markers-top10.tsv", sep="."), sep = "\t")


png(filename=paste("HeatMap_GO_marker.png", sep="."), width=3000, height=1500, bg = "white", res = 150)
DoHeatmap(Seurat_object_subset, genes.use = unique(top10$gene), assay.type = "non_imputed_GO", 
    slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE)
dev.off()



grep.genes <- grep("GRANULATION", rownames(GO_score))

rownames(GO_score[grep.genes,])

#pathway <- "WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1"

pathway <- rownames(GO_score[grep.genes,])

unwanted_genes <- pathway[grep("GSE", pathway)]

pathway <- pathway[!pathway %in% unwanted_genes]

pathway <- unique(top10$gene)

for (p in pathway) {
    Seurat_object_subset <- ReorderIdent(Seurat_object_subset, feature = p, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
    png(filename=paste(p,"RidgePlot_Markers.png", sep="."), width = 800, height = 500, res=150)
    RidgePlot <- RidgePlot(Seurat_object_subset, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
    dev.off()
  }


##In OIR

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/CellType/Functions/Selection/OIR"))

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Condition %in% c("OIR")))

Seurat_object_subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

Seurat_object.GO.markers <- FindAllMarkers(object =Seurat_object_subset, only.pos = FALSE, assay.type = "non_imputed_GO", test.use = "wilcox")

Seurat_object.GO.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

write.table(top10, file =paste("CellType_GO.Markers-top10.tsv", sep="."), sep = "\t")


png(filename=paste("HeatMap_GO_marker.png", sep="."), width=3000, height=1500, bg = "white", res = 150)
DoHeatmap(Seurat_object_subset, genes.use = unique(top10$gene), assay.type = "non_imputed_GO", 
    slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE)
dev.off()



grep.genes <- c(grep("GRANULATION", rownames(GO_score)), grep("GRANULOCYTE", rownames(GO_score)))

rownames(GO_score[grep.genes,])

#pathway <- "WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1"

pathway <- rownames(GO_score[grep.genes,])

unwanted_genes <- pathway[grep("GSE", pathway)]

pathway <- pathway[!pathway %in% unwanted_genes]

pathway <- unique(top10$gene)

for (p in pathway) {
    Seurat_object_subset <- ReorderIdent(Seurat_object_subset, feature = p, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
    png(filename=paste(p,"RidgePlot_Markers.png", sep="."), width = 800, height = 500, res=150)
    RidgePlot <- RidgePlot(Seurat_object_subset, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
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

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/Aligned/Subclustering/VascularEndothelium/DifferentialAnalysis/Conditions/Functions/DEFs/OIRvsNORM")

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("NORM.P17.CD31.WT.S129", "OIR.P17.CD31.WT.S129")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = FALSE, do.scale = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE)

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
  cells_in_this_cluster_Rep@data <- GetAssayData(cells_in_this_cluster_Rep, assay.type = "non_imputed_GO", slot = "data")
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

#Only NORM

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Rytvela/Rytvela_OIRvsNORM_P14_2run/Clustering/Not_Aligned/DifferentialAnalysis/Conditions/Functions/Alex_Dubrac/NORM")

Subseted_cells <- rownames(subset(Seurat_object@meta.data, orig.ident %in% c("NORM.P14")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = FALSE, do.scale = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE)


for (p in rownames(GO_score)) {
    RidgePlot <- RidgePlot(Seurat_object_Subset, p, do.return = TRUE, size.title.use = 10)
    ggsave(RidgePlot, filename=paste(project_name, res, Dim, perp, p, "RidgePlot_Markers.NORM.P14.png", sep="."), width = 15, height = 15, dpi=150, units = "cm")
    VlnPlot <- VlnPlot(Seurat_object_Subset, p, ident.include = NULL, nCol = NULL,
    do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
    size.y.use = 16, size.title.use = 10, adjust.use = 1,
    point.size.use = -1, cols.use = NULL, group.by = NULL, y.log = FALSE,
    x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
    single.legend = TRUE, remove.legend = FALSE, do.return = TRUE)
    ggsave(VlnPlot, filename=paste(project_name, res, Dim, perp, p, "VlnPlot_Markers.NORM.P14.png", sep="."), width = 20, height = 15, dpi=150, units = "cm")
}



###FOR OIR
#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Rytvela/Rytvela_OIRvsNORM_P14_2run/Clustering/Not_Aligned/DifferentialAnalysis/Conditions/Functions/Alex_Dubrac/OIR")

Subseted_cells <- rownames(subset(Seurat_object@meta.data, orig.ident %in% c("OIR.P14")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = FALSE, do.scale = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE)

for (p in rownames(GO_score)) {
    RidgePlot <- RidgePlot(Seurat_object_Subset, p, do.return = TRUE, size.title.use = 10)
    ggsave(RidgePlot, filename=paste(project_name, res, Dim, perp, p, "RidgePlot_Markers.OIR.P14.png", sep="."), width = 15, height = 15, dpi=150, units = "cm")
    VlnPlot <- VlnPlot(Seurat_object_Subset, p, ident.include = NULL, nCol = NULL,
    do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
    size.y.use = 16, size.title.use = 10, adjust.use = 1,
    point.size.use = -1, cols.use = NULL, group.by = NULL, y.log = FALSE,
    x.lab.rot = TRUE, y.lab.rot = FALSE, legend.position = "right",
    single.legend = TRUE, remove.legend = FALSE, do.return = TRUE)
    ggsave(VlnPlot, filename=paste(project_name, res, Dim, perp, p, "VlnPlot_Markers.OIR.P14.png", sep="."), width = 20, height = 15, dpi=150, units = "cm")
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


####Gating on  enrichment score for specific function

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/DifferentialAnalysis/CellType/Functions/Selection/Global"))

Seurat_object_gated_cells <- WhichCells(Seurat_object, ident = "Cluster 3", ident.remove = NULL, cells.use = NULL,
  subset.name = "NEUTROPHIL_BINETupdated2020", accept.low = 0, accept.high = Inf,
  accept.value = NULL, max.cells.per.ident = Inf, random.seed = 1)

Seurat_object <- StashIdent(Seurat_object, save.name = "ClutserIdent")

Seurat_object <- SetIdent(Seurat_object, cells.use = Seurat_object_gated_cells, ident.use = "CLuster_3_NEUTROPHIL_BINET_POSITIVE")

Seurat_object <- StashIdent(Seurat_object, save.name = "Cell_type")

Subseted_cells <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14", "P17")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells.use = Subseted_cells, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = TRUE, do.scale = TRUE,
  max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE)



grep.genes <- grep("NEUTROPHIL_BINET", rownames(GO_score))

rownames(GO_score[grep.genes,])

#pathway <- "WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1"

pathway <- rownames(GO_score[grep.genes,])

unwanted_genes <- pathway[grep("GSE", pathway)]

pathway <- pathway[!pathway %in% unwanted_genes]


for (p in pathway) {
    Seurat_object_Subset <- SetAllIdent(Seurat_object_Subset, id="Cell_type")
    Seurat_object_Subset <- ReorderIdent(Seurat_object_Subset, feature = p, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
    png(filename=paste(p,"RidgePlot_Markers.poisson.05.new.png", sep="."), width = 1300, height = 500, res=150)
    RidgePlot <- RidgePlot(Seurat_object_Subset, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
    dev.off()
    #print in svg
    svg(paste(p,"RidgePlot_Markers.05.new.poisson.svg", sep="."), width = 10, height = 5)
    RidgePlot <- RidgePlot(Seurat_object_Subset, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
    dev.off()
  }

for (p in pathway) {
    Seurat_object_Subset <- SetAllIdent(Seurat_object_Subset, id="ClutserIdent")
    Seurat_object_Subset <- ReorderIdent(Seurat_object_Subset, feature = p, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
    svg(paste(p,"RidgePlot_Markers.05.new.ClutserIdent.poisson.svg", sep="."), width = 10, height = 5)
    RidgePlot <- RidgePlot(Seurat_object_Subset, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
    dev.off()
  }



cluster_counts <- Seurat_object_Subset@meta.data %>%
  group_by(ClutserIdent) %>%
  dplyr::count(Cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$Cell_type <- factor(cluster_counts$Cell_type,levels=unique(mixedsort(cluster_counts$Cell_type)))

cluster_portions <-  ggplot(cluster_counts,aes(ClutserIdent,freq,fill=ClutserIdent)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_type, scale="free") +
  theme_light()

ggsave(cluster_portions,file="Cluster_proportions_CellType.microglia.P14_P17.ClutserIdent.05new.png", width = 12, height = 5)

write.table(cluster_counts, "cluster_counts.CellType.microglia.P14_P17.ClutserIdent.txt")

cluster_portions_full <-  ggplot(cluster_counts,aes(Cell_type,freq,fill=ClutserIdent)) +
  geom_bar(stat="identity", col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.microglia.P14_P17.ClutserIdent.05new.png")


Seurat_object.GO.markers <- FindMarkers(Seurat_object_Subset, ident.1="CLuster_3_NEUTROPHIL_BINET_POSITIVE", ident.2 = NULL, genes.use = NULL,
  logfc.threshold = 0, test.use = "wilcox", min.pct = 0,
  min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
  min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
  assay.type = "RNA")


top20 <- head(Seurat_object.GO.markers, 20) 

png("DotPlot_Top20_NEUTROPHIL_BINET_POSITIVE.05.new.png", width = 1500, height = 700, res=150)
DotPlot(Seurat_object_Subset, unique(rownames(top20)), cols.use = c("blue", "red"),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 7,
        scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
        plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE)
dev.off()

Seurat_object.GO.Allmarkers <- FindAllMarkers(Seurat_object_Subset)

Seurat_object.GO.Allmarkers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20

top20.neutro <- subset(top20, cluster=="CLuster_3_NEUTROPHIL_BINET_POSITIVE")


png("DotPlot_Top20_markers_NEUTROPHIL_BINET_POSITIVE.05.new.png", width = 1500, height = 700, res=150)
DotPlot(Seurat_object_Subset, unique(top20.neutro$gene), cols.use = c("blue", "red"),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 7,
        scale.by = "radius", scale.min = NA, scale.max = NA, group.by="ident",
        plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE)
dev.off()



geneset <- read.csv(file.path(GENE_LIST_DIR, "Neutrophil_Binet_gene_set.csv"), header=FALSE)
GOI <- intersect(rownames(Seurat_object_Subset@data),toupper(geneset$V1))


broadset_selected <- "sup_h.c2.c5.c6.c7.all.v6.2.symbols"

broadSet.custom <- getGmt(paste(file.path(GMT_DIR, ""), broadset_selected, ".gmt", sep=""),
              geneIdType=SymbolIdentifier())

genesetrows <- grep("BINETupdated2020", names(broadSet.custom))

genesetnames <- names(broadSet.custom[genesetrows])

geneIDs <- geneIds(broadSet.custom[[genesetnames]])

GO_term <- intersect(rownames(Seurat_object_Subset@data),GOI)

GO_term_bis <- c("IL1RN", "CXCL1", "CLEC4D", "CXCR2", "LY6G6D", "SOD2")

png("DotPlot_NEUTROPHIL_BINET_POSITIVE_Neutro_Binet_genes.05.new_updated2020_short.png", width = 1500, height = 700, res=150)
DotPlot(Seurat_object_Subset, unique(c(GO_term,GO_term_bis)), cols.use = c("blue", "red"),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 7,
        scale.by = "radius", scale.min = NA, scale.max = NA, group.by="Cell_type",
        plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE)
dev.off()

svg("DotPlot_NEUTROPHIL_BINET_POSITIVE_Neutro_Binet_genes.05.new_updated2020_short.svg", width = 12, height = 4)
DotPlot(Seurat_object_Subset, unique(c(GO_term,GO_term_bis)), cols.use = c("blue", "red"),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 8,
        scale.by = "radius", scale.min = NA, scale.max = NA, group.by="Cell_type",
        plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE)
dev.off()









###Fgsea

GOI <- genesetnames

celltypes <- names(summary(Seurat_object_Subset@ident))

#celltypes <- "Horizontal cells"

###Start the loop for fgsea on each cell type

##Subset NORM vs OIR


Conditions <- c("NORM", "OIR")

for (Condition_subset in Conditions) {

  Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Condition %in% c(Condition_subset)))

  retina_subset <- SubsetData(Seurat_object_Subset, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = 50,
  random.seed = 1)


  df_total = NULL
  for (celltypename in celltypes) {
    retina_CTL_celltype <- SubsetData(retina_subset, cells.use = NULL, subset.name =NULL, ident.use = celltypename,
      ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
      do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
      random.seed = 1)

    print(summary(retina_CTL_celltype@ident))

    retina_CTL_celltype.scaled <- as.data.frame(retina_CTL_celltype@scale.data)

    CellBarCode_name <- unique(colnames(retina_CTL_celltype.scaled))
    
    gsea=NULL
    
    for (sel_CellBarCode in CellBarCode_name) {

    #CellBarCode_cell <- rownames(subset(retina_CTL_celltype@meta.data, CellBarCode == sel_CellBarCode))

    #retina_CellBarCode_CTL <- SubsetData(retina_CTL_celltype, cells.use = CellBarCode_cell, subset.name =NULL, ident.use = NULL,
    #ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
    #do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
    #random.seed = 1)

    #Plot enrichment graph for the cell type (average)

      #if (length(rownames(retina_CTL_celltype@meta.data)) >1) {

      #average_cell_Type <- AverageExpression(retina_CTL_celltype, genes.use = NULL, return.seurat = FALSE,
      #add.ident = NULL, use.scale = TRUE, use.raw = FALSE,
      #show.progress = TRUE)

      #sel_CellBarCode <- "NORM_P17_WR_Joyal_r1_CTAAGCTAAAAA"

      average_cell_Type <- retina_CTL_celltype.scaled %>% dplyr::select(paste(sel_CellBarCode))

      gene_names <- rownames(average_cell_Type)

      average_cell_Type <- unlist(average_cell_Type, recursive = TRUE, use.names = TRUE)

      names(average_cell_Type) <- gene_names

    
      #plot <- plotEnrichment(GO_term[[1]], average_cell_Type) + labs(title=paste(celltypename,GOI,"_EnrichmentPlot.scaled.png", sep="_"))
    
      #ggsave(plot, file=paste(celltypename,"EnrichmentPlot.scaled.png", sep="_"))


      fgseaRes <- fgsea(pathways = GO_term, stats = average_cell_Type,
                    minSize=1,
                    maxSize=500,
                    nperm=1000,
                    nproc=16,
                    gseaParam=1)
      #}
      #else {
      #fgseaRes <- fgsea(pathways = GO_term, stats = retina_CTL_celltype.scaled,
      #              minSize=15,
      #              maxSize=500,
      #              nperm=1000,
      #              nproc=16)
      #}
      fgseaRes <-data.frame(fgseaRes, cell_type=celltypename, CellBarCode=sel_CellBarCode)
      gsea <- rbind(gsea,fgseaRes)
    }
  df_total <- rbind(df_total,gsea)
  }


  fwrite(df_total, file =paste("GSEA_CellBarCode", GOI, Condition_subset,"scaled.txt", sep="_"), sep="\t")

  #df_total <- fread("GSEA_Batch_GOautophagy.txt", header=T)

  #Make nice bar plot


  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                        conf.interval=.95, .drop=TRUE) {
      library(plyr)

      # New version of length which can handle NA's: if na.rm==T, don't count them
      length2 <- function (x, na.rm=FALSE) {
          if (na.rm) sum(!is.na(x))
          else       length(x)
      }

      # This does the summary. For each group's data frame, return a vector with
      # N, mean, and sd
      datac <- ddply(data, groupvars, .drop=.drop,
        .fun = function(xx, col) {
          c(N    = length2(xx[[col]], na.rm=na.rm),
            mean = mean   (xx[[col]], na.rm=na.rm),
            sd   = sd     (xx[[col]], na.rm=na.rm)
          )
        },
        measurevar
      )

      # Rename the "mean" column    
      datac <- rename(datac, c("mean" = measurevar))

      datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

      # Confidence interval multiplier for standard error
      # Calculate t-statistic for confidence interval: 
      # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
      ciMult <- qt(conf.interval/2 + .5, datac$N-1)
      datac$ci <- datac$se * ciMult

      return(datac)
  }


  tgc <- summarySE(df_total, measurevar="NES", groupvars="cell_type")

  tgc_padj <- summarySE(df_total, measurevar="padj", groupvars="cell_type")

  tgc2 <- merge(tgc,tgc_padj,by="cell_type")

  tgc2$cell_type <- factor(tgc2$cell_type)


  png(paste("BarPlotNES.Retina", GOI, Condition_subset, "SEM_black.scaled.png", sep="_"), width=3000, height=2000, res=300)
  p <- ggplot(tgc2, aes(x=reorder(cell_type,-NES), y=NES, fill=NULL)) + 
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      xlab("Cell Type") +
      ylab("Mean of Normalized Enrichement Score (± SEM)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      scale_fill_distiller(palette = "YlOrRd", direction=-1, trans = "log10", breaks = c(0, 0.01, 0.05, 0.2, 0.8, 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      coord_flip()+ 
      theme_classic() +
      theme(legend.position="bottom") 
  print(p)
  dev.off()

  png(paste("BarPlotNES.Retina", GOI, Condition_subset, "SEM.scaled.colors.png", sep="_"), width=2000, height=2000, res=300)
  p <- ggplot(tgc2, aes(reorder(cell_type, NES),NES)) +
      geom_bar(aes(fill = padj), stat="identity") +
      geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9),
                    color="brown")+
      xlab("Cell Type") +
      ylab("Mean of Normalized Enrichement Score (± SEM)")+
      scale_fill_distiller(palette = "Reds", direction=-1, trans = "log10") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      coord_flip()+ 
      theme_classic() +
      theme(legend.position="bottom") 
  print(p)
  dev.off()

  tgc2$Condition <- Condition_subset

  tgc2$cell_type <- paste(Condition_subset, tgc2$cell_type, sep="_")

  assign(paste("tgc2", Condition_subset, sep="."),tgc2)


  #rbPal <- colorRampPalette(c('red','blue'))

  #tgc2.OIR$Colour <- rbPal(100)[as.numeric(cut(tgc2.OIR$padj,breaks = 100))]



  write.table(tgc2, paste("tgc2", GOI, Condition_subset, "celltype_NES_summary.scaled.txt", sep="_"), sep="\t")


  aov.dfb <- aov(NES ~ cell_type, data=df_total)

  posthoc <- TukeyHSD(x=aov.dfb, "cell_type",  conf.level=0.95)


  summary(aov.dfb)

  print(model.tables(aov.dfb,"means"),digits=3)


  posthoc.save <- as.data.frame(posthoc$cell_type)

  write.table(posthoc.save, paste("Anova_Tukey_Retina", GOI, Condition_subset, "cellbarcode.scaled.txt", sep="_"), sep="\t")

}

###Mergin NORM and OIR

tgc.merged <- rbind(tgc2.NORM,tgc2.OIR)


tgc2.NORM <- tgc2.NORM[order(tgc2.NORM$NES,decreasing=TRUE),]

tgc2.NORM_cell <- as.character(tgc2.NORM$cell_type)
tgc2.NORM_NES <- tgc2.NORM$NES

#
tgc2.OIR_cell <- as.character(tgc2.OIR$cell_type)
tgc2.OIR_NES <- tgc2.OIR$NES

#
order_cell <- insert(tgc2.NORM_cell, "",ats = 2:length(tgc2.NORM_cell))
order_n <- insert(tgc2.NORM_NES, "",ats = 2:length(tgc2.NORM_NES))
order_cell[length(order_cell)+1]<-""
order_n[length(order_n)+1]<-""

#
i <- 1
while(1){
  tmp <- tgc2.OIR_cell[grep(sub("_","",str_match(order_cell[i],"_.*")),tgc2.OIR_cell)]
  order_cell[i+1] <- tmp
  order_n[i+1] <- tgc2.OIR$NES[which(tgc2.OIR$cell_type==tmp)]
  i <- i + 2
  if(i>length(order_cell)) break
}

#
ordered.df <- as.data.frame(cbind(order_cell,order_n))

ordered.df[,2] <- factor(ordered.df[,2], levels = unique(ordered.df[,2]))

tgc.merged$cell_type <- factor(tgc.merged$cell_type, level=rev(unique(ordered.df[,1])))

tgc.merged$NES <- as.numeric(tgc.merged$NES)


png(paste("BarPlotNES.Retina_", GOI, "_NORMvsOIR_SEM.scaled.colors.new.png", sep=""), width=2000, height=2000, res=300)
ggplot(tgc.merged, aes(cell_type, NES, fill=Condition)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9),
                  color="darkgrey")+
    xlab("Cell Type P17") +
    ylab("Mean of Normalized Enrichement Score (± SEM)")+ 
    #geom_col(fill = tgc.merged$Colour)+
    scale_fill_manual(values=c("limegreen", "firebrick1"))+
    #scale_fill_distiller(palette = "Reds", direction=-1, trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    coord_flip()+ 
    theme_classic() +
    theme(legend.position="bottom")+ 
    ggtitle(paste(GOI))
dev.off()



q("no")
