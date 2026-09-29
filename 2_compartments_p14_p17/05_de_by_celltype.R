# ---------------------------------------------------------------------------
# Step 05 - Differential expression between cell types.
#
# ESTABLISHES
#     Which genes and pathways distinguish the annotated cell types from one
#     another, pooled across condition.
#
# WHY THIS WAY
#     Cell-type contrasts are run before condition contrasts (step 06) because
#     a gene that separates cell types will also appear in a condition contrast
#     wherever composition shifts; having the cell-type axis characterised
#     first is what makes the condition result interpretable.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 03
#
# INPUTS
#     <B_DIR>/CellTypeAnnotation/<project_name>.<res>.<Dim>.<perp>.Seurat_object.annotated.rds
#     <GMT_DIR> MSigDB v6.2 collections
#
# OUTPUTS
#     DE tables and figures under <B_DIR>/DifferentialAnalysis/CellType/
#
# USAGE
#     Rscript 05_*.R
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
library(limma)

### Downgrading Seurat from v3.1.0 to v2.3.4
##source("https://z.umn.edu/archived-seurat")


setwd(file.path(B_DIR, "CellTypeAnnotation/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))

summary(Seurat_object@meta.data)

summary(Seurat_object@ident)

summary(as.factor(Seurat_object@meta.data$cell_type))

batch_long <- Seurat_object@meta.data
batch_assigned <- batch_long 
batch_assigned <- batch_assigned %>%
   mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
batch_assigned <- as.data.frame(batch_assigned[,"Cond_TimePoint"])
 
rownames(batch_assigned) <- rownames(batch_long)
 
colnames(batch_assigned) <- "Cond_TimePoint"
 
Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned, col.name = "Cond_TimePoint")

###Re edit Seurat Object idents

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Quiescent Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Activated Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Mki67 Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Early Muller Glia"), new.ident.name = "Muller glia")

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Basal cells 1", "Basal cells 2", "Red blood cells"), accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


#Subset from P12 to P17

Subseted_cells_P12_P14_P17 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P12", "P14", "P17")))

Seurat_object.P12_P14_P17 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P12_P14_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


#Subset P14

Subseted_cells_P14 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14")))

Seurat_object.P14 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P14, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

#Subset P17

Subseted_cells_P17 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P17")))

Seurat_object.P17 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)



#Subset NORM

Subseted_cells_P17_NORM <- rownames(subset(Seurat_object.P17@meta.data, Condition %in% c("NORM")))

Subseted_cells_P17_NORM <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17_NORM, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

Subseted_cells_P14_NORM <- RenameIdent(Subseted_cells_P14_NORM, old.ident.name = c("Quiescent Muller Glia"), new.ident.name = "Muller glia")
Subseted_cells_P14_NORM <- RenameIdent(Subseted_cells_P14_NORM, old.ident.name = c("Activated Muller Glia"), new.ident.name = "Muller glia")
Subseted_cells_P14_NORM <- RenameIdent(Subseted_cells_P14_NORM, old.ident.name = c("Mki67 Muller Glia"), new.ident.name = "Muller glia")
Subseted_cells_P14_NORM <- RenameIdent(Subseted_cells_P14_NORM, old.ident.name = c("Early Muller Glia"), new.ident.name = "Muller glia")

Subseted_cells_P14_NORM <- SubsetData(Subseted_cells_P14_NORM, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Basal cells 1", "Basal cells 2", "Red blood cells"), accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

#Subset OIR

Subseted_cells_P17_OIR <- rownames(subset(Seurat_object.P17@meta.data, Condition %in% c("OIR")))

Subseted_object_P17_OIR <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17_OIR, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


####Subclustering
### Make distribution plot for P12 to P17

setwd(file.path(B_DIR, "DifferentialAnalysis/CellType/Pathways"))



png("Seurat_object_DAND5.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 150)
SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=c("DAND5"),
  cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png("Seurat_object_HMGCL_OXCT1.P14_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 150)
SplitDotPlotGG(Seurat_object.P14, grouping.var="Condition", genes.plot=c("HMGCL", "OXCT1"),
  cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png("Seurat_object_IL1Signalling.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 200)
SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=rev(c("IL1A", "IL1B", "IL1R1", "IL1RN", "IL1RAP", "TNF", "SEMA3A")),
  cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png("Seurat_object_IL1Signalling_alt.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 200)
SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=rev(c("IL6", "IL33")),
  cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

Seurat_object.P17_NORM <- ReorderIdent(Seurat_object.P17_NORM, feature = c("HMGCS2","HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1"), 
  rev = FALSE, aggregate.fxn = mean,
  reorder.numeric = FALSE)

png("Seurat_object_HMGCL_HMGCS2.P17_NORM.png", width=1000, height=1000, bg = "white", res = 200)
DotPlot(Seurat_object.P17_NORM, genes.plot=rev(c("HMGCS2","HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1")),
  cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5, scale.by = "radius",
  dot.min = 0, dot.scale = 8, group.by= "Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

cluster_counts <- Seurat_object.P12_P14_P17@meta.data %>%
  group_by(Cond_TimePoint) %>%
  count(cell_type) %>%
  mutate(freq = n / sum(n))
  
cluster_counts$cell_type <- factor(cluster_counts$cell_type,levels=unique(mixedsort(cluster_counts$cell_type)))

png("Cluster_proportions_CellType.P12toP17.Cond_TimePoint.png", width=2000, height=1500, bg = "white", res = 100)
ggplot(cluster_counts,aes(Cond_TimePoint,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ cell_type, scale="free") +
  theme_light()
dev.off()


ggsave(cluster_portions,file="Cluster_proportions_CellType.P12toP17.Cond_TimePoint.png", width = 6, height = 3)

cluster_portions_full <-  ggplot(cluster_counts,aes(cell_type,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.P12toP17.Cond_TimePoint.png")


####Subclustering cell types

Subseted_cells_ECs <- rownames(subset(Subseted_cells_P17_NORM@meta.data, cell_type %in% c("Endothelial cells")))

Seurat_object.P17_ECs <- SubsetData(Subseted_cells_P17_NORM, cells.use = Subseted_cells_ECs, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)




png("Seurat_object_IL1Signalling.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 200)
SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=rev(c("IL1A", "IL1B", "IL1R1", "IL1RN", "IL1RAP", "TNF", "SEMA3A")),
  cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()


#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17_ECs, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png("Seurat_object_HMGCL.P12_P17_OIRvsNORM.Cond_TimePoint.ECs.png", width=1000, height=1000, bg = "white", res = 200)
#~ DotPlot(Seurat_object.P12_P14_P17_ECs, genes.plot=rev(c("HMGCL")),
#~   cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 8, group.by= "Cond_TimePoint", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()


#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17_ECs, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png("Seurat_object_HMGCL.P12_P17_OIRvsNORM.Cond_TimePoint.ECs.png", width=1000, height=1000, bg = "white", res = 200)
#~ DotPlot(Seurat_object.P12_P14_P17_ECs, genes.plot=rev(c("HMGCL")),
#~   cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 8, group.by= "Cond_TimePoint", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()

Subseted_cells_ImmuneCells <- rownames(subset(Seurat_object.P12_P14_P17@meta.data, cell_type %in% c("Immune cells")))

Seurat_object.P12_P14_P17_ImmuneCells <- SubsetData(Seurat_object.P12_P14_P17, cells.use = Subseted_cells_ImmuneCells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

png("Seurat_object_IL1Signalling.P12_P17_OIRvsNORM.Cond_TimePoint.ImmuneCells.png", width=1000, height=1000, bg = "white", res = 200)
DotPlot(Seurat_object.P12_P14_P17_ImmuneCells, genes.plot=rev(c("IL1A", "IL1B")),
  cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by= "Cond_TimePoint", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()


### Analysis of ontological functions

GO_score=data.frame(fread(file.path(B_DIR, "GSVA/NotImputed/Subset1000CellperIdent/gsva.exprs.Full_eset.NonImputed.scaled.sup_h.c2.c5.c6.c7.270519.csv"), sep=",", header=TRUE), row.names=1)

#GO_score <- read.table(file.path(B_DIR, "GSVA/NotImputed/retina_p17_data_gsva_out_senescence_May2020.csv"), sep = ",", header = T, row.names=1, stringsAsFactors=F)

Seurat_object_subset <- SubsetData(Seurat_object.P17, cells.use = colnames(GO_score), subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

GO_score <- GO_score[,colnames(Seurat_object_subset@data)]

Seurat_object_subset <- SetAssayData(Seurat_object_subset, assay.type = "GO", slot = "raw.data", new.data = GO_score)

Seurat_object_subset <- NormalizeData(Seurat_object_subset, assay.type = "GO", normalization.method = "genesCLR")

Seurat_object_subset <- ScaleData(Seurat_object_subset, assay.type = "GO", display.progress = FALSE, model.use = "negbinom")

###Run Dim reduction on subset data


Seurat_object_subset <- RunTSNE(Seurat_object_subset,
                               reduction.use = "cca.aligned",
                               dims.use = DIM_nb, do.fast = T, dim_embed=2, perplexity=perp)

setwd(file.path(B_DIR, "DifferentialAnalysis/CellType/Mapping"))

png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "ident.cca.P17.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
TSNEPlot(object = Seurat_object_subset, do.return = T, no.legend = T, do.label = T)
dev.off()
#
png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "Cond_Sorting.cca.png", sep="."), width=1500, height=1000, bg = "white", res = 150)
TSNEPlot(object =Seurat_object_subset, do.return = T, no.legend = F, do.label = F, group.by = "Cond_Sorting")
dev.off()

Seurat_object_subset <- RunUMAP(Seurat_object_subset, reduction.use = "cca.aligned", dims.use = DIM_nb,
  genes.use = NULL, assay.use = "RNA", max.dim = 2L,
  reduction.name = "umap", reduction.key = "UMAP", n_neighbors = 50L,
  min_dist = 0.7, metric = "correlation", seed.use = 42)


png(filename=paste("UmapPlot", project_name, res, Dim, perp, "ident.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object_subset, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()

png(filename=paste("UmapPlot", project_name, res, Dim, perp, "Cond_Sorting.cca.png", sep="."), width=1500, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object_subset, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "Cond_Sorting", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()



##Find the max expressed functions

setwd(file.path(B_DIR, "DifferentialAnalysis/CellType/Functions/Global"))

test <- Seurat_object_subset

test@data <- GetAssayData(Seurat_object_subset, assay.type = "GO", slot = "data")


test.average <- AverageExpression(test, genes.use = NULL, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE,
  show.progress = TRUE)

dataframe <- head(test.average[order(test.average$Endothelial, decreasing= T),], n = 50)

my_palette <- colorRampPalette(c("blue", "white", "red"))

dataframe <- as.matrix(dataframe)


png(filename=paste("DoHeatmap_top50_pathway_ECs.png", sep="."), width=2000, height=2000, res = 150)  
heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
  Rowv=TRUE,
  Colv=TRUE,
  srtCol=45,
  revC=TRUE, key.title = NA, density.info="none", lhei=c(1,8), lwid=c(1,6), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,35), key.ylab=NA)
dev.off()




##For conditional regulation on a cell type basis in NORM or OIR

setwd(file.path(B_DIR, "DifferentialAnalysis/CellType/Functions/OIR"))

condition <- "OIR"

Subseted_cells_P17_condition <- rownames(subset(Seurat_object_subset@meta.data, Condition %in% c(condition)))

Seurat_object_subset_P17_condition <- SubsetData(Seurat_object_subset, cells.use = Subseted_cells_P17_condition, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)


test <- Seurat_object_subset_P17_condition

test@data <- GetAssayData(Seurat_object_subset_P17_condition, assay.type = "GO", slot = "data")


test.average <- AverageExpression(test, genes.use = NULL, return.seurat = FALSE,
  add.ident = NULL, use.scale = FALSE, use.raw = FALSE,
  show.progress = TRUE)

dataframe <- head(test.average[order(test.average$Endothelial, decreasing= T),], n = 50)

my_palette <- colorRampPalette(c("blue", "white", "red"))

dataframe <- as.matrix(dataframe)


png(filename=paste("DoHeatmap_top50_pathway_ECs.P17",condition,"png", sep="."), width=2000, height=2000, res = 150)  
heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
  Rowv=TRUE,
  Colv=TRUE,
  srtCol=45,
  revC=TRUE, key.title = NA, density.info="none", lhei=c(1,8), lwid=c(1,6), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,35), key.ylab=NA)
dev.off()



png(filename=paste("TSNEPlot", project_name, res, Dim, perp, "Seurat_object_subset_P17",condition,"ident.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
TSNEPlot(object = Seurat_object_subset_P17_condition, do.return = T, no.legend = T, do.label = T)
dev.off()

svg(paste("TSNEPlot", project_name, res, Dim, perp, "Seurat_object_subset_P17",condition,"ident.cca.svg", sep="."), width=7, height=7, bg = "white")
TSNEPlot(object = Seurat_object_subset_P17_condition, do.return = T, no.legend = T, do.label = T)
dev.off()



png(filename=paste("UmapPlot", project_name, res, Dim, perp, "Seurat_object_subset_P17",condition,"ident.cca.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
DimPlot(Seurat_object_subset, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()


svg(filename=paste("UmapPlot", project_name, res, Dim, perp, "Seurat_object_subset_P17",condition,"ident.svg.png", sep="."), width=7, height=7, bg = "white")
DimPlot(Seurat_object_subset, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = TRUE, label.size = 4, no.legend = TRUE,
  coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
  plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
  sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
dev.off()





grep.genes <- grep("SENESCENCE", rownames(GO_score))


grep.genes <- c(grep(c("MIR191"), rownames(GO_score)),
  grep(c("MIR126"), rownames(GO_score)),
  grep(c("MIR143"), rownames(GO_score)),
  grep(c("MIR340"), rownames(GO_score)),
  grep(c("MIR30B"), rownames(GO_score)),
  grep(c("MIR335"), rownames(GO_score)),
  grep(c("MIR30E"), rownames(GO_score)),
  grep(c("MIR382"), rownames(GO_score)),
  grep(c("LET7C"), rownames(GO_score)),
  grep(c("MIR211"), rownames(GO_score)),
  grep(c("LET7A"), rownames(GO_score)),
  grep(c("LET7E"), rownames(GO_score)),
  grep(c("MIR9B"), rownames(GO_score)),
  grep(c("MIR181"), rownames(GO_score)),
  grep(c("MIR96"), rownames(GO_score)),
  grep(c("LET7B"), rownames(GO_score)),
  grep(c("MIR30C"), rownames(GO_score)),
  grep(c("MIR135"), rownames(GO_score)),
  grep(c("MIR190"), rownames(GO_score)),
  grep(c("MIR30E"), rownames(GO_score)),
  grep(c("MIR30A"), rownames(GO_score)),
  grep(c("MIR125"), rownames(GO_score))
  )

grep.genes <- c(grep("BETA_OXIDATION", rownames(GO_score)),
                grep(c("FATTY_ACID_OXIDATION"), rownames(GO_score)),
                grep(c("RESPIRATION"), rownames(GO_score))
               )

rownames(GO_score[unique(grep.genes),])

#pathway <- "WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1"

pathway <- rownames(GO_score[unique(grep.genes),])


for (p in pathway) {
    Seurat_object_subset_P17_condition <- ReorderIdent(Seurat_object_subset_P17_condition, feature = p, rev = FALSE, aggregate.fxn = mean,
    reorder.numeric = FALSE)
    png(filename=paste(p,"RidgePlot_Markers.P17",condition,"scaled.new.png", sep="."), width = 800, height = 500, res=150)
    RidgePlot <- RidgePlot(Seurat_object_subset_P17_condition, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
    dev.off()
    svg(filename=paste(p,"RidgePlot_Markers.P17",condition,"scaled.new.svg", sep="."), width = 8, height = 5)
    RidgePlot <- RidgePlot(Seurat_object_subset_P17_condition, p, do.return = TRUE, size.title.use = 8)
    print(RidgePlot)
    dev.off()
  }


color <- colorRampPalette(rev(c("red", "yellow", "green")))(n = 10)

#color <- colorRampPalette(viridis(100))(n = 10)


for (p in pathway) {
    png(filename=paste(p,"FeaturePlot_tsne_Markers.P17",condition,"scaled.new.png", sep="."), width = 800, height = 800, res=150)
    FeaturePlot <- FeaturePlot(Seurat_object_subset_P17_condition, p, min.cutoff = NA, max.cutoff = NA,
    dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
    cols.use = color, pch.use = 16, overlay = FALSE,
    do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
    reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
    no.axes = FALSE, no.legend = FALSE)
    print(FeaturePlot)
    dev.off()
    svg(filename=paste(p,"FeaturePlot_tsne_Markers.P17",condition,"scaled.new.svg", sep="."), width = 7, height = 7)
    FeaturePlot <- FeaturePlot(Seurat_object_subset_P17_condition, p, min.cutoff = NA, max.cutoff = NA,
    dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
    cols.use = color, pch.use = 16, overlay = FALSE,
    do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
    reduction.use = "tsne", use.imputed = FALSE, nCol = NULL,
    no.axes = FALSE, no.legend = FALSE)
    print(FeaturePlot)
    dev.off()
  }


for (p in pathway) {
    png(filename=paste(p,"FeaturePlot_umap_Markers.P17_OIR",condition,"scaled.new.png", sep="."), width=2000, height=2000, bg = "white", res = 300)
    FeaturePlot <- FeaturePlot(Seurat_object_subset_P17_condition, p, min.cutoff = NA, max.cutoff = NA,
    dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
    cols.use = color, pch.use = 16, overlay = FALSE,
    do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
    reduction.use = "umap", use.imputed = FALSE, nCol = NULL,
    no.axes = FALSE, no.legend = FALSE)
    print(FeaturePlot)
    dev.off()
    #
    svg(filename=paste(p,"FeaturePlot_umap_Markers.P17_OIR",condition,"scaled.new.svg", sep="."), width = 7, height = 7)
    FeaturePlot <- FeaturePlot(Seurat_object_subset_P17_condition, p, min.cutoff = NA, max.cutoff = NA,
    dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
    cols.use = color, pch.use = 16, overlay = FALSE,
    do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
    reduction.use = "umap", use.imputed = FALSE, nCol = NULL,
    no.axes = FALSE, no.legend = FALSE)
    print(FeaturePlot)
    dev.off()
    #
    png(filename=paste("UmapPlot", project_name, res, Dim, perp, "ident.cca.P17_OIR.png", sep="."), width=2000, height=2000, bg = "white", res = 300)
    DimPlot(Seurat_object_subset_P17_condition, reduction.use = "umap", dim.1 = 1, dim.2 = 2,
    cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
    cols.use = NULL, group.by = "ident", pt.shape = NULL,
    do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
    do.label = TRUE, label.size = 4, no.legend = TRUE,
    coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE,
    plot.order = NULL, cells.highlight = NULL, cols.highlight = "red",
    sizes.highlight = 1, plot.title = NULL, vector.friendly = FALSE)
    dev.off()

  }

### Calculating frequency of positive cell for specific GO term


Seurat_object_gated_pos_cells <- WhichCells(Seurat_object_subset_P17_condition, ident = NULL, ident.remove = NULL, cells.use = NULL,
  subset.name = "FRIDMAN_SENESCENCE_UP", accept.low = 0, accept.high = Inf,
  accept.value = NULL, max.cells.per.ident = Inf, random.seed = 1)


Seurat_object_subset_P17_condition <- SetIdent(Seurat_object_subset_P17_condition, cells.use = Seurat_object_gated_pos_cells, ident.use = "FRIDMAN_SENESCENCE_UP_Positive")

Seurat_object_gated_neg_cells <- WhichCells(Seurat_object_subset_P17_condition, ident = NULL, ident.remove = NULL, cells.use = NULL,
  subset.name = "FRIDMAN_SENESCENCE_UP", accept.low = -Inf, accept.high = 0,
  accept.value = NULL, max.cells.per.ident = Inf, random.seed = 1)

Seurat_object_subset_P17_condition <- SetIdent(Seurat_object_subset_P17_condition, cells.use = Seurat_object_gated_neg_cells, ident.use = "FRIDMAN_SENESCENCE_UP_Negative")


Seurat_object_subset_P17_condition <- StashIdent(Seurat_object_subset_P17_condition, save.name = "FRIDMAN_SENESCENCE_UP")

table(Seurat_object_subset_P17_condition@meta.data$FRIDMAN_SENESCENCE_UP)



cluster_counts <- Seurat_object_subset_P17_condition@meta.data %>%
  group_by(cell_type) %>%
  count(FRIDMAN_SENESCENCE_UP) %>%
  mutate(freq = n / sum(n))

write.csv(cluster_counts, "cluster_counts.FRIDMAN_SENESCENCE_UP.csv")

cluster_counts <- subset(cluster_counts, FRIDMAN_SENESCENCE_UP != "FRIDMAN_SENESCENCE_UP_Negative")


#cluster_counts <- rbind(as.data.frame(cluster_counts),c("Rods","SENESCENCE_SAPIEHA_Positive",0,0))

cluster_counts$cell_type <- factor(cluster_counts$cell_type,levels=cluster_counts$cell_type[order(cluster_counts$freq, decreasing=FALSE)])


#cluster_counts$cell_type <- factor(cluster_counts$cell_type,levels=rev(c("Pericytes", "Endothelial cells", "Astrocytes", "Muller glia","Horizontal cells", 
#																	"Immune cells", "Retinal ganglion cells", "Bipolar cells", "Amacrine cells", "Cones", "Rods")))

#cluster_counts$freq <- as.numeric(cluster_counts$freq)

png("Cluster_proportions_full_CellType.OIR_P17.FRIDMAN_SENESCENCE_UP.png", width=2000, height=1000, bg = "white", res = 300)
ggplot(cluster_counts,aes(cell_type,freq)) +
  geom_bar(stat="identity", colour="black", fill="darkmagenta")+
 coord_flip()
dev.off()





#END
q("no")
