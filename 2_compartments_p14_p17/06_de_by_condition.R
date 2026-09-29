# ---------------------------------------------------------------------------
# Step 06 - Differential expression between OIR and normoxia.
#
# ESTABLISHES
#     The OIR-vs-normoxia contrast within each annotated cell type, at P14 and
#     at P17.
#
# WHY THIS WAY
#     Tested within cell type rather than across the whole retina, for the same
#     reason as the later analysis: OIR changes cell-type composition, and a
#     whole-retina contrast cannot separate a change in proportion from a
#     change in expression.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 03, and step 05 for the cell-type axis it is read against
#
# INPUTS
#     <B_DIR>/CellTypeAnnotation/<project_name>.<res>.<Dim>.<perp>.Seurat_object.annotated.rds
#     <GMT_DIR> MSigDB v6.2 collections
#
# OUTPUTS
#     DE tables and figures under <B_DIR>/DifferentialAnalysis/Condition/
#
# USAGE
#     Rscript 06_*.R
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
#library(plotly)
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
#library(Rmagic)
library(readr)
library(phateR)
library(viridis)
library(doubletFinder)
library(magrittr)
library(harmony)
library(limma)
library(dendextend)


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

summary(as.factor(Seurat_object@meta.data$Cond_Sorting))

batch_long <- Seurat_object@meta.data
batch_assigned <- batch_long 
batch_assigned <- batch_assigned %>%
   mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
batch_assigned <- as.data.frame(batch_assigned[,"Cond_TimePoint"])
 
rownames(batch_assigned) <- rownames(batch_long)
 
colnames(batch_assigned) <- "Cond_TimePoint"
 
Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned, col.name = "Cond_TimePoint")


###Re edit object ident

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Quiescent Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Activated Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Mki67 Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Early Muller Glia"), new.ident.name = "Muller glia")

##OPtional
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Muller glia"), new.ident.name = "Neuroglia cells")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Astrocytes"), new.ident.name = "Neuroglia cells")

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Basal cells 1", "Basal cells 2", "Red blood cells"), accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)



#Subset from P14 to P17

Subseted_cells_P14_P17 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14", "P17")))

Seurat_object.P14_P17 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P14_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


#Subset P14

Subseted_cells_P14 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14")))

Seurat_object.P14 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P14, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

#Subset P17

Subseted_cells_P17 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P17")))

Seurat_object.P17 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

Seurat_object.P17.endo_immune <- SubsetData(Seurat_object.P17, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = c("Endothelial cells", "Immune cells"),
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)


#Subset NORM

Subseted_cells_P17_NORM <- rownames(subset(Seurat_object.P17@meta.data, Condition %in% c("NORM")))

Seurat_object.P17_NORM <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17_NORM, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

#Subset OIR

Subseted_cells_P17_OIR <- rownames(subset(Seurat_object.P17@meta.data, Condition %in% c("OIR")))

Seurat_object.P17_OIR <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17_OIR, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

### Make distribution plot

setwd(file.path(B_DIR, "DifferentialAnalysis/Condition/Pathways"))

genesofinterest <- c("Ntn1", "Ntn3", "Ntn4", "Ntn5", "Unc5b", "Dcc", "Nrp1", "Nrp2", "Sema4a", "Sema4b", "Sema4c", "Sema4d", "Sema3a", "Sema3b", "Sema3c", "Sema3d", 
"Sema6a", "Sema6b", "Sema6c", "Sema6d", "Flrt1", "Flrt2", "Flrt3", "Plxna1", "Plxna2", "Plxna3", "Plxna4", "Plxnb1", "Plxnb2", "Plxnb3", "Plxnd1")

genelist <- read.csv(file.path(GENE_LIST_DIR, "SASP_Sawchyn_UP.csv"), header=FALSE)

genesofinterest <- intersect(rownames(Seurat_object.P17@data),genelist$V1)

genesofinterest_sel <- genesofinterest[c(-6, -16,-17, -31, -34)]

png("Seurat_object_SASP_Sawchyn_UP.P17_OIRvsNORM_full_new.png", width=3500, height=2000, bg = "white", res = 300)
SplitDotPlotGG(Seurat_object.P17, grouping.var="Condition", genes.plot=genesofinterest,
  cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

geneofinterest <- rownames(Seurat_object.P17@data[grep("IL6", rownames(Seurat_object.P17@data)),])

png("Seurat_object_PDE6.P14_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 150)
SplitDotPlotGG(Seurat_object.P14, grouping.var="Condition", genes.plot=rev(c(geneofinterest, "RGS9")),
  cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

png("Seurat_object_TGFB_pathway.scacolors.P17_OIRvsNORM.png", width=1000, height=1000, bg = "white", res = 150)
SplitDotPlotGG(Seurat_object.P17, grouping.var="Condition", genes.plot=c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3"),
  cols.use = c("deepskyblue4", "orangered"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png("Seurat_object_IL1Signalling.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 200)
#~ SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=rev(c("IL1A", "IL1B", "IL1R1", "IL1RN", "IL1RAP", "TNF", "SEMA3A")),
#~   cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png("Seurat_object_IL1Signalling_alt.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 200)
#~ SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=rev(c("IL6", "IL33")),
#~   cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()

Seurat_object.P14_P17 <- ReorderIdent(Seurat_object.P14_P17, feature = c("HMGCS2","HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1"), 
  rev = FALSE, aggregate.fxn = mean,
  reorder.numeric = FALSE)

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png("Seurat_object_HMGCL_HMGCS2.P14_P17_NORMvsOIR.png", width=1000, height=1000, bg = "white", res = 200)
#~ SplitDotPlotGG(Seurat_object.P14_P17, grouping.var="Condition", genes.plot=rev(c("HMGCS2","HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1")),
#~   cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 8, group.by= "Ident", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()
#~
#~ cluster_counts <- Seurat_object.P12_P14_P17@meta.data %>%
#~   group_by(Cond_TimePoint) %>%
#~   count(cell_type) %>%
#~   mutate(freq = n / sum(n))
#~
#~ cluster_counts$cell_type <- factor(cluster_counts$cell_type,levels=unique(mixedsort(cluster_counts$cell_type)))
#~
#~ png("Cluster_proportions_CellType.P12toP17.Cond_TimePoint.png", width=2000, height=1500, bg = "white", res = 100)
#~ ggplot(cluster_counts,aes(Cond_TimePoint,freq,fill=Cond_TimePoint)) +
#~   geom_bar(stat="identity",col="black") +
#~   facet_wrap(~ cell_type, scale="free") +
#~   theme_light()
#~ dev.off()


ggsave(cluster_portions,file="Cluster_proportions_CellType.P12toP17.Cond_TimePoint.png", width = 6, height = 3)

cluster_portions_full <-  ggplot(cluster_counts,aes(cell_type,freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  theme_light()

ggsave(cluster_portions_full,file="Cluster_proportions_full_CellType.P12toP17.Cond_TimePoint.png")


####Subclustering cell types

Subseted_cells_ECs <- rownames(subset(Seurat_object.P14_P17@meta.data, cell_type %in% c("Endothelial cells")))

Seurat_object.P14_P17_ECs <- SubsetData(Seurat_object.P14_P17, cells.use = Subseted_cells_ECs, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  random.seed = 1)

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ png("Seurat_object_IL1Signalling.P12_P17_OIRvsNORM.png", width=1200, height=1500, bg = "white", res = 200)
#~ SplitDotPlotGG(Seurat_object.P12_P14_P17, grouping.var="Condition", genes.plot=rev(c("IL1A", "IL1B", "IL1R1", "IL1RN", "IL1RAP", "TNF", "SEMA3A")),
#~   cols.use = c("blue", "darkgreen"), col.min = -2.5, col.max = 2.5,
#~   dot.min = 0, dot.scale = 8, group.by="Ident", plot.legend = TRUE,
#~   do.return = FALSE, x.lab.rot = TRUE)
#~ dev.off()


png("Seurat_object_HMGCL_SLC16A6_OXCT1.SLC16A7.BDH1.BDH2.P14_P17_OIRvsNORM.Cond_TimePoint.ECs.png", width=1500, height=1000, bg = "white", res = 200)
DotPlot(Seurat_object.P14_P17_ECs, genes.plot=rev(c("HMGCL","OXCT1", "SLC16A6",  "SLC16A7", "BDH1", "BDH2")),
  cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by= "Cond_TimePoint", plot.legend = TRUE,
  do.return = FALSE, x.lab.rot = TRUE)
dev.off()


png("Seurat_object_HMGCS2.SLC16A1.P14_P17_OIRvsNORM.Cond_TimePoint.ECs.png", width=800, height=1000, bg = "white", res = 200)
DotPlot(Seurat_object.P14_P17_ECs, genes.plot=rev(c("HMGCS2","SLC16A1")),
  cols.use = c("blue", "red"), col.min = -2.5, col.max = 2.5,
  dot.min = 0, dot.scale = 8, group.by= "Cond_TimePoint", plot.legend = TRUE,
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



Seurat_object.P14_P17_ECs <- SetAllIdent(Seurat_object.P14_P17_ECs, id = "Condition")

Markers <- FindMarkers(Seurat_object.P14_P17_ECs, ident.1="NORM", ident.2 = "OIR", 
  genes.use = c("HMGCL","OXCT1", "SLC16A6",  "SLC16A7", "BDH1", "BDH2","HMGCS2","SLC16A1"),
  logfc.threshold = -Inf, test.use = "wilcox", min.pct = -Inf,
  min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
  min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
  assay.type = "RNA")




#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ Subseted_cells_ImmuneCells <- rownames(subset(Seurat_object.P12_P14_P17@meta.data, cell_type %in% c("Immune cells")))

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This block operates on Seurat_object.P12_P14_P17, which is never loaded or assigned anywhere in
#~ this script. It was pasted from a session that had that object in memory,
#~ belonging to a different dataset. Disabled rather than deleted.
#~ ---------------------------------------------------------------------------
#~ Seurat_object.P12_P14_P17_ImmuneCells <- SubsetData(Seurat_object.P12_P14_P17, cells.use = Subseted_cells_ImmuneCells, subset.name =NULL, ident.use = NULL,
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

GO_score=data.frame(fread(file.path(B_DIR, "GSVA/NotImputed/Subset1000CellperIdent/gsva.exprs.Full_eset.NonImputed.scaled.c3.all.v6.2.symbols.csv"), sep=",", header=TRUE), row.names=1)

GO_score= read.table(file.path(B_DIR, "GSVA/NotImputed/retina_p17_data_gsva_out_senescence_feb2020.csv"), sep = ",", header = T, row.names=1, stringsAsFactors=F)

Seurat_object_subset <- SubsetData(Seurat_object.P17, cells.use = colnames(GO_score), subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

GO_score <- GO_score[,colnames(Seurat_object_subset@data)]

Seurat_object_subset <- SetAssayData(Seurat_object_subset, assay.type = "GO", slot = "raw.data", new.data = GO_score)

Seurat_object_subset <- NormalizeData(Seurat_object_subset, assay.type = "GO", normalization.method = "genesCLR")

Seurat_object_subset <- ScaleData(Seurat_object_subset, assay.type = "GO", display.progress = FALSE, model.use = "negbinom")


##For conditional regulation on a cell type basis between NORM and OIR from P12 to P17

setwd(file.path(B_DIR, "DifferentialAnalysis/Condition/Functions/P17"))

Selected.genes <- c(grep(c("KEGG_REGULATION_OF_AUTOPHAGY"), rownames(GO_score)),
  grep(c("GO_ENDOLYSOSOME_MEMBRANE"), rownames(GO_score)),
  grep(c("REACTOME_LYSOSOME_VESICLE_BIOGENESIS"), rownames(GO_score)),
  grep(c("GO_NEGATIVE_REGULATION_OF_MACROAUTOPHAGY"), rownames(GO_score)),
  grep(c("GO_AUTOPHAGY"), rownames(GO_score)),
  grep(c("NFAT_"), rownames(GO_score))
  )

Selected.genes <- c(grep(c("HALLMARK_GLYCOLYSIS"), rownames(GO_score)),
  grep(c("HALLMARK_FATTY_ACID_METABOLISM"), rownames(GO_score)),
  grep(c("REACTOME_FATTY_ACID"), rownames(GO_score)),
  grep(c("GO_REGULATION_OF_FATTY_ACID_TRANSPORT"), rownames(GO_score)),
  grep(c("GO_POSITIVE_REGULATION_OF_FATTY_ACID_BIOSYNTHETIC"), rownames(GO_score)),
  grep(c("MOOTHA_GLYCOLYSIS"), rownames(GO_score)),
  grep(c("GO_RESPONSE_TO_FOOD"), rownames(GO_score))
  )

Selected.genes <- c(grep(c("GLYCOLYSIS"), rownames(GO_score)),
  grep(c("FATTY_ACID"), rownames(GO_score)),
  grep(c("RESPIRATION"), rownames(GO_score)))


Selected.genes <- c(grep(c("KETO"), rownames(GO_score)))

Selected.genes <- c(grep(c("FRIDMAN_SENESCENCE_UP"), rownames(GO_score)),
  grep(c("GO_CELLULAR_SENESCENCE"), rownames(GO_score)),
  grep(c("SENESCENCE_MALLETTE"), rownames(GO_score)),
  grep(c("SENESCENCE_DEMARIA"), rownames(GO_score)),
  #grep(c("GO_REPLICATIVE_SENESCENCE"), rownames(GO_score)),
  grep(c("SENESCENCE_SASP"), rownames(GO_score)),
  grep(c("^GO_REGULATION_OF_CYTOKINE_SECRETION$"), rownames(GO_score))#,
  #grep(c("CYTOKINE_ACTIVITY"), rownames(GO_score)),
  #grep(c("^GO_REGULATION_OF_CYTOKINE_PRODUCTION$"), rownames(GO_score)),
  #grep(c("GO_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS"), rownames(GO_score))
  )

Selected.genes <- c(grep(c("BILD_HRAS_ONCOGENIC_SIGNATURE"), rownames(GO_score)),
  grep(c("SWEET_KRAS_ONCOGENIC_SIGNATURE"), rownames(GO_score)),
  grep(c("SWEET_LUNG_CANCER_KRAS_UP"), rownames(GO_score)),
  grep(c("HALLMARK_KRAS_SIGNALING_UP"), rownames(GO_score)),
  #grep(c("KRAS.LUNG_UP.V1_UP"), rownames(GO_score)),
  #grep(c("KRAS.PROSTATE_UP.V1_UP"), rownames(GO_score)),
  #grep(c("KRAS.300_UP.V1_UP"), rownames(GO_score)),
  grep(c("P53"), rownames(GO_score))
  )

Selected.genes <- c(grep(c("KRAS"), rownames(GO_score)),
  grep(c("HRAS"), rownames(GO_score))
  )

Selected.genes <- c(grep(c("MIR191"), rownames(GO_score)),
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


Selected.genes <- c("GLOBAL_SENESCENCE_SAPIEHA_UP", "FRIDMAN_SENESCENCE_UP",
  "MALLETTE_SENESCENCE_UP", "SENESCENCE_SASP_Sawchyn_UP", "GO_REGULATION_OF_CYTOKINE_SECRETION")#, "GO_CELLULAR_SENESCENCE"


##Plot heatmap

Selected.genes <- unique(c(Selected.genes))

Selected.genes <- rownames(GO_score[Selected.genes,])


#this_cluster <- "Cones"

cell_type_list <- unique(Seurat_object_subset@ident)

dataframe <- data.frame()

dataframe <- dataframe[1:length(Selected.genes),]

rownames(dataframe) <- Selected.genes

for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = "Condition")
  
  #Perform DGE analysis using one of the model above
  cells_in_this_cluster@scale.data <- GetAssayData(cells_in_this_cluster, assay.type = "GO", slot = "scale.data")
  #average_cells_in_this_cluster_seurat <- AverageExpression(cells_in_this_cluster, genes.use = NULL, return.seurat = FALSE,
  #add.ident = NULL, use.scale = FALSE, use.raw = FALSE, show.progress = TRUE, assay.type = "GO")
  average_cells_in_this_cluster <- AverageExpression(cells_in_this_cluster, genes.use = Selected.genes, return.seurat = FALSE,
  add.ident = NULL, use.scale = TRUE, use.raw = FALSE, show.progress = TRUE, assay.type = "GO")
  average_cells_in_this_cluster <- data.matrix(average_cells_in_this_cluster, rownames.force = TRUE)


  #design <- model.matrix(~0+average_cells_in_this_cluster_seurat@ident)

  design <- model.matrix(~0+colnames(average_cells_in_this_cluster))

  colnames(design) <- c("NORM","OIR")

  contrast_dir <- paste(colnames(design)[2], "vs", colnames(design)[1], sep="")

  contrast <- makeContrasts(OIR - NORM, levels = design)

  fit <- lmFit(average_cells_in_this_cluster, design)

  fit <- contrasts.fit(fit, contrast)

  sel_FC <- as.data.frame(fit$coefficients[Selected.genes,])

  colnames(sel_FC) <- paste(this_cluster)

  dataframe <- cbind(dataframe,sel_FC)
}

my_palette <- colorRampPalette(c("blue", "white", "red"))

dataframe <- as.matrix(dataframe)


Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
   set("branches_k_color", k = 3) %>% set("branches_lwd", 3) %>%
   ladderize %>% rotate_DendSer

#png("Dend.png")
#plot(Colv)
#dev.off()


Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
   set("branches_k_color", k = 5) %>% set("branches_lwd", 5) %>%
   ladderize


png(filename=paste("DoHeatmap_Selection_Senescence_Updated_noCellSEN_P17_OIRvsNORM.png", sep="."), width=2000, height=1000, res = 150)  
heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
  Rowv=Rowv,
  #Rowv=TRUE,
  Colv=Colv,
  #Colv=TRUE,
  srtCol=45,
  revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(1,6), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,35), key.ylab=NA, main="Senescence expression 
  in OIR vs NORM in p17 retina")
dev.off()



setEPS()
postscript("DoHeatmap_Selection_Senescence_Updated_noCellSEN_P17_OIRvsNORM.eps", height = 5, width = 10)
heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
  Rowv=Rowv,
  #Rowv=TRUE,
  Colv=Colv,
  #Colv=TRUE,
  srtCol=45,
  revC=TRUE, key.title = NA, density.info="none", lhei=c(1,3), lwid=c(1,8), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(7,25), key.ylab=NA, main=paste0("Fold values OIR vs Norm"),
  key = F)
dev.off()
 



svg(paste("DoHeatmap_Selection_HRAS_KRAS_selection_P17_OIRvsNORM_Markers_new_selection.svg", sep="."), width=12, height=6)  
heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
  Rowv=Rowv,
  Colv=Colv,
  srtCol=45,
  revC=TRUE, key.title = NA, density.info="none", lhei=c(1,5), lwid=c(1,6), 
  key.par=list(mar=c(5, 1, 1, 0)), 
  margins=c(10,35), key.ylab=NA)
dev.off()


Selected.genes <- c("GO_ENDOLYSOSOME_MEMBRANE")


for(this_cluster in cell_type_list){

  #####this_cluster <- "Pericytes"

  ## Print status for which identity is being processedcc
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object,
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



#END
q("no")
