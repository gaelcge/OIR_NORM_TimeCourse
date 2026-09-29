# ---------------------------------------------------------------------------
# Step 08 - GSVA senescence scores, P17, all retained cell types.
#
# ESTABLISHES
#     Per-cell senescence gene-set scores at P17.
#
# WHY THIS WAY
#     kcdf='Gaussian' is correct here: the input is log-normalised expression,
#     not counts. Restricted to P17, the timepoint at which the senescence
#     question is asked.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 07
#
# INPUTS
#     the prepared expression matrix
#     senescence gene sets (see PROVENANCE)
#
# OUTPUTS
#     <B_DIR>/GSVA/NotImputed/retina_p17_data_gsva_out_senescence_May2020.csv
#
# USAGE
#     Rscript 08_*.R
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

library(GSVA)
library(GSEABase)
library(parallel)
library(sigPathway)
library(Seurat)
library(limma)


setwd(file.path(B_DIR, "CellTypeAnnotation/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))

retina <- UpdateSeuratObject(Seurat_object)
retina <- AddMetaData(retina, paste(retina@meta.data$cell_type, retina@meta.data$Condition, sep="_"), col.name = "cell_type_Condition")
retina <- subset(retina, subset = cell_type == 'Bipolar cells' | cell_type == 'Rods' | cell_type=='Cones' | cell_type=='Muller glia' | cell_type=='Amacrine cells' | cell_type == 'Pericytes' | cell_type == 'Retinal ganglion cells' | cell_type == 'Astrocytes' | cell_type == 'Horizontal cells' | cell_type == 'Endothelial cells' | cell_type == 'Immune cells')

p17 <- subset(retina, subset = TimePoint == 'P17')

p17 <- NormalizeData(p17)

exp_mat <- as.matrix(p17@assays$RNA@data)

###Run GSVA

setwd(file.path(B_DIR, "GSVA/NotImputed"))


gmt <- getGmt(file.path(CUSTOM_GENESET_DIR, "senescence_gene_sets.gmt"), geneIdType=SymbolIdentifier())

#geneset <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "MALLETTE_SENESCENCE_UP.gmx"), verbose = TRUE)

#geneset <- GeneSet(geneset[[1]]$probes, geneIdType=SymbolIdentifier(), setName="MALLETTE_SENESCENCE_UP")

#gsc <- GeneSetCollection(geneset)

#gsc_final <- GeneSetCollection(c(gsc, gmt))

parallel.sz <- 15

out.file <- "retina_p17_data_gsva_out_senescence_May2020.csv"

gsva.mat <- gsva(exp_mat, gmt, method='gsva',kcdf='Gaussian', parallel.sz = parallel.sz, verbose = TRUE)

message(gsva.mat[1:5,1:5])

write.table(gsva.mat, file = out.file, sep = ",", row.names = TRUE, col.names = TRUE)

q("no")

