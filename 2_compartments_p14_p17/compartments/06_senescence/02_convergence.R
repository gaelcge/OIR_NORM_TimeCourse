# ---------------------------------------------------------------------------
# 06_senescence - 02_convergence.R
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


setwd(file.path(B_DIR, "CellTypeAnnotation/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))

summary(as.factor(Seurat_object@meta.data$Cond_Sorting))


####loading Subcluster

### April 7 2020

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/CellIdentification_5"))

Microglia <- readRDS("Microglia_subclustered_renamed.rds")

summary(Microglia@ident)

Seurat_object <- SetAllIdent(Seurat_object, id = "cell_type")

new_cell <- WhichCells(Microglia, ident = "Cluster 1")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Immune cells 1")

new_cell <- WhichCells(Microglia, ident = "Cluster 2")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Immune cells 2")

new_cell <- WhichCells(Microglia, ident = "Cluster 3")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Immune cells 3")

new_cell <- WhichCells(Microglia, ident = "Cluster 4")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Immune cells 4")

new_cell <- WhichCells(Microglia, ident = "Cluster 5")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Immune cells 5")

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name = NULL, ident.use = NULL,
  ident.remove = "Immune cells")


### Subset data

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name = NULL, 
				ident.use = c("Immune cells 1", "Immune cells 2", "Immune cells 3", "Immune cells 4",
								"Immune cells 5", "Astrocytes", "Endothelial cells", "Pericytes", "Muller glia"))

##Save combine ident

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/Mapping"))

saveRDS(Seurat_object, paste(project_name, res, Dim, perp, "Seurat_object.annotated.subseted.rds", sep="."))

