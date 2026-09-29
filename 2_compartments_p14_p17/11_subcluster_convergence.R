# ---------------------------------------------------------------------------
# Step 11 - Compartment sub-labels merged back onto the whole retina.
#
# ESTABLISHES
#     One object carrying the fine labels from every compartment, so that
#     subtypes discovered separately can be compared on a single embedding.
#
# WHY THIS WAY
#     Numbered last because it CONSUMES the compartments. Each compartment was
#     subclustered independently - a shared clustering at a resolution fine
#     enough to split microglial states would have shattered the rods - and
#     this step is where those independent label sets are reconciled.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 03, and every compartment under compartments/
#
# INPUTS
#     <B_DIR>/CellTypeAnnotation/<project_name>.<res>.<Dim>.<perp>.Seurat_object.annotated.rds
#     compartments/*/ subclustered objects (NeuroGlia_subclustered_renamed.rds,
#         Microglia_subclustered_renamed.rds, and the rest)
#
# OUTPUTS
#     the converged object and its figures under <B_DIR>/Subclustering/
#
# USAGE
#     Rscript 11_*.R
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

### 2019

setwd(file.path(B_DIR, "Subclustering/NeuroGlialCell/CellIdentification"))

NeuroGlia <- readRDS("NeuroGlia_subclustered_renamed.rds")

summary(NeuroGlia@ident)

new_cell <- WhichCells(NeuroGlia, ident = "Activated Muller Glia")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Activated Muller Glia")

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name = NULL, ident.use = NULL,
  ident.remove = "Muller glia")

### April 7 2020

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/CellIdentification_5"))

Microglia <- readRDS("Microglia_subclustered_renamed.rds")

summary(Microglia@ident)

new_cell <- WhichCells(Microglia, ident = "Activated Muller Glia")

Seurat_object <- SetIdent(Seurat_object, cells.use = new_cell, ident.use = "Activated Muller Glia")

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name = NULL, ident.use = NULL,
  ident.remove = "Muller glia")

##Save combine ident

setwd(file.path(B_DIR, "CellTypeAnnotation/"))

saveRDS(Seurat_object, paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))



#END
q("no")
