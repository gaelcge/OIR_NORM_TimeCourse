# ---------------------------------------------------------------------------
# Step 03 - GSVA enrichment scores per cell, MSigDB v7.1 plus curated
#           senescence sets.
#
# ESTABLISHES
#     A gene-set-by-cell score matrix for the four P14/P17 groups, which step 05
#     loads back onto the object as a "GO" assay. Without this file step 05
#     cannot run.
#
# WHY THIS WAY
#     kcdf = "Gaussian" is correct for the input here: GSVA is given
#     log-normalised expression, not integer counts, and the Poisson kernel is
#     for the latter.
#
#     Restricted to NORM/OIR at P14 and P17 because those are two of the three
#     timepoints where both conditions exist; P12 is present in the object but
#     was not included in this run.
#
# NOT REPRODUCIBLE AS WRITTEN - read this before trusting a re-run
#     The gene-set collection is MSigDB h + c2 + c5 (v7.1, checksums in
#     reference/PROVENANCE.txt) COMBINED WITH three curated .gmx files -
#     MALLETTE_SENESCENCE_UP, GLOBAL_SENESCENCE_LITERATURE_CURATED and
#     SAWCHYN_UNBIASED_UP - that live in a collaborator's project directory and
#     are not readable by this account (permission denied, checked 2026-09-29).
#     They are not committed and their contents are not recorded anywhere in
#     this repository. A re-run using only the MSigDB collections will produce a
#     different score matrix from the published one.
#
# REQUIRES
#     step 02 (reads its annotated object)
#
# INPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#     <MSIGDB_DIR>/{h.all,c2.cp,c5.all}.v7.1.symbols.gmt
#     <CUSTOM_GENESET_DIR>/*.gmx                  UNAVAILABLE, see above
#
# OUTPUTS
#     <GSVA_DIR>/gsva.exprs.MsigDB_h_c2_c5.senescence.all.v7.1.symbols.csv
#
# RUNTIME
#     Requested 48 h across 26 tasks (03_gsva_senescence.slurm.sh). The original
#     wrapper asks for --mem-per-cpu=150G on 26 tasks, which is ~3.9 TB in
#     total; that is almost certainly not what was intended and the request is
#     reproduced here only as a record of what was submitted.
#
# USAGE
#     Rscript 03_gsva_senescence.R
# ---------------------------------------------------------------------------

# ===========================================================================
# 0. parameters
# ===========================================================================
# ---------------------------------------------------------------------------
# cluster locations. The merged Drop-seq matrix sits on the ctb-jsjoyal storage
# allocation; every other input is on def-jsjoyal. Both are read-only here.
# Point these two at your own copies and nothing else in the file changes.
# ---------------------------------------------------------------------------
SEQ_DIR  <- "/project/ctb-jsjoyal/gaelcge/Sequencing/Merging/TimeCourseOIR"
PROJ_DIR <- "/project/def-jsjoyal/gaelcge"

OBJ_DIR  <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering")
GSVA_DIR <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR")

GENE_LIST_DIR     <- file.path(PROJ_DIR, "Gene_list")
MSIGDB_DIR        <- file.path(GENE_LIST_DIR, "Gmt.file/MsigDB/V7.1")
PATHWAY_DIR       <- file.path(GENE_LIST_DIR, "Gmt.file/Pathways")
# Three curated senescence .gmx files held by a collaborator and NOT readable
# by this account (checked 2026-09-29). See reference/PROVENANCE.txt: this step
# cannot be reproduced without them.
CUSTOM_GENESET_DIR <- file.path(PROJ_DIR, "..", "jhowa105/projects/Mike")

# object identity, interpolated into every output filename by the original
# scripts. res=1 is the clustering resolution actually used - confirmed from the
# object, whose clustering column is `integrated_snn_res.1`.
project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"
res  <- 1
Dim  <- 20
DIM_nb <- c(1:20)
perp <- 30

GSVA_METHOD <- "gsva"
GSVA_KCDF   <- "Gaussian"   # input is log-normalised, not counts

library(GSVA)
library(GSEABase)
library(parallel)
library(sigPathway)
library(Seurat)
library(limma)
library(gskb)
library(future)
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(GSVAdata)
data(c2BroadSets)
library(gplots)
library(heatmap3)
library(viridis)
library(msigdbr)
library(org.Hs.eg.db)


##Import matrix
setwd(OBJ_DIR)

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"

res = 1

DIM_nb <- c(1:20)

Dim <- 20

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Cond_TimePoint %in% c("NORM_P14", "NORM_P17","OIR_P14", "OIR_P17")))

Seurat_object <- SubsetData(Seurat_object, cells = Subseted_cells)

Seurat_object

exp_mat <- GetAssayData(object = Seurat_object, assay = "RNA", slot = "data")

rownames(exp_mat) <- toupper(rownames(exp_mat) )

#exp_mat <- exp_mat[,1:500]


##IMport your own geneset

geneset1 <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "MALLETTE_SENESCENCE_UP.gmx"), verbose = TRUE)

geneset1 <- GeneSet(geneset1[[1]]$probes, geneIdType=SymbolIdentifier(), setName="MALLETTE_SENESCENCE_UP")

geneset2 <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "GLOBAL_SENESCENCE_LITERATURE_CURATED_2020.gmx"), verbose = TRUE)

geneset2 <- GeneSet(geneset2[[1]]$probes, geneIdType=SymbolIdentifier(), setName="GLOBAL_SENESCENCE_LITERATURE_CURATED")

geneset3 <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "SAWCHYN_UNBIASED_UP_v4.gmx"), verbose = TRUE)

geneset3 <- GeneSet(geneset3[[1]]$probes, geneIdType=SymbolIdentifier(), setName="SAWCHYN_UNBIASED_UP")


##Use Geneset from MsigDB gmt file
gmt.MsigDB.h <- getGmt(file.path(MSIGDB_DIR, "h.all.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.c2 <- getGmt(file.path(MSIGDB_DIR, "c2.cp.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.c5 <- getGmt(file.path(MSIGDB_DIR, "c5.all.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())


#Use Geneset from R msigdbr tool

h_gene_sets = msigdbr(species = "Mus musculus", category = "H")

### Make collection of geneset

gsc_final <- GeneSetCollection(c(gmt.MsigDB.h, gmt.MsigDB.c2,gmt.MsigDB.c5,geneset1,geneset2,geneset3))


###Run GSVA
setwd(GSVA_DIR)

parallel.sz <- availableCores()-1

#Choose a name for your output file
out.file <- "gsva.exprs.MsigDB_h_c2_c5.senescence.all.v7.1.symbols.csv"


#Run GSVA
print(paste0('parallel.sz = ', parallel.sz))

gsva.mat <- gsva(as.matrix(exp_mat), gsc_final, method='gsva',kcdf='Gaussian', parallel.sz = parallel.sz, verbose = TRUE, parallel.type="SOCK")

#gsva.mat <- gsva(as.matrix(exp_mat), gmt, method='gsva',kcdf='Gaussian', parallel.sz = parallel.sz, verbose = TRUE, parallel.type="SOCK")

message(gsva.mat[1:5,1:5])

##Write GSVA results

write.table(gsva.mat, file = out.file, sep = ",", row.names = TRUE, col.names = TRUE)

q("no")


