# ---------------------------------------------------------------------------
# 01_immune_cells - 04_gsva_prepare.R
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
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA) 
data(c2BroadSets)
library(gplots)
library(heatmap3)
library(viridis)

library(Seurat)

#Set directory of dataset to analyse

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/CellIdentification_5"))

microglia <- readRDS("Microglia_subclustered_renamed.rds")

summary(microglia@ident)
summary(microglia@meta.data)


summary(as.factor(microglia@meta.data$Cond_Sorting))

Seurat_object <- microglia 

# #Imputation

# library(Rmagic)
# library(readr)
# library(phateR)


# Seurat_object.data <- as.data.frame(t(as.matrix(Seurat_object@data)))

# Seurat_object.data <- magic(Seurat_object.data, genes="all_genes")

# #saveRDS(Seurat_object.data, paste(project_name, res, Dim, perp, "imputed.data.rds", sep="."))

# #Seurat_object.data <- readRDS(paste(project_name, res, Dim, perp, "imputed.data.rds", sep="."))

# Seurat_object.data <- data.frame(t(Seurat_object.data$result))

# Seurat_object@imputed <- Seurat_object.data

# saveRDS(Seurat_object, paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))


#Make matrix

setwd(file.path(B_DIR, "Subclustering/ImmuneCells/GSVA"))

#length(x=Seurat_object@var.genes)

#Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name = NULL, ident.use = NULL,
#  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
#  accept.value = NULL, do.center = FALSE, do.scale = FALSE,
#  max.cells.per.ident = 1000, random.seed = 1, do.clean = FALSE)

#exprs <- as.matrix(VascularEndothelium@data[VascularEndothelium@var.genes,])

exprs <- as.matrix(Seurat_object@data)

####subset for testing, optional

#random.cells <- sample(colnames(exprs), 10)

#exprs <- exprs[,random.cells]

##Transform gene symbol into entrez gene ID, OPTIONAL
#require(biomaRt)

#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)
#annots <- getBM(mart=mart, attributes=c("hgnc_symbol", "entrezgene"), filter="hgnc_symbol", values=rownames(exprs), uniqueRows=TRUE)
#annots <- annots[!duplicated(annots[,1]),]
#exprs <- exprs[which(rownames(exprs) %in% annots[,1]),]
#annots <- annots[which(annots[,1] %in% rownames(exprs)),]
#exprs <- exprs[match(annots[,1], rownames(exprs)),]
#rownames(exprs) <- annots[,2]

#exprs <- as.matrix(exprs[!grepl('NA', rownames(exprs)), ])

#anyDuplicated(rownames(exprs))

#rownames(exprs) = make.unique(rownames(exprs))

#exprs <- drop_na(as.data.frame(exprs))

#exprs <- as.matrix(exprs)

#Phenotypic data summarizes information about the samples (e.g., sex, age, and treatment
#status; referred to as ‘covariates’). The information describing the samples can be represented
#as a table with S rows and V columns, where V is the number of covariates. An example of
#phenotypic data can be input with

#pData <- Subset.cells@meta.data

#pData <- Seurat_object@meta.data[random.cells,]

pData <- Seurat_object@meta.data

#Bioconductor’s Biobase package provides a class called AnnotatedDataFrame that conveniently
#stores and manipulates the phenotypic data and its metadata in a coordinated
#fashion. Create and view an AnnotatedDataFrame instance with:


phenoData <- new("AnnotatedDataFrame", data=pData)

#An ExpressionSet object is created by assembling its component parts and callng the ExpressionSet
#constructor:


Full_eset <- ExpressionSet(assayData=exprs, phenoData=phenoData, annotation="org.Hs.eg.db")

#filtered_eset <- nsFilter(Full_eset, require.entrez=FALSE, remove.dupEntrez=TRUE,
#var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE)

#Subset.cells_filtered_eset <- filtered_eset$eset

#exprs.Subset.cells_filtered_eset <- exprs(Subset.cells_filtered_eset) 

exprs.Full_eset <- exprs(Full_eset) 

write.csv(exprs.Full_eset, "exprs.Full_eset.NonImputed.csv")




q("no")
