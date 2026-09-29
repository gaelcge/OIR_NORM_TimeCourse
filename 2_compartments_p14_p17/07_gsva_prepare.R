# ---------------------------------------------------------------------------
# Step 07 - Expression matrices prepared for GSVA.
#
# ESTABLISHES
#     The per-cell expression matrices that steps 08 and 09 score. Also the
#     1000-cells-per-identity subsample that the 'Subset1000CellperIdent'
#     outputs are built from.
#
# WHY THIS WAY
#     MAGIC imputation is present in this script and COMMENTED OUT. That is
#     what 'NotImputed' means in every downstream path: imputation was tried
#     and not used. Scoring gene sets on imputed values inflates within-cell
#     correlation between genes of the same set, which is exactly the quantity
#     GSVA measures - so the un-imputed matrix is the conservative input.
#
#     Subsampling to 1000 cells per identity equalises cell types whose sizes
#     differ by two orders of magnitude, so that a gene-set score distribution
#     is not dominated by rods.
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
#
# OUTPUTS
#     expression matrices under <B_DIR>/GSVA/NotImputed/ and
#     <B_DIR>/GSVA/NotImputed/Subset1000CellperIdent/
#
# USAGE
#     Rscript 07_*.R
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

setwd(file.path(B_DIR, "CellTypeAnnotation/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))

summary(Seurat_object@ident)
summary(Seurat_object@meta.data)



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

dir.create(file.path(B_DIR, "GSVA/NotImputed/Subset1000CellperIdent"))

setwd(file.path(B_DIR, "GSVA/NotImputed/Subset1000CellperIdent"))

length(x=Seurat_object@var.genes)

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name = NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  accept.value = NULL, do.center = TRUE, do.scale = TRUE,
  max.cells.per.ident = 1000, random.seed = 1, do.clean = FALSE)

#exprs <- as.matrix(VascularEndothelium@data[VascularEndothelium@var.genes,])

exprs <- as.matrix(Seurat_object@scale.data)

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

write.csv(exprs.Full_eset, "exprs.Full_eset.NonImputed.scaled.csv")

write.csv(exprs, "exprs.Full_eset.NonImputed.scaled.csv")



q("no")
