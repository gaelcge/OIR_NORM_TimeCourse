# ---------------------------------------------------------------------------
# Step 01 - Integrated Seurat object from the merged Drop-seq matrix.
#
# ESTABLISHES
#     The object every other step in this folder reads: 51,151 cells x 21,705
#     features spanning P5, P7, P10, P12, P14 and P17, integrated across the two
#     sorted fractions (WR and Cd73ft) by SCT anchors.
#
# WHY THIS WAY
#     Integration is anchored on `Sorting`, not on condition or timepoint. The
#     fractions are the technical split that has to be reconciled before cells
#     can be compared; anchoring on condition would regress out the OIR effect
#     this study is measuring. Note the consequence in the design: P7 and P10
#     were sequenced as WR only, so those timepoints enter the integration
#     through a single fraction.
#
#     percent.crystal < 0.025 is an unusual QC filter and it is deliberate.
#     Dissected retina carries lens fragments, and crystallin-high cells are
#     contamination rather than a retinal population. CRYA and CRYB genes are
#     summed for the filter; CRYG genes are collected and then not used, which
#     is a dead line left in place because this script produced published
#     numbers.
#
#     Cell-cycle S and G2M scores are regressed out inside SCTransform rather
#     than corrected afterwards, so the integration anchors are computed on
#     residuals that are already free of cycle structure - the retina at P5-P17
#     is still proliferating and cycle would otherwise dominate the early
#     timepoints.
#
# REQUIRES
#     nothing in this repository. The merged matrix is built upstream of it;
#     see the root README for the Drop-seq alignment and merge.
#
# INPUTS
#     <SEQ_DIR>/Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt     4.67 GB, tab-delimited
#     <CC_GENES>                                            98 lines, see PROVENANCE
#
# OUTPUTS
#     <OBJ_DIR>/DGE.data.sparse.<project_name>.rds
#     <OBJ_DIR>/<project_name>.Seurat_object.integrated.rds
#     VlnPlotQC.{Cond_Sorting,Conditions,Sorting}.Init.png    pre-filter QC
#
# RUNTIME
#     ~7 h on 16 cores / 180 GB (see 01_integrate_by_sorting.slurm.sh).
#
# USAGE
#     Rscript 01_integrate_by_sorting.R
#     To resume from the saved object instead, start at step 02.
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
MERGED_DGE <- file.path(SEQ_DIR, "Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt")
CC_GENES   <- file.path(PROJ_DIR, "Seurat_ressource/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

# QC thresholds, applied in one subset() call
MIN_FEATURES  <- 100
MAX_FEATURES  <- 6000
MAX_COUNTS    <- 10000
MAX_PCT_MITO  <- 0.10
MAX_PCT_CRYST <- 0.025

N_INT_FEATURES <- 3000   # SelectIntegrationFeatures

# object identity, interpolated into every output filename by the original
# scripts. res=1 is the clustering resolution actually used - confirmed from the
# object, whose clustering column is `integrated_snn_res.1`.
project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"
res  <- 1
Dim  <- 20
DIM_nb <- c(1:20)
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
#library(DoubletFinder)
library(harmony)
library(sctransform)
library(future)
plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 50000 * 1024^2)

# Load the dataset

setwd(OBJ_DIR)

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"

res = 1

DIM_nb <- c(1:20)

Dim <- 20

perp = 30


DGE.data=data.frame(fread(MERGED_DGE, sep="\t", header=TRUE), row.names=1)

# Subset dataset and matrix

#retina.data <- retina.data[ ,sample(names(retina.data), 500)]

DGE.data.Mccarrol <- DGE.data %>% dplyr:: select(grep("Mccarrol_r3", names(DGE.data)),grep("Mccarrol_r5", names(DGE.data)))

#write.table(DGE.data.Mccarrol, "Dropseq_p14_retina_Mccarroll.txt", quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)


#DGE.data.Mccarrol <- DGE.data.Mccarrol[ ,sample(names(DGE.data.Mccarrol), 6000)]

corner(DGE.data.Mccarrol)
dim(DGE.data.Mccarrol)

DGE.data.CHUSJ <- DGE.data %>% dplyr:: select(grep("Joyal", names(DGE.data)))

#write.table(DGE.data.CHUSJ, "Dropseq_p14_retina_Joyal.txt", quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)

#DGE.data.CHUSJ_OIR <- DGE.data %>% dplyr:: select(grep("OIR.P14.WR.Joyal.r1", names(DGE.data)))

#write.table(DGE.data.CHUSJ_OIR, "Dropseq_p14_retina_retinopathy_Joyal.txt", quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)

corner(DGE.data.CHUSJ)
dim(DGE.data.CHUSJ)

DGE.data <- merge(DGE.data.Mccarrol, DGE.data.CHUSJ, by=0, all=TRUE)

rownames(DGE.data) <- DGE.data$Row.names


DGE.data <- DGE.data %>% dplyr:: select(grep("NORM", names(DGE.data)), grep("OIR", names(DGE.data)))

colnames(DGE.data) = gsub(".", "_", colnames(DGE.data), fixed = TRUE)

# Look at the data matrix
corner(DGE.data)
dim(DGE.data)

#Save the P14 and P17 dataset into sparse matrix
DGE.data_P14_P17 <- DGE.data %>% dplyr:: select(grep("P14", names(DGE.data)), grep("P17", names(DGE.data)))
DGE.sparse <- Matrix(as.matrix(DGE.data_P14_P17), sparse = TRUE)
saveRDS(DGE.sparse, "DGE.data.sparse.RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting.rds")


###Optional = downsize to 2000

#DGE.data <- DGE.data[ ,sample(names(DGE.data), 2000)]

#Initialize the Seurat object with the raw (non-normalized data)
# Keep all genes expressed in >= 3 cells, keep all cells with >= 100 genes
Seurat_object <- CreateSeuratObject(DGE.data, project = paste0(project_name), assay = "RNA",
  min.cells = 3, min.features = 100)



#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
percent.mito <- Matrix::colSums(Seurat_object@assays$RNA[mito.genes, ])/Matrix::colSums(Seurat_object@assays$RNA)

CRYA.genes <- grep(pattern = "^CRYA", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
CRYB.genes <- grep(pattern = "^CRYB", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
CRYG.genes <- grep(pattern = "^CRYG", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
crystal.genes <- c(CRYA.genes, CRYB.genes)
percent.crystal <- Matrix::colSums(Seurat_object@assays$RNA[crystal.genes, ])/Matrix::colSums(Seurat_object@assays$RNA)


Seurat_object[["percent.crystal"]] <- PercentageFeatureSet(Seurat_object, pattern = crystal.genes)

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.mito, col.name = "percent.mito")

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal, col.name = "percent.crystal")

#Assigning batch to data
batch_long <- data.frame("merged_name"=colnames(DGE.data))

batch_assigned <- batch_long %>% separate(merged_name,c("Condition", "TimePoint", "Sorting", "Labo", "Replicate", "CellBarCode"),sep="_")

batch_assigned <- batch_assigned %>%
  mutate("Batch"=paste(Condition,TimePoint,Sorting,Replicate,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Cond_Sorting"=paste(Condition,Sorting,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
rownames(batch_assigned) <- batch_long$merged_name

Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned)

png(filename="VlnPlotQC.Cond_Sorting.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Cond_Sorting")
dev.off()

png(filename="VlnPlotQC.Conditions.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Condition")
dev.off()

png(filename="VlnPlotQC.Sorting.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Sorting")
dev.off()

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well

plot1 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.png", width=900, height=600, bg = "white", res = 50)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#We filter out cells that have unique gene counts over 6000
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)


Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA < 10000 & percent.mito < 0.10 & percent.crystal < 0.025)


#Cell Cycle scoring

cc.genes <- readLines(con=CC_GENES)

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]


Seurat_object <- CellCycleScoring(Seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

### Subset for testing

#Seurat_object <- SubsetData(Seurat_object, assay = NULL, cells = NULL,
#  subset.name = NULL, ident.use = NULL, ident.remove = NULL,
#  low.threshold = -Inf, high.threshold = Inf, accept.value = NULL,
#  max.cells.per.ident = 50, random.seed = 1)

###

Seurat_object.list <- SplitObject(Seurat_object, split.by = "Sorting")

for (i in 1:length(Seurat_object.list)) {
    Seurat_object.list[[i]] <- SCTransform(Seurat_object.list[[i]], verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch", "S.Score", "G2M.Score"))
}

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list, nfeatures = 3000)

Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list, anchor.features = Seurat_object.features, 
    verbose = FALSE)

#Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set  normalization.method = 'SCT':

Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list, normalization.method = "SCT", 
    anchor.features = Seurat_object.features, verbose = FALSE, dims = DIM_nb)

Seurat_object.integrated <- IntegrateData(anchorset = Seurat_object.anchors, normalization.method = "SCT", 
    verbose = FALSE)

plot1 <- FeatureScatter(Seurat_object.integrated, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(Seurat_object.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.integrated.png", width=900, height=600, bg = "white", res = 50)
CombinePlots(plots = list(plot1, plot2))
dev.off()


saveRDS(Seurat_object.integrated, paste(project_name, "Seurat_object.integrated.rds", sep="."))



q("no")
