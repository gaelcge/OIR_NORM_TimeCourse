# ---------------------------------------------------------------------------
# Step 01 - Raw UMI count matrix for P14/P17, and the first QC view of it.
#
# ESTABLISHES
#     The count matrix deposited at GEO as GSE150703, and the initial per-batch
#     QC violin over unfiltered data.
#
# WHY THIS WAY
#     This step exists to produce a deposited artifact, not to feed the
#     analyses: neither 1_timecourse_p5_p17 nor 2_compartments_p14_p17 reads its
#     output. Both read the merged matrix directly. It is numbered 01 because it
#     is where a reader starting from the public deposit enters the project.
#
# NOT REPRODUCIBLE - read before re-running
#     The McCarroll-lab columns are reduced with
#     sample(names(DGE.data.Mccarrol), 6000) and NO SEED IS SET. The exact 6000
#     cells that went into the deposited matrix cannot be recovered from this
#     code. They are recoverable from the deposited matrix itself, which is why
#     the GEO record rather than this script is the authoritative definition of
#     the P14/P17 cell set.
#
#     This affects only this step. The two analysis folders read all McCarroll
#     cells from the merged matrix (37,766 Joyal and 13,385 McCarroll in the
#     P5-P17 object), so neither depends on this subsample.
#
# SEURAT VERSION
#     Seurat v2 API throughout - raw.data, min.genes, MakeSparse(),
#     features.plot, x.lab.rot. None of these exist in Seurat 3 or later. See
#     SOFTWARE_VERSIONS.md.
#
# INPUTS
#     <SEQ_DIR>/Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt     4.67 GB
#
# OUTPUTS
#     retina_NORM_OIR_P14_P17_C57_WR_CD73FT_RawUMI_Count_DGEmatrix.txt
#     VlnPlotQC.Cond_Sorting.Init.png
#
# USAGE
#     Rscript 01_make_dge_matrix.R
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
MERGED_DGE <- file.path(SEQ_DIR, "Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt")

N_MCCARROLL_CELLS <- 6000   # drawn WITHOUT a seed - see the header
MIN_CELLS         <- 3
MIN_GENES         <- 100

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"
res  <- 3
DIM_nb <- c(1:17)
Dim  <- 17
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
library(DoubletFinder)
library(harmony)

# Load the dataset


project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30


DGE.data=data.frame(fread(MERGED_DGE, sep="\t", header=TRUE), row.names=1)

# Subset dataset and matrix

#retina.data <- retina.data[ ,sample(names(retina.data), 500)]

DGE.data.Mccarrol <- DGE.data %>% dplyr:: select(grep("Mccarrol", names(DGE.data)))

DGE.data.Mccarrol <- DGE.data.Mccarrol[ ,sample(names(DGE.data.Mccarrol), 6000)]

corner(DGE.data.Mccarrol)
dim(DGE.data.Mccarrol)

DGE.data.CHUSJ <- DGE.data %>% dplyr:: select(grep("Joyal", names(DGE.data)))

corner(DGE.data.CHUSJ)
dim(DGE.data.CHUSJ)

DGE.data <- merge(DGE.data.Mccarrol, DGE.data.CHUSJ, by=0, all=TRUE)

rownames(DGE.data) <- DGE.data$Row.names


DGE.data <- DGE.data %>% dplyr:: select(grep("NORM", names(DGE.data)), grep("OIR", names(DGE.data)))

colnames(DGE.data) = gsub(".", "_", colnames(DGE.data), fixed = TRUE)


###Write DGE 

setwd(file.path(TC_DIR, ""))

DGE.data_OIR_NORM_P14_P17 <- DGE.data %>% dplyr:: select(grep("NORM_P14", names(DGE.data)), 
                                        grep("NORM_P17", names(DGE.data)), 
                                        grep("OIR_P14", names(DGE.data)), 
                                        grep("OIR_P17", names(DGE.data))
                                        )

df <- as.matrix(DGE.data_OIR_NORM_P14_P17)

head(colnames(df))

head(rownames(df))

write.table(df, "retina_NORM_OIR_P14_P17_C57_WR_CD73FT_RawUMI_Count_DGEmatrix.txt", sep="\t",
                                      quote=FALSE,
                                      row.names=TRUE,
                                      col.names=TRUE)


# Look at the data matrix
corner(DGE.data)
dim(DGE.data)


###Optional = downsize to 2000

#DGE.data <- DGE.data[ ,sample(names(DGE.data), 2000)]

#Initialize the Seurat object with the raw (non-normalized data)
# Keep all genes expressed in >= 3 cells, keep all cells with >= 100 genes
Seurat_object <- CreateSeuratObject(raw.data=DGE.data_OIR_NORM_P14_P17, min.cells = 3, min.genes = 100, project = project_name)

class(x = Seurat_object@raw.data)
Seurat_object <- MakeSparse(object = Seurat_object)
class(x = Seurat_object@raw.data)


#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_object@data), value = TRUE)
percent.mito <- Matrix::colSums(Seurat_object@raw.data[mito.genes, ])/Matrix::colSums(Seurat_object@raw.data)

CRYA.genes <- grep(pattern = "^CRYA", x = rownames(x = Seurat_object@data), value = TRUE)
CRYB.genes <- grep(pattern = "^CRYB", x = rownames(x = Seurat_object@data), value = TRUE)
CRYG.genes <- grep(pattern = "^CRYG", x = rownames(x = Seurat_object@data), value = TRUE)
crystal.genes <- c(CRYA.genes, CRYB.genes)
percent.crystal <- Matrix::colSums(Seurat_object@raw.data[crystal.genes, ])/Matrix::colSums(Seurat_object@raw.data)

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.mito, col.name = "percent.mito")

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal, col.name = "percent.crystal")

#Assigning batch to data
batch_long <- data.frame("merged_name"=colnames(DGE.data_OIR_NORM_P14_P17))

batch_assigned <- batch_long %>% separate(merged_name,c("Condition", "TimePoint", "Sorting", "Labo", "Replicate", "CellBarCode"),sep="_")

batch_assigned <- batch_assigned %>%
  mutate("Batch"=paste(Condition,TimePoint,Sorting,Replicate,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Cond_Sorting"=paste(Condition,Sorting,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
rownames(batch_assigned) <- batch_long$merged_name

Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned, col.name = "Batch")


png(filename="VlnPlotQC.Cond_Sorting.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features.plot = c("nGene"), group.by = "Batch",  nCol = 2, x.lab.rot = TRUE)
dev.off()


q("no")
