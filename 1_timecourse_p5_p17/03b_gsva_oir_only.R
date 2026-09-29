# ---------------------------------------------------------------------------
# Step 03b - GSVA on the OIR cells only, MSigDB collections alone.
#
# ESTABLISHES
#     Nothing that any later step consumes. Kept because it ran and its output
#     exists on the cluster, so a reader comparing the repository against the
#     directory would otherwise find an unexplained file.
#
# WHY THIS WAY / HOW IT DIFFERS FROM STEP 03
#     This is not a faster variant of step 03 despite the directory it sat in
#     being named Cpu. It is a different analysis:
#       - cells:     OIR_P14 and OIR_P17 only, no normoxic comparison
#       - gene sets: MSigDB h + c2 + c5 only, without the three curated
#                    senescence sets
#       - output:    gsva.exprs.MsigDB.csv, which nothing reads
#     Step 05 loads the step 03 output, not this one.
#
# REQUIRES
#     step 02
#
# INPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#     <MSIGDB_DIR>/{h.all,c2.cp,c5.all}.v7.1.symbols.gmt
#
# OUTPUTS
#     <GSVA_CPU_DIR>/gsva.exprs.MsigDB.csv        not consumed downstream
#
# USAGE
#     Rscript 03b_gsva_oir_only.R
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
GSVA_CPU_DIR <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR/Cpu")

OBJ_DIR  <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering")
GSVA_DIR <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR")
MSIGDB_DIR <- file.path(PROJ_DIR, "Gene_list/Gmt.file/MsigDB/V7.1")

# object identity, interpolated into every output filename by the original
# scripts. res=1 is the clustering resolution actually used - confirmed from the
# object, whose clustering column is `integrated_snn_res.1`.
project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"
res  <- 1
Dim  <- 20
DIM_nb <- c(1:20)
perp <- 30

GSVA_METHOD <- "gsva"
GSVA_KCDF   <- "Gaussian"

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


##Import matrix
setwd(OBJ_DIR)

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"

res = 1

DIM_nb <- c(1:20)

Dim <- 20

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Cond_TimePoint %in% c("OIR_P14", "OIR_P17")))

Seurat_object <- SubsetData(Seurat_object, cells = Subseted_cells)


exp_mat <- GetAssayData(object = Seurat_object, assay = "RNA", slot = "data")

rownames(exp_mat) <- toupper(rownames(exp_mat) )

#exp_mat <- exp_mat[,1:500]


gmt.MsigDB.h <- getGmt(file.path(MSIGDB_DIR, "h.all.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.c2 <- getGmt(file.path(MSIGDB_DIR, "c2.cp.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.c5 <- getGmt(file.path(MSIGDB_DIR, "c5.all.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gsc_final <- GeneSetCollection(c(gmt.MsigDB.h, gmt.MsigDB.c2,gmt.MsigDB.c5))


###Run GSVA
setwd(GSVA_CPU_DIR)

parallel.sz <- availableCores()-1

out.file <- "gsva.exprs.MsigDB.csv"

print(paste0('parallel.sz = ', parallel.sz))

gsva.mat <- gsva(as.matrix(exp_mat), gsc_final, method='gsva',kcdf='Gaussian', parallel.sz = parallel.sz, verbose = TRUE, parallel.type="SOCK")

#gsva.mat <- gsva(as.matrix(exp_mat), gmt, method='gsva',kcdf='Gaussian', parallel.sz = parallel.sz, verbose = TRUE, parallel.type="SOCK")

message(gsva.mat[1:5,1:5])

write.table(gsva.mat, file = out.file, sep = ",", row.names = TRUE, col.names = TRUE)

q("no")


