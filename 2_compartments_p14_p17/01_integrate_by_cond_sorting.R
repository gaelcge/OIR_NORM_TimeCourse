# ---------------------------------------------------------------------------
# Step 01 - Integrated object, anchored on condition x sorted fraction.
#
# ESTABLISHES
#     The P14/P17 object this whole folder reads: CCA integration of the four
#     Condition x Sorting groups, without the variance-ratio cell filter that
#     step 01b applies.
#
# WHY THIS WAY
#     Integration is grouped by Cond_Sorting - condition crossed with sorted
#     fraction - rather than by fraction alone. This is the choice that most
#     distinguishes this analysis from the later P5-P17 one, which anchors on
#     Sorting only. Anchoring on a grouping that includes condition aligns the
#     OIR and normoxic subspaces to each other, which makes cell types
#     comparable across condition at the cost of shrinking the condition
#     effect the study is measuring - the reason the later analysis dropped it.
#
#     'noDiscarding' names what it does not do: step 01b additionally removes
#     cells whose variance is better explained by PCA than by CCA. Downstream
#     steps read THIS object, not that one.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     nothing in this repository; the merged matrix is built upstream.
#
# INPUTS
#     <SEQ_DIR>/Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt
#     <CC_GENES>
#
# OUTPUTS
#     <B_DIR>/Clustering_noDiscarding/<project_name>.Seurat_object.rds
#
# USAGE
#     Rscript 01_*.R
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
CC_GENES   <- file.path(PROJ_DIR, "Seurat_ressource/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

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
library(DoubletFinder)
library(harmony)

# Load the dataset

setwd(file.path(B_DIR, "Clustering_noDiscarding/"))

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

# Look at the data matrix
corner(DGE.data)
dim(DGE.data)


###Optional = downsize to 2000

#DGE.data <- DGE.data[ ,sample(names(DGE.data), 2000)]

#Initialize the Seurat object with the raw (non-normalized data)
# Keep all genes expressed in >= 3 cells, keep all cells with >= 100 genes
Seurat_object <- CreateSeuratObject(raw.data=DGE.data, min.cells = 3, min.genes = 100, project = project_name)

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
batch_long <- data.frame("merged_name"=colnames(DGE.data))

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
VlnPlot(object = Seurat_object, features.plot = c("nGene", "nUMI", "percent.mito", "percent.crystal"), group.by = "Cond_Sorting",  nCol = 2, x.lab.rot = TRUE)
dev.off()

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
png(filename="percentMitoPlot.png")
GenePlot(object=Seurat_object, gene1 = "nUMI", gene2 = "percent.mito")
dev.off()
png(filename="nGenePlot.png")
GenePlot(object=Seurat_object, gene1 ="nUMI", gene2 ="nGene")
dev.off()

#We filter out cells that have unique gene counts over 6000
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)

Seurat_object <- FilterCells(object = Seurat_object, subset.names = c("nGene", "percent.mito", "nUMI", "percent.crystal"), 
    low.thresholds = c(100, -Inf, -Inf, -Inf), high.thresholds = c(6000, 0.10, 10000, 0.025))

# Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)

Seurat_object <- NormalizeData(object = Seurat_object, normalization.method = "LogNormalize", 
    scale.factor = 10000)

# Detection of variable genes across the single cells
# Seurat calculates highly variable genes and focuses on these for downstream analysis. MeanVarPlot(), which works by calculating the average expression and dispersion for each gene, placing these genes into bins, and then calculating a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko et al.), but new methods for variable gene expression identification are coming soon. We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here identify ~2,000 variable genes, and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.

#Regress out unwanted sources of variation
# Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. Seurat implements a basic version of this by constructing linear models to predict gene expression based on user-defined variables. Seurat stores the z-scored residuals of these models in the scale.data slot, and they are used for dimensionality reduction and clustering.
# We typically regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a ‘cell-cycle’ score (as in Macosko et al) and Regress this out as well. Here, we simply regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content an example. note that this overwrites @scale.data. Therefore, if you intend to use RegressOut, you can set do.scale=F and do.center=F in the original object to save some time.
#item{model.use}{Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'}


#Cell Cycle scoring

cc.genes <- readLines(con=CC_GENES)

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]

Seurat_object <- CellCycleScoring(Seurat_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)


##Regression

Seurat_object <- ScaleData(object = Seurat_object, vars.to.regress = c("nUMI", "percent.mito", "Batch", "S.Score", "G2M.Score"))


png(filename="percentMitoPlot.regressed.png")
GenePlot(object=Seurat_object, gene1 = "nUMI", gene2 = "percent.mito")
dev.off()
png(filename="nGenePlot.regressed.png")
GenePlot(object=Seurat_object, gene1 ="nUMI", gene2 ="nGene")
dev.off()


###Subset the different dataset

test <- unique(Seurat_object@meta.data$Cond_Sorting)

for (i in test) {
	subset_cells <- rownames(subset(Seurat_object@meta.data, Cond_Sorting %in% i))
	subset_retina <- SubsetData(Seurat_object, cells.use = subset_cells, subset.name = NULL, ident.use = NULL,
  	ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  	do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
  	random.seed = 1)
	subset_retina <- FindVariableGenes(object = subset_retina, do.plot = T, display.progress = F, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = 0.5)
	print(length(x=subset_retina@var.genes))
  assign(paste("subset_retina",i,sep="_"), subset_retina)
}


# Determine genes to use for CCA, must be highly variable in at least 2 datasets


ob.list <- list()

for (i in test) {
  vecname<-paste("subset_retina",i,sep="_")
  ob.list[[i]]<-get(vecname)
	}

genes.use <- c()

for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 3000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

length(genes.use)

#Perform linear dimensional reduction
#Perform PCA on the scaled data. By default, the genes in object@var.genes are used as input, but can be defined using pc.genes. We have typically found that running dimensionality reduction on genes with high-dispersion can improve performance. However, with UMI data - particularly after using RegressOut, we often see that PCA returns similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome.
Seurat_object.intregrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 30)
    
png(filename="CCA_1v2_before_alignment.Cond_Sorting.png", width=1500, height=700, res=150)
VlnPlot(object = Seurat_object.intregrated, features.plot = "CC1", group.by = "Cond_Sorting", 
        do.return = TRUE, x.lab.rot = TRUE)
dev.off()


png(filename="MetageneBicorPlot.Seurat_object.intregrated.png", width=3500, height=1000, res=150)
MetageneBicorPlot(Seurat_object.intregrated, grouping.var = "Cond_Sorting", dims.eval = 1:30)
dev.off()

png(filename="DimHeatmap.cca.Seurat_object.intregrated.png")
DimHeatmap(object = Seurat_object.intregrated, reduction.type = "cca", cells.use = 500, 
    dim.use = 1:12, do.balanced = TRUE)
dev.off()

saveRDS(Seurat_object.intregrated, paste(project_name, "Seurat_object.rds", sep="."))



q("no")
