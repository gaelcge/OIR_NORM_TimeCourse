# ---------------------------------------------------------------------------
# Step 05 - Differential gene-set activity, OIR vs normoxia, per cell type.
#
# ESTABLISHES
#     Which gene-set programmes - senescence among them - differ between OIR and
#     normoxia within each cell type, by testing the GSVA scores from step 03
#     the same way step 04 tests genes.
#
# WHY THIS WAY
#     The GSVA score matrix is loaded back onto the Seurat object as an assay
#     named "GO" and then CLR-normalised. Putting the scores in an assay is what
#     lets FindMarkers, the same routine used for genes in step 04, be used for
#     gene sets, so the two analyses are comparable rather than two different
#     tests. CLR is the normalisation Seurat uses for feature sets that are not
#     counts.
#
# REQUIRES
#     step 03 (reads its GSVA csv - this is a hard dependency, and the reason
#              GSVA is numbered before this step rather than after it)
#     step 02 (reads its annotated object)
#
#     Because step 03 is not reproducible without the three unavailable curated
#     gene sets, this step is not reproducible either. See
#     reference/PROVENANCE.txt.
#
# INPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#     <GSVA_DIR>/gsva.exprs.MsigDB_h_c2_c5.senescence.all.v7.1.symbols.csv
#     <GENE_LIST_DIR>/SIRT3_TARGET_GENES.txt                 58 lines
#     <GENE_LIST_DIR>/WP_NAD_METABOLISM_SIRTUINS_AND_AGING.txt 10 lines
#
# OUTPUTS
#     per-cell-type gene-set DGE tables, Nebulosa density plots, and the
#     clustered heatmaps of gene-set activity
#
# USAGE
#     Rscript 05_de_functions.R
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
OUT_DIR <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/DifferentialExpression/Functions")

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

GO_ASSAY_NORM <- "CLR"   # gene-set scores are not counts
DE_MIN_PCT    <- 0.1

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
library(magrittr)
library(harmony)
library(RColorBrewer)
library(tidyverse)
library(future)
plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)


#Set directory of dataset to analyse
setwd(OBJ_DIR)

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"

res = 1

DIM_nb <- c(1:20)

Dim <- 20

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))


##SUbset cells of interest

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Cond_TimePoint %in% c("NORM_P14", "NORM_P17", "OIR_P14", "OIR_P17")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells = Subseted_cells)

Seurat_object_Subset <- RenameIdents(object = Seurat_object_Subset, 
                                      "Glycinergic_amacrine_cells" ='Amacrine_cells',
                                      'Bipolar_cells_1' = 'Bipolar_cells', 
                                      'Bipolar_cells_2' = 'Bipolar_cells', 
                                      'Early_rods' = 'Rods', 
                                      'Neuronal_progenitor_cells' = 'Amacrine_cells', 
                                      'Early_muller_glia' = 'Muller_glia', 
                                      'Early_bipolar_cells' = 'Bipolar_cells', 
                                      'Early_amacrine_cells' = 'Amacrine_cells', 
                                      'Activated_muller_glia' = 'Muller_glia', 
                                      'Cholinergic_amacrine_cells' = 'Amacrine_cells')

Seurat_object_Subset[["General_CellType"]] <- Idents(object = Seurat_object_Subset)

Seurat_object_Subset <- SubsetData(Seurat_object_Subset, ident.remove = "Opticin_cells")


Seurat_object_Subset <- RunUMAP(Seurat_object_Subset, dims = DIM_nb, reduction.name = "umap")

colorsplot <- DiscretePalette(21, palette = "polychrome")

pdf("umap.annotated.ident.P14_P17.splitbyCondition.pdf", width=15, height=8, bg = "white")
DimPlot(Seurat_object_Subset, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "ident", split.by="Condition", cols = rev(colorsplot))
dev.off()


setEPS(bg = "white", family = "Times", width=15, height=8)
postscript("umap.annotated.ident.P14_P17.splitbyCondition.eps")
DimPlot(Seurat_object_Subset, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "ident", split.by="Condition", cols = rev(colorsplot))
dev.off()

png(filename="umap.annotated.Cell_Type.P14_P17.int.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object_Subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "Cell_Type") + NoLegend()
dev.off()


###R 4.0
library("Nebulosa")
png("plot_density_Tuft.png", width = 4000, height =2000, res = 150)
plot_density(Seurat_object, 
  c("MAL", "NKD1", "PLVAP", "SAT1", "COL18A1", "EML1", "COL15A1", "PRSS23", "SIX3", "SIX3OS1",  "AQP1"), 
  joint = TRUE, combine = TRUE, reduction = "umap", method = c("wkde"), adjust = 2) + 
plot_layout(ncol = 4)
dev.off()

#### For Diff functions between conditions globally

##Analysis global
setwd(OUT_DIR) 

###Add GO function to Seurat object

GO_score= read.table(file.path(GSVA_DIR, "gsva.exprs.MsigDB_h_c2_c5.senescence.all.v7.1.symbols.csv"), sep = ",", header = T, row.names=1, stringsAsFactors=F)

Seurat_object_Subset[["GO"]] <- CreateAssayObject(counts = GO_score)

Seurat_object_Subset <- NormalizeData(Seurat_object_Subset, assay = "GO", normalization.method = "CLR")

Seurat_object_Subset <- ScaleData(Seurat_object_Subset, assay = "GO", verbose = FALSE, model.use = "negbinom")

Seurat_object_Subset@meta.data$Cell_SubType <- factor(Seurat_object_Subset@meta.data$Cell_SubType, 
                                                  levels = c("Arterial ECs", "Capillary II", "Capillary I", "Tuft", "Proliferative Stalk ECs", "Vein ECs", "Tip ECs")
                                                  )

### Make ridge plot for specific functions acrross cell types


rowf <- grep("SIRT", rownames(Seurat_object_Subset[["GO"]]))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

for (thisfunction in selected_functions) {

    Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Cond_TimePoint")

    loop_list <- unique(Seurat_object_loop@active.ident)

    for(this_cluster in loop_list){

        cells_in_this_cluster <- SubsetData(Seurat_object_loop,
                                          ident.use=this_cluster)

        ridgeplot <- RidgePlot(
        cells_in_this_cluster,
        assay = "GO",
        features=thisfunction,
        cols = viridis(11),
        idents = NULL,
        sort = "decreasing",
        group.by = "General_CellType",
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "scale.data",
        combine = TRUE) + NoLegend()
        
        ggsave(ridgeplot, filename=paste(this_cluster,thisfunction,"eps",sep="."), device = "eps", dpi=300, width = 5, height=4)

    }
}


###Make heatmap for specific functions acrross cell types

rowf <- unique(c(grep("GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
          grep("GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]]))
          ))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

#selected_functions <- rownames(Seurat_object_Subset[["GO"]])

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(selected_functions),]

rownames(dataframe.all) <- selected_functions

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Age")

timepoint <- "P17"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=timepoint)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Dataset")

Dataset <- "OIR_WT"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=Dataset)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Cell_SubType")

for(this_cluster in unique(Seurat_object_loop@active.ident)) {

  if(length(grep(this_cluster, Seurat_object_loop@active.ident,fixed = TRUE)) > 3) {

    print(paste("Working on cluster #",this_cluster,sep=""))

    DGE_test_P14 <- FindMarkers(Seurat_object_loop, ident.1=this_cluster, logfc.threshold = -Inf, min.pct=0.1, 
          features=selected_functions, assay = "GO", only.pos=FALSE)

    #DGE_test_P14 <- subset(DGE_test_P14, p_val_adj < 0.05)

    #DGE_test_P14 <- subset(DGE_test_P14, avg_logFC > 0.4 | avg_logFC < -0.4 )

    DGE_test_P14_avg_logFC <- as.data.frame(DGE_test_P14$avg_logFC)

    rownames(DGE_test_P14_avg_logFC) <- rownames(DGE_test_P14)

    colnames(DGE_test_P14_avg_logFC) <- paste(this_cluster, "avg_logFC", timepoint, sep="_")

    dataframe.all <- merge(dataframe.all,DGE_test_P14_avg_logFC, by=0, all=TRUE)

    rownames(dataframe.all) <- dataframe.all$Row.names
   
    dataframe.all$Row.names <- NULL
 
  }
   
}

dataframe <- dataframe.all

##Select for specific subtypes
dataframe <- dataframe[complete.cases(dataframe[ ,7]),]

#Remove row containing only NA
ind <- apply(dataframe, 1, function(x) all(is.na(x)))

dataframe <- dataframe[ !ind, ]

#### make matrix for heatmap

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

dataframe <- as.matrix(dataframe)


Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
           ladderize %>% rotate_DendSer

        #png("Dend.png")
        #plot(Colv)
        #dev.off()


Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", 2) %>% set("branches_lwd", 2) %>%
           ladderize


dataframe[is.na(dataframe)] <- 0
dataframe[is.nan(dataframe)] <- 0

dim(dataframe)

png(filename=paste("EC_markers_",timepoint,Dataset,"_GLYCOLYSIS_ECs_heatmap.png", sep=""), width=4500, height=4500, res = 300)  
heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=TRUE,
          srtCol=45,
          revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,5), 
          key.par=list(mar=c(5, 8, 8, 8)), 
          margins=c(65,45), key.ylab=NA, main="OIR_WT P17")
dev.off()



###Make heatmap of OIR vs NORM for specific functions acrross cell types

rowf <- unique(c(grep("BETA-OX", rownames(Seurat_object_Subset[["GO"]])),
          grep("BETA-OX", rownames(Seurat_object_Subset[["GO"]]))
          ))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(selected_functions),]

rownames(dataframe.all) <- selected_functions

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Age")

timepoint <- "P17"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=timepoint)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Cell_SubType")                                      

loop_list <- unique(Seurat_object_loop@active.ident)

loop_list <- c("Arterial ECs", "Capillary II", "Capillary I", "Vein ECs", "Tip ECs")


#loop_list <- loop_list[-5]

for(this_cluster in loop_list){

  #####this_cluster <- "Tip_ECs"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_loop,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Treatment")

  DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR", ident.2 = "NORM", logfc.threshold = -Inf, min.pct=-Inf, 
    features=selected_functions, assay = "GO")

  #DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)

  DGE_test_P14_avg_logFC <- as.data.frame(DGE_test_P14$avg_logFC)

  rownames(DGE_test_P14_avg_logFC) <- rownames(DGE_test_P14)

  colnames(DGE_test_P14_avg_logFC) <- paste(this_cluster, "avg_logFC_OIR_vs_NORM", timepoint, sep="_")

  dataframe.all <- merge(dataframe.all,DGE_test_P14_avg_logFC, by=0, all=TRUE)

  rownames(dataframe.all) <- dataframe.all$Row.names
 
  dataframe.all$Row.names <- NULL
}

dataframe <- dataframe.all

##Select for specific subtypes
#dataframe <- dataframe[complete.cases(dataframe[ ,5]),]

#Remove row containing only NA
#ind <- apply(dataframe, 1, function(x) all(is.na(x)))

#dataframe <- dataframe[ !ind, ]

#### make matrix for heatmap

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

dataframe <- as.matrix(dataframe)


Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
           ladderize %>% rotate_DendSer

        #png("Dend.png")
        #plot(Colv)
        #dev.off()


Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", 2) %>% set("branches_lwd", 2) %>%
           ladderize


dataframe[is.na(dataframe)] <- 0
dataframe[is.nan(dataframe)] <- 0

png(filename=paste("OIR_vs_NORM_",timepoint,"_BETA_OX_ECs_heatmap.png", sep=""), width=3500, height=3000, res = 200)  
heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=TRUE,
          srtCol=45,
          revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,5), 
          key.par=list(mar=c(5, 8, 8, 8)), 
          margins=c(25,55), key.ylab=NA, main="OIR vs NORM P17")
dev.off()


### Add module score


module_features <- read.table(file.path(GENE_LIST_DIR, "WP_NAD_METABOLISM_SIRTUINS_AND_AGING.txt"), stringsAsFactors=F)

Seurat_object_Subset <- AddModuleScore(
  object = Seurat_object_Subset,
  features = list(module_features$V1),
  ctrl = 5,
  name = 'WP_NAD_METABOLISM_SIRTUINS_AND_AGING'
)

module_features <- read.table(file.path(GENE_LIST_DIR, "SIRT3_TARGET_GENES.txt"), stringsAsFactors=F)

Seurat_object_Subset <- AddModuleScore(
  object = Seurat_object_Subset,
  features = list(module_features$V1),
  ctrl = 5,
  name = 'SIRT3_TARGET_GENES'
)


cells_in_this_cluster <- SetIdent(Seurat_object_Subset, value = "Cond_TimePoint")

thisfunction <- 'SIRT3_TARGET_GENES1'
timepoint <- "OIR_P14"

ridgeplot <- RidgePlot(
        cells_in_this_cluster,
        features= thisfunction,
        cols = viridis(11),
        idents = NULL,
        sort = "decreasing",
        group.by = "General_CellType",
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "data",
        combine = TRUE) + NoLegend()
        
ggsave(ridgeplot, filename=paste(thisfunction,timepoint,"P17.eps",sep="."), device = "eps", dpi=300, width = 5, height=7)


q("no")
