# ---------------------------------------------------------------------------
# Step 04 - Differential expression, OIR vs normoxia, per cell type.
#
# ESTABLISHES
#     Per-cell-type gene-level OIR-vs-NORM contrasts at P14 and at P17, and
#     pooled across both, with the volcano plots that accompany them.
#
# WHY THIS WAY
#     The contrast is run within cell type rather than on all cells, because the
#     composition itself shifts with OIR (see step 02) and a whole-retina
#     contrast would confound a change in cell proportion with a change in
#     expression.
#
#     Cells are collapsed from 21 Cell_Type labels to 10 General_CellType labels
#     before testing - the map is in annotation/celltype_grouping.csv, which is
#     the same ten pairs used in steps 05 and 06 (checked: identical in all
#     three). Opticin_cells are dropped from this step.
#
# READ THIS BEFORE USING THE OUTPUT TABLES
#     Every table written by this script is filtered as
#     subset(DGE_test, p_val < 0.05) - the UNADJUSTED p-value - and the filename
#     nevertheless reads `p0.05`. The volcano plots in the same script colour
#     points on p_val_adj < 0.05, so the tables and the figures do not use the
#     same threshold. The adjusted p-value is present as a column in every
#     table, so nothing is lost, but `p0.05` in a filename here does not mean
#     FDR 0.05. Not reconciled, because this script produced published numbers;
#     recorded in METHODS_CHOICES.md.
#
#     One gene-set file it reads, Gmt.file/Pathways/GO_HEXOSE_CATABOLIC_PROCESS.gmt,
#     is EMPTY on the cluster (0 lines). The section that indexes
#     gsc_final[["GO_HEXOSE_CATABOLIC_PROCESS"]] therefore cannot have
#     contributed any genes. Left in place and flagged rather than removed.
#
# REQUIRES
#     step 02
#
# INPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#     <MSIGDB_DIR>/{h.all,c2.cp,c5.all}.v7.1.symbols.gmt
#     <PATHWAY_DIR>/GO_HEXOSE_CATABOLIC_PROCESS.gmt      EMPTY on the cluster
#     <CUSTOM_GENESET_DIR>/*.gmx                         UNAVAILABLE
#
# OUTPUTS
#     wholeEC_OIRvsNORM_{P14,P17}_WT_DGE_test_p0.05.txt   unadjusted-p filtered
#     wholeEC_OIRvsNORM_WT_DGE_test_p0.05.txt
#     per-cell-type DGE tables and volcano plots
#
# USAGE
#     Rscript 04_de_genes.R
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
OUT_DIR <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/DifferentialExpression/Genes")

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

# FindMarkers settings, as run
DE_LOGFC_MIN <- 0.5     # for the per-timepoint contrasts
DE_MIN_PCT   <- 0.2
DE_P_FILTER  <- 0.05    # applied to p_val, NOT p_val_adj - see the header

rm(list = ls())
library(Seurat)
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
library(Rmagic)
library(readr)
#library(phateR)
library(viridis)
library(magrittr)
library(harmony)
library(RColorBrewer)
library(tidyverse)
library(future)
library(GSVA)
library(GSEABase)
library(parallel)
#library(sigPathway)
library(limma)
#library(gskb)
library(Biobase)
library(genefilter)
library(GSVAdata)
data(c2BroadSets)
library(heatmap3)
plan("multicore", workers = (availableCores()-1))
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

Subseted_cells <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P14", "P17")))

Seurat_object_Subset <- subset(Seurat_object, cells = Subseted_cells)

## Renamed general cell type
Seurat_object_Subset <- RenameIdents(object = Seurat_object_Subset, 
                                      "Glycinergic_amacrine_cells" ='Amacrine_cells',
                                      'Bipolar_cells_1' = 'Bipolar_cells', 
                                      'Bipolar_cells_2' = 'Bipolar_cells', 
                                      'Early_rods' = 'Rods', 
                                      'Neuronal_progenitor_cells' = 'Amacrine_cells', 
                                      'Early_muller_glia' = 'Muller_glia', 
                                      'Early_bipolar_cells' = 'Bipolar_cells', 
                                      'Early_amacrine_cells' = 'Amacrine_cells', 
                                      #'Activated_muller_glia' = 'Muller_glia', 
                                      'Cholinergic_amacrine_cells' = 'Amacrine_cells')

Seurat_object_Subset[["General_CellType"]] <- Idents(object = Seurat_object_Subset)

Seurat_object_Subset[["General_CellType_Cond_TimePoint"]] <- paste(Seurat_object_Subset$General_CellType,Seurat_object_Subset$Cond_TimePoint, sep="_")


Seurat_object_Subset <- subset(Seurat_object_Subset, ident = "Opticin_cells", invert=T)


##Or remove unwanted cell types

Seurat_object_Subset <- SubsetData(Seurat_object_Subset, 
  ident.remove = c("Opticin_cells", "Neuronal_progenitor_cells", 'Early_amacrine_cells','Early_bipolar_cells', 'Early_muller_glia',
  "Early_rods"))

### Run Umap
Seurat_object_Subset <- RunUMAP(Seurat_object_Subset, dims = DIM_nb, reduction.name = "umap")

colorsplot <- DiscretePalette(21, palette = "polychrome")


setEPS(bg = "white", family = "Times", width=8, height=6)
postscript("umap.annotated.ident.OIR_P14_P17.eps")
DimPlot(Seurat_object_Subset, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "General_CellType", cols = rev(colorsplot))
dev.off()



setwd(OUT_DIR)


##Find DEGs in whole ECs

cells_in_this_cluster <- SetIdent(Seurat_object_Subset, value = "Condition")

  
#Perform DGE analysis using one of the model above

DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P14_WT", ident.2 = "NORM_P14_WT",logfc.threshold = 0.5, min.pct=0.2)
DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)
write.table(DGE_test_P14, paste("wholeEC_OIRvsNORM_P14_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


top.1 <- rownames(head(DGE_test_P14[order(DGE_test_P14$avg_logFC),], 5))
top.2 <- rownames(head(DGE_test_P14[order(-DGE_test_P14$avg_logFC),], 5))

top10 <- c(top.1,top.2)  

dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("wholeEC_top10_DE_genes.Wilcox.OIRvsNORM_P14.png", sep=""), dpi=300, width = 10, height=5)


DGE_test_P17 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P17_WT", ident.2 = "NORM_P17_WT",logfc.threshold = 0.5, min.pct=0.2)
DGE_test_P17 <- subset(DGE_test_P17, p_val < 0.05)
write.table(DGE_test_P17, paste("wholeEC_OIRvsNORM_P17_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


top.1 <- rownames(head(DGE_test_P17[order(DGE_test_P17$avg_logFC),], 5))
top.2 <- rownames(head(DGE_test_P17[order(-DGE_test_P17$avg_logFC),], 5))
top10 <- c(top.1,top.2)  
dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))


ggsave(dotplot, filename=paste("wholeEC_top10_DE_genes.Wilcox.OIRvsNORM_P17.png", sep=""), dpi=300, width = 10, height=5)


#Time point confounded

cells_in_this_cluster <- SetIdent(Seurat_object_Subset, value = "Condition")

DGE_test <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR", ident.2 = "NORM",logfc.threshold = -Inf, min.pct=-Inf)
DGE_test <- subset(DGE_test, p_val < 0.05)
write.table(DGE_test, paste("wholeEC_OIRvsNORM_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


top.1 <- rownames(head(DGE_test[order(DGE_test$avg_logFC),], 5))
top.2 <- rownames(head(DGE_test[order(-DGE_test$avg_logFC),], 5))

top10 <- c(top.1,top.2)  

dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("wholeEC_top10_DE_genes.Wilcox.OIRvsNORM.png", sep=""), dpi=300, width = 10, height=5)

#Make volcano plot


theme_Publication <- function(base_size=14, base_family="helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

volcanoplot <- ggplot(DGE_test,aes(avg_logFC,-log(p_val))) +
      geom_point(size=2) +
      geom_point(data = subset(DGE_test,(avg_logFC > 0.5) & (p_val < 0.05) & (p_val_adj < 0.05)),col="red") +
      geom_point(data = subset(DGE_test,(avg_logFC < -0.5 ) & (p_val < 0.05) & (p_val_adj < 0.05)),col="green") +
      geom_text_repel(
        data = subset(DGE_test, (avg_logFC > 0.5 | avg_logFC < -0.5) & (p_val < 0.05) & (p_val_adj < 0.05)),
        aes(label = rownames(subset(DGE_test, (avg_logFC > 0.5 | avg_logFC < -0.5) & (p_val < 0.05) & (p_val_adj < 0.05)))),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      theme_Publication() 
ggsave(volcanoplot, filename=paste("volcanoplot_wholeEC.Wilcox.OIRvsNORM.png", sep=""), dpi=300, width = 7, height=7)


##Find DEGs per subtype
#cell_type_list <- sort(unique(Seurat_object@ident))

#cell_type_list <- c("Tip ECs")                 

cell_type_list <- unique(Seurat_object@active.ident)

for(this_cluster in cell_type_list){

  #####this_cluster <- "Tip ECs"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Condition")


        if(length(grep("OIR_P14_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P14_WT", cells_in_this_cluster@meta.data$Condition))>2) {
  
        #Perform DGE analysis using one of the model above

        DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P14_WT", ident.2 = "NORM_P14_WT",logfc.threshold = 0.5, min.pct=0.2)
        DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)
        write.table(DGE_test_P14, paste(this_cluster,"OIRvsNORM_P14_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


        top.1 <- rownames(head(DGE_test_P14[order(DGE_test_P14$avg_logFC),], 5))
        top.2 <- rownames(head(DGE_test_P14[order(-DGE_test_P14$avg_logFC),], 5))

        top10 <- c(top.1,top.2)  

        dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))

        ggsave(dotplot, filename=paste("cluster_",this_cluster,"_top10_DE_genes.Wilcox.OIRvsNORM_P14.png", sep=""), dpi=300, width = 10, height=5)

        }

        if(length(grep("OIR_P17_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P17_WT", cells_in_this_cluster@meta.data$Condition))>2) {

        DGE_test_P17 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P17_WT", ident.2 = "NORM_P17_WT",logfc.threshold = 0.5, min.pct=0.2)
        DGE_test_P17 <- subset(DGE_test_P17, p_val < 0.05)
        write.table(DGE_test_P17, paste(this_cluster,"OIRvsNORM_P17_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


        top.1 <- rownames(head(DGE_test_P17[order(DGE_test_P17$avg_logFC),], 5))
        top.2 <- rownames(head(DGE_test_P17[order(-DGE_test_P17$avg_logFC),], 5))

        top10 <- c(top.1,top.2)  

        dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))


        ggsave(dotplot, filename=paste("cluster_",this_cluster,"_top10_DE_genes.Wilcox.OIRvsNORM_P17.png", sep=""), dpi=300, width = 10, height=5)

        }




    ## Make venn diag for all differentially expressed genes containing testing results

        library(VennDiagram)
        flog.threshold(ERROR)

        if(length(grep("OIR_P14_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P14_WT", cells_in_this_cluster@meta.data$Condition))>2
          & length(grep("OIR_P17_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P17_WT", cells_in_this_cluster@meta.data$Condition))>2) {

        venn.diagram(x = list(P14=rownames(DGE_test_P14), P17=rownames(DGE_test_P17)), filename= paste("cluster_",this_cluster,"_VennDiagram_OIRvsNORM.png"), 
        height = 2000, width = 2000, resolution =500, 
        imagetype = "png", units = "px", compression ="lzw", 
        na = "stop", main = NULL, sub = NULL, main.pos= c(0.5, 1.05), main.fontface = "plain",
        main.fontfamily = "serif", main.col = "black",
        main.cex = 2, main.just = c(0.5, 1), sub.pos = c(0.5,
        1.05), sub.fontface = "plain", sub.fontfamily ="serif", sub.col = "blue", sub.cex = 1)
        }
}


##########
##########

##Heatmap number of DEGs per subtype per timepoint



#cell_type_list <- c("Muller_glia", "Astrocytes")                 

cell_type_list <- unique(Seurat_object_Subset@active.ident)

timepoint_list <- unique(Seurat_object_Subset[["TimePoint"]]$TimePoint)


names <- cell_type_list
finaltable <- data.frame()
for (k in names) finaltable[[k]] <- as.integer()


for(this_timepoint in timepoint_list){

## Print status for which identity is being processed
  print(paste("Working on timepoint #",this_timepoint,sep=""))


  cells_in_this_time_point <- SetIdent(Seurat_object_Subset, value = "TimePoint")

  cells_in_this_time_point <- SubsetData(cells_in_this_time_point,
                                      ident.use=this_timepoint)


  intertable <- data.frame()
  intertable <- intertable[1:length(this_timepoint),]
  rownames(intertable) <- this_timepoint

  for(this_cluster in cell_type_list){

    

    #####this_cluster <- "Tip ECs"

    ## Print status for which identity is being processed
    print(paste("Working on cluster #",this_cluster,sep=""))

    ## Subset Seurat object to only contain cells from this cluster

    cells_in_this_cluster <- SetIdent(cells_in_this_time_point, value = "Cell_Group")

    cells_in_this_cluster <- SubsetData(cells_in_this_cluster,
                                        ident.use=this_cluster)

    ## Check whether there are cells in both groups, otherwise skip this cluster
          
    cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Condition")


          if(length(grep("OIR", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM", cells_in_this_cluster@meta.data$Condition))>2) {
    
          #Perform DGE analysis using one of the model above

          DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR", ident.2 = "NORM", logfc.threshold = 0.5, min.pct=0.1)
          DGE_test_P14 <- subset(DGE_test_P14, p_val_adj < 0.05)
          #write.table(DGE_test_P14, paste(this_cluster,"OIRvsNORM_P14_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)
          table <- data.frame(length(rownames(DGE_test_P14)))
          rownames(table) <- this_timepoint
          colnames(table) <- this_cluster

          intertable <- merge(intertable,table, by="row.names", all.x=TRUE)

          rownames(intertable) <- intertable$Row.names
   
          intertable$Row.names <- NULL
          }
      }

  finaltable <- dplyr::bind_rows(finaltable, intertable)

}

finaltable1 <- finaltable

finaltable <- t(finaltable1)

finaltable <- finaltable[, c("P12", "P14", "P17")]

my_palette <- colorRampPalette(c("green", "yellow", "red"))


setEPS(bg = "white", family = "Times", width=10, height=15)
postscript(paste("heatmap_","_OIR_NORM_P12_P14.P17.DEGs.eps",sep=""))
heatmap.2(as.matrix(finaltable), col=my_palette(50), symbreak=FALSE, trace='none', cexRow=2, cexCol= 2,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=FALSE,
          srtCol=90,
          revC=FALSE, key.title = NA, density.info="none", 
          key.ylab=NA, 
          main="P17", 
          margins=c(25,25),
          key.par=list(mar=c(5, 5, 10, 5))
          )
dev.off()





###Check expression of specific genes

png(paste("FeaturePlot_ANGPT2_ESM1.png", sep=""), res=300, width = 4000, height=1000)
FeaturePlot(
  Seurat_object_Subset,
  c("ESM1", "ANGPT2"),
  dims = c(1, 2),
  cells = NULL,
  cols = c("red", "green"),
  pt.size = NULL,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  shape.by = NULL,
  slot = "data",
  blend = TRUE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = FALSE,
  combine = TRUE
)
dev.off()



pdf(paste("DotPlot_REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS.bycondition.pdf", sep=""), width = 10, height=10)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("ACADL",
"ACADM",
"ACADS",
"ACADVL",
"ACSM3",
"ACSM6",
"ECHS1",
"HADH",
"HADHA",
"HADHB",
"MECR", "DECR1", "ECI1", "CPT1A", "SIRT3"),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "General_CellType_Cond_TimePoint",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()


png(paste("DotPlot_tuft_markers.bycondition.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("AQP1",
"PLVAP",
"HADHA", "HBEGF","TFPI"),
  cols = c("darkblue", "orange", "purple", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = "Condition",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
dev.off()

png(paste("DotPlot_TREM2.bycondition.png", sep=""), res=300, width = 3000, height=4000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("TREM2"),
  cols = c("darkblue", "orange", "purple", "red", "green", "pink"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = "Cond_TimePoint",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
dev.off()


png(paste("DotPlot_Tipcells_glycolysis_markers_OIRvsNORM_P14_P17.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("PKM",
"ALDOA",
"ENO1", "TPI1","PFKP", "LDHA", "CHST1", "CXCR4", "MIF", "PLOD1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

png(paste("DotPlot_Tipcells_glyco_ppp_hexo_markers_OIRvsNORM_P14_P17.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("FBP1", "EXT1", "PFKP", "IDH1", "EGLN3", "SLC16A3", "PGK1", "HK1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

png(paste("DotPlot_Tipcells_glyco_ppp_hexo_markers_OIRvsNORM_P14_P17.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("FBP1", "EXT1", "PFKP", "IDH1", "EGLN3", "SLC16A3", "PGK1", "HK1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

###Get genes from gene list

geneset1 <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "MALLETTE_SENESCENCE_UP.gmx"), verbose = TRUE)

geneset1 <- GeneSet(geneset1[[1]]$probes, geneIdType=SymbolIdentifier(), setName="MALLETTE_SENESCENCE_UP")

geneset2 <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "GLOBAL_SENESCENCE_LITERATURE_CURATED_2020.gmx"), verbose = TRUE)

geneset2 <- GeneSet(geneset2[[1]]$probes, geneIdType=SymbolIdentifier(), setName="GLOBAL_SENESCENCE_LITERATURE_CURATED")

geneset3 <- importGeneSets(file.path(CUSTOM_GENESET_DIR, "SAWCHYN_UNBIASED_UP_v4.gmx"), verbose = TRUE)

geneset3 <- GeneSet(geneset3[[1]]$probes, geneIdType=SymbolIdentifier(), setName="SAWCHYN_UNBIASED_UP")

geneset4 <- getGmt(file.path(PATHWAY_DIR, "GO_HEXOSE_CATABOLIC_PROCESS.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.h <- getGmt(file.path(MSIGDB_DIR, "h.all.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.c2 <- getGmt(file.path(MSIGDB_DIR, "c2.cp.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gmt.MsigDB.c5 <- getGmt(file.path(MSIGDB_DIR, "c5.all.v7.1.symbols.gmt"), geneIdType=SymbolIdentifier())

gsc_final <- GeneSetCollection(c(gmt.MsigDB.h, gmt.MsigDB.c2,gmt.MsigDB.c5,geneset1,geneset2,geneset3,geneset4))

###

geneset <- names(gsc_final)[unique(c(
	grep("REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS$", names(gsc_final))
	))]

geneset <- sort(unique(geneIds(gsc_final[[geneset]])))

Seurat_object_Subset <- ReorderIdent(Seurat_object_Subset, unique(geneIds(gsc_final[[geneset]])))

dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = geneset,
  cols = c("gold1", "red2"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = "TimePoint",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS.OIR_P14_P17.png", sep=""), 
  dpi=150, width = 10, height=6)


library("Nebulosa")
png("plot_density_FAO.png", width = 4000, height =1500, res = 150)
plot_density(Seurat_object_Subset, 
  c("HADHA", "ACADL"), 
  joint = TRUE, combine = TRUE, reduction = "umap", method = c("wkde"), adjust = 2) + 
plot_layout(ncol = 3)
dev.off()

##SUbset conditions

cond4subset <- c("NORM_P14")

Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Cond_TimePoint %in% cond4subset))

Seurat_object_cond <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)

Seurat_object_cond <- ReorderIdent(Seurat_object_cond, c("HMGCS2", "HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1"))

dotplot <- DotPlot(
  Seurat_object_cond,
  assay = NULL,
  features = c("HMGCS2", "HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1"),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_ketonegenes.",cond4subset,".png", sep=""), 
  dpi=150, width = 8, height=4)

ggsave(dotplot, filename=paste("DotPlot_ketonegenes.",cond4subset,".eps", sep=""), 
  dpi=300, width = 8, height=4)


Seurat_object_cond <- RunUMAP(Seurat_object_cond, dims = DIM_nb, reduction.name = "umap")

colorsplot <- DiscretePalette(21, palette = "polychrome")


setEPS(bg = "white", family = "Times", width=8, height=6)
postscript("umap.annotated.ident.OIR_P14.eps")
DimPlot(Seurat_object_cond, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "General_CellType", cols = rev(colorsplot))
dev.off()


setEPS(bg = "white", family = "Times", width=12, height=6)
postscript("plot_density_FAO.P14_OIR.eps")
plot_density(Seurat_object_cond, 
  c("HADHA", "ACADL"), 
  joint = TRUE, combine = TRUE, reduction = "umap", method = c("wkde"), adjust = 1) + 
plot_layout(ncol = 3)
dev.off()


####    Make heatmap of gene expression

cond4subset <- c("NORM_P14", "NORM_P17", "OIR_P14","OIR_P17")

genes <- c("SIRT3")

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(unique(Seurat_object_Subset@active.ident)),]

rownames(dataframe.all) <- unique(Seurat_object_Subset@active.ident)


#Start for loop

for(thiscond in cond4subset) {

    print(paste("Working on", thiscond, sep=" "))

    Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Cond_TimePoint %in% thiscond))

    Seurat_object_cond <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)

    AverageExpression.data <- AverageExpression(
      Seurat_object_cond,
      assays = "RNA",
      features = genes,
      return.seurat = FALSE,
      add.ident = NULL,
      slot = "scale.data",
      use.scale = FALSE,
      use.counts = FALSE,
      verbose = TRUE)

    df <- t(AverageExpression.data$RNA)

    colnames(df) <- paste(colnames(df), thiscond, sep="_")

    dataframe.all <- merge(dataframe.all,df, by=0, all=TRUE)

    rownames(dataframe.all) <- dataframe.all$Row.names
       
    dataframe.all$Row.names <- NULL 
}

dataframe <- dataframe.all[-9,]

#dataframe = scale(t(dataframe))

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))



setEPS(bg = "white", family = "Times", width=10, height=15)
postscript(paste("heatmap_",genes,"_OIR_NORM_P14.P17.scale.data.eps",sep=""))
heatmap.2(as.matrix(dataframe), col=my_palette(50), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=FALSE,
          srtCol=90,
          revC=FALSE, key.title = NA, density.info="none", 
          key.ylab=NA, 
          main="P17", 
          margins=c(25,15),
          key.par=list(mar=c(5, 5, 10, 5))
          )
dev.off()


####    Make heatmap of gene expression ratio between OIR and NORM

#genelist <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7", "NAMPT", "HADHA", "NT5C3")

genelist <- c("SIRT3", "CPT1A", "HADHA")

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(genelist),]

rownames(dataframe.all) <- genelist


#Start for loop

timepoints <- c("P14")

#thistimepoint <- "P14"

celltypelist <- unique(Seurat_object_Subset@active.ident)

#thiscelltype <- "Endothelial_cells"


for(thistimepoint in timepoints) {

  print(paste("Working on", thistimepoint, sep=" "))

  Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, TimePoint %in% thistimepoint))

  Seurat_object_timepoints <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)


        for(thiscelltype in celltypelist) {

            print(paste("Working on", thiscelltype, sep=" "))

            Subseted_cells <- rownames(subset(Seurat_object_timepoints@meta.data, General_CellType %in% thiscelltype))

            Seurat_object_celltype <- SubsetData(Seurat_object_timepoints, cells = Subseted_cells)

            Seurat_object_celltype <- SetIdent(Seurat_object_celltype, value = "Condition")

            DGE_test <- FindMarkers(Seurat_object_celltype, ident.1 = "OIR", ident.2 = "NORM", logfc.threshold = -Inf, min.pct=-Inf, features = genelist, only.pos=FALSE)


            df <- as.data.frame(DGE_test[,"avg_logFC"], row.names=rownames(DGE_test))

            colnames(df) <- paste(thiscelltype, thistimepoint, "OIRvsNORM", sep="_")

            dataframe.all <- merge(dataframe.all,df, by=0, all=TRUE)

            rownames(dataframe.all) <- dataframe.all$Row.names
               
            dataframe.all$Row.names <- NULL 
        }
    }


dataframe <- t(dataframe.all)

#dataframe = scale(t(dataframe))

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))



setEPS(bg = "white", family = "Times", width=10, height=15)
postscript(paste("heatmap_","Sirt3_Cpt1a_Hadha","_OIR_NORM_P14.OIRvsNORM_ratio.eps",sep=""))
heatmap.2(as.matrix(dataframe), col=my_palette(50), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=FALSE,
          srtCol=90,
          revC=FALSE, key.title = NA, density.info="none", 
          key.ylab=NA, 
          main="P17", 
          margins=c(25,25),
          key.par=list(mar=c(5, 5, 10, 5))
          )
dev.off()


####    Make heatmap of gene expression ratio between P14 and P17

genelist <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7", "NAMPT", "HADHA", "NT5C3")

genelist <- c("SIRT3", "CPT1A", "HADHA")


dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(genelist),]

rownames(dataframe.all) <- genelist


#Start for loop

Conditions <- c("OIR")

#thistimepoint <- "P14"

celltypelist <- unique(Seurat_object_Subset@active.ident)

#thiscelltype <- "Endothelial_cells"


for(thisCondition in Conditions) {

  print(paste("Working on", thisCondition, sep=" "))

  Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Condition %in% thisCondition))

  Seurat_object_condition <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)


        for(thiscelltype in celltypelist) {

            print(paste("Working on", thiscelltype, sep=" "))

            Subseted_cells <- rownames(subset(Seurat_object_condition@meta.data, General_CellType %in% thiscelltype))

            Seurat_object_celltype <- SubsetData(Seurat_object_condition, cells = Subseted_cells)

            Seurat_object_celltype <- SetIdent(Seurat_object_celltype, value = "TimePoint")

            DGE_test <- FindMarkers(Seurat_object_celltype, ident.1 = "P17", ident.2 = "P14", logfc.threshold = -Inf, min.pct=-Inf, features = genelist, only.pos=FALSE)


            df <- as.data.frame(DGE_test[,"avg_logFC"], row.names=rownames(DGE_test))

            colnames(df) <- paste(thiscelltype, thisCondition, "P17vsP14", sep="_")

            dataframe.all <- merge(dataframe.all,df, by=0, all=TRUE)

            rownames(dataframe.all) <- dataframe.all$Row.names
               
            dataframe.all$Row.names <- NULL 
        }
    }

dataframe <- dataframe.all[c(6,7,8),-7]


dataframe <- t(dataframe)

#dataframe = scale(t(dataframe))

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

x <- as.matrix(dataframe)

setEPS(bg = "white", family = "Times", width=10, height=15)
postscript(paste("heatmap_","Sirt3_Cpt1a_Hadha","_OIR_P17vsP14_ratio.eps",sep=""))
heatmap.2(x, col=my_palette(50), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=FALSE,
          srtCol=90,
          revC=FALSE, key.title = NA, density.info="none", 
          key.ylab=NA, 
          main=NULL, 
          margins=c(25,25),
          key.par=list(mar=c(5, 5, 10, 5)),
          distfun = function(x) dist(x, method = 'canberra')
          )
dev.off()


###FInd markers using gene list

geneset <- names(gsc_final)[unique(c(
	grep("HALLMARK_GLYCOLYSIS", names(gsc_final)),
	grep("KEGG_PENTOSE_PHOSPHATE_PATHWAY", names(gsc_final)),
	grep("REACTOME_PENTOSE_PHOSPHATE_PATHWAY", names(gsc_final)),
	grep("GO_PENTOSE_PHOSPHATE_SHUNT_NON_OXIDATIVE_BRANCH", names(gsc_final)),
	grep("GO_HEXOSE_CATABOLIC_PROCESS", names(gsc_final))
	))]


genelist <- unique(c(geneIds(gsc_final[["HALLMARK_GLYCOLYSIS"]]),
		 geneIds(gsc_final[["KEGG_PENTOSE_PHOSPHATE_PATHWAY"]]),
		 geneIds(gsc_final[["REACTOME_PENTOSE_PHOSPHATE_PATHWAY"]]),
		 geneIds(gsc_final[["GO_PENTOSE_PHOSPHATE_SHUNT_NON_OXIDATIVE_BRANCH"]]),
		 geneIds(gsc_final[["GO_HEXOSE_CATABOLIC_PROCESS"]])
	   ))


genelist <- intersect(genelist,rownames(Seurat_object_Subset))


###Subset condition of interest

cond4subset <- c("OIR_P14_WT")

Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Condition %in% c(cond4subset)))

Seurat_object_Subset_Condition <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)

DGE_test <- FindMarkers(Seurat_object_Subset_Condition, ident.1 = "Tip ECs_OIR_P14_WT", logfc.threshold = 0.25, min.pct=0.1, features = genelist, only.pos=TRUE)

DGE_test <- subset(DGE_test, p_val < 0.05)


dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = rownames(DGE_test),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_Tip_P14_Markers_GLYCO_PPP_HEXOSE.",cond4subset,".png", sep=""), dpi=150, width = 20, height=6)


dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = rownames(DGE_test),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_Tip_P14_Markers_GLYCO_PPP_HEXOSE.byTreatment.png", sep=""), dpi=150, width = 20, height=8)



q("no")
