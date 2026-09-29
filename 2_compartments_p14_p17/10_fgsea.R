# ---------------------------------------------------------------------------
# Step 10 - fgsea on ranked OIR-vs-normoxia gene lists.
#
# ESTABLISHES
#     Pathway enrichment per cell type from the step 06 rankings.
#
# WHY THIS WAY
#     Rankings are taken unfiltered, as in the sibling analysis: fgsea needs a
#     complete ranked list and pre-filtering by fold change truncates it.
#
# WHICH OBJECT THIS RUNS ON
#     RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30 - P14 and P17
#     only, integrated by Cond_Sorting, resolution 3, 17 dims, built under
#     Seurat v2 / R 3.5.0. This is NOT the P5-P17 object used by
#     1_timecourse_p5_p17; the two are different integrations of the same
#     merged matrix and their cluster numbers are unrelated.
#
# REQUIRES
#     step 03, and step 06 for the contrast definition
#
# INPUTS
#     the annotated object and the MSigDB v6.2 collections
#
# OUTPUTS
#     fgsea tables and NES figures under <B_DIR>/fgsea/
#
# USAGE
#     Rscript 10_*.R
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
library(devtools)
library(fgsea)
library(useful)
library(plotly)
library(ggrepel)
library(gtools)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gplots)
library(Matrix)
library(Seurat)
library(GSEABase)
library(GSVAdata)
library(GSVA) 
data(c2BroadSets)
library(tidyverse)
library(R.utils)



#Select dataset
setwd(file.path(B_DIR, "CellTypeAnnotation/"))

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting"

res = 3

DIM_nb <- c(1:17)

Dim <- 17

perp = 30

Seurat_object <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.annotated.rds", sep="."))

###Re edit object ident

Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Quiescent Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Activated Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Mki67 Muller Glia"), new.ident.name = "Muller glia")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Early Muller Glia"), new.ident.name = "Muller glia")

##OPtional
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Muller glia"), new.ident.name = "Neuroglia cells")
Seurat_object <- RenameIdent(Seurat_object, old.ident.name = c("Astrocytes"), new.ident.name = "Neuroglia cells")

Seurat_object <- SubsetData(Seurat_object, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = c("Basal cells 1", "Basal cells 2", "Red blood cells"), accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)

## subet P17

Subseted_cells_P17 <- rownames(subset(Seurat_object@meta.data, TimePoint %in% c("P17")))

Seurat_object.P17 <- SubsetData(Seurat_object, cells.use = Subseted_cells_P17, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = Inf,
  random.seed = 1)



#Select pathway

GOI <- "SENESCENCE_SASP_Sawchyn_UP"

dir.create(paste(file.path(B_DIR, "fgsea/"), GOI, sep=""))

setwd(paste(file.path(B_DIR, "fgsea/"), GOI, sep=""))


broadset_selected <- "sup_h.c2.c5.c6.c7.all.270519.symbols"

broadSet.custom <- getGmt(paste(file.path(GMT_DIR, ""), broadset_selected, ".gmt", sep=""),
              geneIdType=SymbolIdentifier())

GO_term <- geneIds(broadSet.custom[[GOI]])

#or via gene set list

#GO_term <- fread(paste(file.path(GENE_LIST_DIR, ""), GOI, ".csv", sep=""), header=F, fill=TRUE)

#GO_term <- GO_term$V1

# Prepare the gene list and run GSEA


GO_term <- list(intersect(rownames(Seurat_object@data),GO_term))

celltypes <- names(summary(Seurat_object@ident))

#celltypes <- "Horizontal cells"

###Start the loop for fgsea on each cell type

##Subset NORM vs OIR


Conditions <- c("NORM", "OIR")

for (Condition_subset in Conditions) {

  Subseted_cells <- rownames(subset(Seurat_object.P17@meta.data, Condition %in% c(Condition_subset)))

  retina_subset <- SubsetData(Seurat_object.P17, cells.use = Subseted_cells, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = 50,
  random.seed = 1)


  df_total = NULL
  for (celltypename in celltypes) {
    retina_CTL_celltype <- SubsetData(retina_subset, cells.use = NULL, subset.name =NULL, ident.use = celltypename,
      ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
      do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
      random.seed = 1)

    print(summary(retina_CTL_celltype@ident))

    retina_CTL_celltype.scaled <- as.data.frame(retina_CTL_celltype@scale.data)

    CellBarCode_name <- unique(colnames(retina_CTL_celltype.scaled))
    
    gsea=NULL
    
    for (sel_CellBarCode in CellBarCode_name) {

    #CellBarCode_cell <- rownames(subset(retina_CTL_celltype@meta.data, CellBarCode == sel_CellBarCode))

    #retina_CellBarCode_CTL <- SubsetData(retina_CTL_celltype, cells.use = CellBarCode_cell, subset.name =NULL, ident.use = NULL,
    #ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
    #do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
    #random.seed = 1)

    #Plot enrichment graph for the cell type (average)

      if (length(rownames(retina_CTL_celltype@meta.data)) >1) {

      #average_cell_Type <- AverageExpression(retina_CTL_celltype, genes.use = NULL, return.seurat = FALSE,
      #add.ident = NULL, use.scale = TRUE, use.raw = FALSE,
      #show.progress = TRUE)

      #sel_CellBarCode <- "NORM_P17_WR_Joyal_r1_CTAAGCTAAAAA"

      average_cell_Type <- retina_CTL_celltype.scaled %>% dplyr::select(paste(sel_CellBarCode))

      gene_names <- rownames(average_cell_Type)

      average_cell_Type <- unlist(average_cell_Type, recursive = TRUE, use.names = TRUE)

      names(average_cell_Type) <- gene_names

    
      #plot <- plotEnrichment(GO_term[[1]], average_cell_Type) + labs(title=paste(celltypename,GOI,"_EnrichmentPlot.scaled.png", sep="_"))
    
      #ggsave(plot, file=paste(celltypename,"EnrichmentPlot.scaled.png", sep="_"))


      fgseaRes <- fgsea(pathways = GO_term, stats = average_cell_Type,
                    minSize=15,
                    maxSize=500,
                    nperm=1000,
                    nproc=16,
                    gseaParam=1)
      }
      else {
      fgseaRes <- fgsea(pathways = GO_term, stats = retina_CTL_celltype.scaled,
                    minSize=15,
                    maxSize=500,
                    nperm=1000,
                    nproc=16)
      }
      fgseaRes <-data.frame(fgseaRes, cell_type=celltypename, CellBarCode=sel_CellBarCode)
      gsea <- rbind(gsea,fgseaRes)
    }
  df_total <- rbind(df_total,gsea)
  }


  fwrite(df_total, file =paste("GSEA_CellBarCode", GOI, Condition_subset,"scaled.txt", sep="_"), sep="\t")

  #df_total <- fread("GSEA_Batch_GOautophagy.txt", header=T)

  #Make nice bar plot


  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                        conf.interval=.95, .drop=TRUE) {
      library(plyr)

      # New version of length which can handle NA's: if na.rm==T, don't count them
      length2 <- function (x, na.rm=FALSE) {
          if (na.rm) sum(!is.na(x))
          else       length(x)
      }

      # This does the summary. For each group's data frame, return a vector with
      # N, mean, and sd
      datac <- ddply(data, groupvars, .drop=.drop,
        .fun = function(xx, col) {
          c(N    = length2(xx[[col]], na.rm=na.rm),
            mean = mean   (xx[[col]], na.rm=na.rm),
            sd   = sd     (xx[[col]], na.rm=na.rm)
          )
        },
        measurevar
      )

      # Rename the "mean" column    
      datac <- rename(datac, c("mean" = measurevar))

      datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

      # Confidence interval multiplier for standard error
      # Calculate t-statistic for confidence interval: 
      # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
      ciMult <- qt(conf.interval/2 + .5, datac$N-1)
      datac$ci <- datac$se * ciMult

      return(datac)
  }


  tgc <- summarySE(df_total, measurevar="NES", groupvars="cell_type")

  tgc_padj <- summarySE(df_total, measurevar="padj", groupvars="cell_type")

  tgc2 <- merge(tgc,tgc_padj,by="cell_type")

  tgc2$cell_type <- factor(tgc2$cell_type)


  png(paste("BarPlotNES.Retina", GOI, Condition_subset, "SEM_black.scaled.png", sep="_"), width=3000, height=2000, res=300)
  ggplot(tgc2, aes(x=reorder(cell_type,-NES), y=NES, fill=NULL)) + 
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      xlab("Cell Type") +
      ylab("Mean of Normalized Enrichement Score (± SEM)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      scale_fill_distiller(palette = "YlOrRd", direction=-1, trans = "log10", breaks = c(0, 0.01, 0.05, 0.2, 0.8, 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      coord_flip()+ 
      theme_classic() +
      theme(legend.position="bottom") 
  dev.off()

  png(paste("BarPlotNES.Retina", GOI, Condition_subset, "SEM.scaled.colors.png", sep="_"), width=2000, height=2000, res=300)
  ggplot(tgc2, aes(reorder(cell_type, NES),NES)) +
      geom_bar(aes(fill = padj), stat="identity") +
      geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9),
                    color="brown")+
      xlab("Cell Type NORM P17") +
      ylab("Mean of Normalized Enrichement Score (± SEM)")+
      scale_fill_distiller(palette = "Reds", direction=-1, trans = "log10") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      coord_flip()+ 
      theme_classic() +
      theme(legend.position="bottom") 
  dev.off()


  tgc2$Condition <- Condition_subset

  tgc2$cell_type <- paste(Condition_subset, tgc2$cell_type, sep="_")

  assign(paste("tgc2", Condition_subset, sep="."),tgc2)


  #rbPal <- colorRampPalette(c('red','blue'))

  #tgc2.OIR$Colour <- rbPal(100)[as.numeric(cut(tgc2.OIR$padj,breaks = 100))]



  write.table(tgc2, paste("tgc2", GOI, Condition_subset, "celltype_NES_summary.scaled.txt", sep="_"), sep="\t")


  aov.dfb <- aov(NES ~ cell_type, data=df_total)

  posthoc <- TukeyHSD(x=aov.dfb, "cell_type",  conf.level=0.95)


  summary(aov.dfb)

  print(model.tables(aov.dfb,"means"),digits=3)


  posthoc.save <- as.data.frame(posthoc$cell_type)

  write.table(posthoc.save, paste("Anova_Tukey_Retina", GOI, Condition_subset, "cellbarcode.scaled.txt", sep="_"), sep="\t")

}

###Mergin NORM and OIR

tgc.merged <- rbind(tgc2.NORM,tgc2.OIR)


tgc2.NORM <- tgc2.NORM[order(tgc2.NORM$NES,decreasing=TRUE),]

tgc2.NORM_cell <- as.character(tgc2.NORM$cell_type)
tgc2.NORM_NES <- tgc2.NORM$NES

#
tgc2.OIR_cell <- as.character(tgc2.OIR$cell_type)
tgc2.OIR_NES <- tgc2.OIR$NES

#
order_cell <- insert(tgc2.NORM_cell, "",ats = 2:length(tgc2.NORM_cell))
order_n <- insert(tgc2.NORM_NES, "",ats = 2:length(tgc2.NORM_NES))
order_cell[length(order_cell)+1]<-""
order_n[length(order_n)+1]<-""

#
i <- 1
while(1){
  tmp <- tgc2.OIR_cell[grep(sub("_","",str_match(order_cell[i],"_.*")),tgc2.OIR_cell)]
  order_cell[i+1] <- tmp
  order_n[i+1] <- tgc2.OIR$NES[which(tgc2.OIR$cell_type==tmp)]
  i <- i + 2
  if(i>length(order_cell)) break
}

#
ordered.df <- as.data.frame(cbind(order_cell,order_n))

ordered.df[,2] <- factor(ordered.df[,2], levels = unique(ordered.df[,2]))

tgc.merged$cell_type <- factor(tgc.merged$cell_type, level=rev(unique(ordered.df[,1])))

tgc.merged$NES <- as.numeric(tgc.merged$NES)


png(paste("BarPlotNES.Retina_", GOI, "_NORMvsOIR_SEM.scaled.colors.new.png", sep=""), width=2000, height=2000, res=300)
ggplot(tgc.merged, aes(cell_type, NES, fill=Condition)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9),
                  color="darkgrey")+
    xlab("Cell Type P17") +
    ylab("Mean of Normalized Enrichement Score (± SEM)")+ 
    #geom_col(fill = tgc.merged$Colour)+
    scale_fill_manual(values=c("limegreen", "firebrick1"))+
    #scale_fill_distiller(palette = "Reds", direction=-1, trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    coord_flip()+ 
    theme_classic() +
    theme(legend.position="bottom")+ 
    ggtitle(paste(GOI))
dev.off()


q("no")
