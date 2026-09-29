# ---------------------------------------------------------------------------
# Step 06 - fgsea on OIR-vs-normoxia ranked gene lists, per cell type.
#
# ESTABLISHES
#     Normalised enrichment scores per cell type for selected pathways, with the
#     per-cell-type NES summary and the ANOVA/Tukey comparison across cell types.
#
# WHY THIS WAY
#     The ranking statistic comes from FindMarkers run with
#     logfc.threshold = 0 and min.pct = 0 - deliberately unfiltered. fgsea needs
#     a complete ranked list; pre-filtering by fold change or detection rate
#     would truncate the ranking and bias the enrichment score. This is the
#     opposite setting from step 04, where filtering is appropriate because the
#     output is a gene table rather than a ranking.
#
# REQUIRES
#     step 02
#
# INPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#     <MSIGDB_DIR>/{h.all,c2.cp,c5.all}.v7.1.symbols.gmt
#     <GENE_LIST_DIR>/ CLEAR_NETWORK.csv, TFEB_TARGET_LIST.txt,
#                      GO_MITOPHAGY.txt, REACTOME_PINK_PARKIN_MEDIATED_MITOPHAGY.txt,
#                      GO_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION*.txt,
#                      WP_NAD_METABOLISM_SIRTUINS_AND_AGING.txt
#
# OUTPUTS
#     GSEA_<pathway>.txt                     per-cell-type fgsea results
#     tgc2<pathway>_*.celltype_NES_summary.txt
#     Anova_Tukey_Retina_<pathway>.txt
#
# WHAT WAS REMOVED FROM THIS FILE
#     Source lines 353-596 of the cluster original were a second analysis
#     contrasting ident.1 = "VldlrKO" against ident.2 = "CTL" on a per-cell-type
#     basis. Those identities do not exist in this object, and none of that
#     section's outputs (*Vldlr* files, barplotNESPvalue.PR.*) are present in the
#     fgsea directory on the cluster, so it cannot have run on these cells. It is
#     removed here so the published script is the code that produced the
#     published figures. See METHODS_CHOICES.md.
#
#     Not removed, and not a mistake: the gene set named JOYAL_FAO_VldlrKO_Retina
#     defined around line 193. That is a fatty-acid-oxidation signature derived
#     from a Vldlr-KO retina study and used here as one of the pathways tested in
#     OIR. The name refers to the signature's origin, not to a contrast.
#
#     Also left in place: one write.table inside the OIR section names its output
#     tgc2<GOI>_VldlrKOvsCTL.celltype_NES_summary.txt, which mislabels an
#     OIR-vs-NORM result. No such file exists on the cluster, so no published
#     result carries that name.
#
# SURVIVING OUTPUTS ON THE CLUSTER
#     Two .eps barplots per pathway, four pathways, eight files in total:
#     <GOI>/1000Cells/BarPlotNES.Retina_<GOI>OIRvsNORM_{P14,P17}.eps
#     The GSEA_<GOI>.txt, tgc2* and Anova_Tukey_* tables this script writes are
#     not present in any pathway directory. Either only the sections up to the
#     barplot were executed or the tables were removed afterwards; the files do
#     not say which, and the analysis was not re-run to find out.
#
# USAGE
#     Rscript 06_fgsea.R
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
OUT_DIR <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/DifferentialExpression/fgsea")

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

# GOI is set once, low in the body of the script, and the analysis was run
# once per pathway by editing that line between runs - four runs, four output
# directories. The four values swept are recorded here because the script alone
# would tell a reader only one pathway was examined. This vector is
# documentation: the body is unchanged and still uses its own single GOI.
GOI_SWEPT <- c("GO_AEROBIC_RESPIRATION",
               "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS",
               "SIRT3_TARGET_GENES",
               "WP_NAD_METABOLISM_SIRTUINS_AND_AGING")

FGSEA_LOGFC_MIN <- 0      # unfiltered on purpose: fgsea needs a full ranking
FGSEA_MIN_PCT   <- 0
FGSEA_TEST      <- "wilcox"

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
library(fgsea)
library(GSEABase)
library(GSVAdata)
library(GSVA) 
data(c2BroadSets)
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

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Cond_TimePoint %in% c("OIR_P17", "NORM_P17")))

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

#Seurat_object_Subset <- SubsetData(Seurat_object_Subset, ident.remove = "Opticin_cells")


###Max 1000 cells per ID
retina_subset <- SubsetData(Seurat_object_Subset, cells.use = NULL, subset.name =NULL, ident.use = NULL,
  ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
  do.center = TRUE, do.scale = TRUE, max.cells.per.ident = 1000,
  random.seed = 1)


#Select pathway

GOI <- "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS"

dir.create(file.path(OUT_DIR, GOI))

setwd(file.path(OUT_DIR, GOI))
dir.create("1000Cells")
setwd("1000Cells")




###Load gene set collection

broadset_selected <- "sup_h.c2.c5.c6.c7.all.270519.symbols"

broadSet.sup_h.c2.c5.c6.c7.all.270519 <- getGmt(paste(file.path(GENE_LIST_DIR/Gmt.file, ""), broadset_selected, ".gmt", sep=""),
              geneIdType=SymbolIdentifier())



broadset_selected <- "c2.cp.v7.1.symbols"

broadSet.c2.cp <- getGmt(paste(file.path(MSIGDB_DIR, ""), broadset_selected, ".gmt", sep=""),
              geneIdType=SymbolIdentifier())


broadset_selected <- "c5.all.v7.1.symbols"

broadSet.c5.all <- getGmt(paste(file.path(MSIGDB_DIR, ""), broadset_selected, ".gmt", sep=""),
              geneIdType=SymbolIdentifier())



broadset_selected <- "h.all.v7.1.symbols"

broadSet.h.all <- getGmt(paste(file.path(MSIGDB_DIR, ""), broadset_selected, ".gmt", sep=""),
              geneIdType=SymbolIdentifier())


##Create custom geneset collection or via gene set list


custom_geneset = c("ACADVL", #JOYAL_FAO_VldlrKO_Retina
"CPT1A", "PPARA", "ALDH9A1", "ECHS1", "ACAA2", "CPT2", "ECI2", "ACADL", "MLYCD", "AKT1", "PPARGC1A", 
"FABP1", "FABP2","FABP3","FABP4","FABP5","FABP6","SLC27A1","SLC27A2","SLC27A3","SLC27A4","SLC27A5","SLC27A6",
"HMGCL", "CD36"
#"ACSL1","ACSL2","ACSL3","ACSL4","ACSL5","ACSL6"
)

genesetcol_1 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="JOYAL_FAO_VldlrKO_Retina"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "GO_MITOPHAGY.txt"), header=F, fill=TRUE)$V1[-(1:2)]

genesetcol_2 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="GO_MITOPHAGY"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "REACTOME_PINK_PARKIN_MEDIATED_MITOPHAGY.txt"), header=F, fill=TRUE)$V1[-(1:2)]

genesetcol_3 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="REACTOME_PINK_PARKIN_MEDIATED_MITOPHAGY"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "GO_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION_IN_RESPONSE_TO_MITOCHONDRIAL_DEPOLARIZATION.txt"), header=F, fill=TRUE)$V1[-(1:2)]

genesetcol_4 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="GO_REGULATION_OF_AUTOPHAGY_OF_MITO_DEPOLARIZATION"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "GO_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION.txt"), header=F, fill=TRUE)$V1[-(1:2)]

genesetcol_5 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="GO_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "CLEAR_NETWORK.csv"), header=F, fill=TRUE)$V1

genesetcol_6 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="CLEAR_NETWORK"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "TFEB_TARGET_LIST.txt"), header=F, fill=TRUE)$V1[-(1:2)]

genesetcol_7 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="TFEB_TARGET_GENES"))


custom_geneset <- fread(file.path(GENE_LIST_DIR, "WP_NAD_METABOLISM_SIRTUINS_AND_AGING.txt"), header=F, fill=TRUE)$V1[-(1:2)]

genesetcol_8 <- GeneSetCollection(GeneSet(custom_geneset, geneIdType=SymbolIdentifier(), setName="WP_NAD_METABOLISM_SIRTUINS_AND_AGING"))

##Make genesetcollection by combining geneset

gsc_final <- GeneSetCollection(c(broadSet.h.all,broadSet.c2.cp,broadSet.c5.all,genesetcol_8))


###Select genes from genesetcollection

GO_term <- grep("REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS", names(gsc_final))
GO_term <- names(gsc_final)[GO_term]
GO_term <- geneIds(gsc_final[[GO_term]])

GO_NAD_DEPENDENT_PROTEIN_DEACETYLASE_ACTE_ACTIVITY
GO_NADPLUS_BINDING
GO_NAD_DEPENDENT_PROTEIN_DEACETYLASE_ACTIVITY


#or via specific gene set

GO_term <- read.table(paste(file.path(GENE_LIST_DIR, ""), GOI, ".txt", sep=""), sep="\t")
GO_term <- GO_term$V1


# Prepare the gene list for GSEA


GO_term <- intersect(rownames(retina_subset),GO_term)


######


###############################

###Perform fgsea on ranked dataset of fold-change between OIR vs NORM


celltypes <- names(summary(retina_subset@active.ident))


df_total = NULL


for (celltypename in celltypes) {
  #celltypename <- "Cones"
  retina_subset_celltype <- SubsetData(retina_subset, cells.use = NULL, subset.name =NULL, ident.use = celltypename,
    ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
    do.center = FALSE, do.scale = FALSE, max.cells.per.ident = Inf,
    random.seed = 1)

  print(summary(retina_subset_celltype@active.ident))

  retina_subset_celltype <- SetIdent(retina_subset_celltype, value="Condition")

  DGEs <- FindMarkers(retina_subset_celltype, ident.1="OIR", ident.2 = "NORM", genes.use = NULL,
  logfc.threshold = 0, test.use = "wilcox", min.pct = 0,
  min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
  max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
  min.cells.gene = 3, min.cells.group = 3, pseudocount.use = 1,
  assay.type = "RNA")

  DGEs_df <- data.frame(DGEs$avg_logFC)
  rownames(DGEs_df) <- rownames(DGEs)

  gene_names <- rownames(DGEs_df)

  DGEs_df <- unlist(DGEs_df, recursive = TRUE, use.names = TRUE)

  names(DGEs_df) <- gene_names

  #plot <- plotEnrichment(GO_term[[1]], DGEs_df) + labs(title=paste(celltypename,GOI,"_EnrichmentPlot.scaled.png", sep="_"))
  
  #ggsave(plot, file=paste(celltypename,"EnrichmentPlot.png", sep="_"))


  GO_term_genes <- list(intersect(names(DGEs_df),GO_term))

  names(GO_term_genes) <- GOI

  fgseaRes <- fgsea(pathways = GO_term_genes, stats = DGEs_df,
                  minSize=2,
                  maxSize=1000,
                  nperm=1000,
                  nproc=16)
  
  fgseaRes <- data.frame(fgseaRes, cell_type=celltypename)
  gsea=NULL
  gsea <- rbind(gsea,fgseaRes)

df_total <- rbind(df_total,gsea)
}


fwrite(df_total, file =paste("GSEA_", GOI, ".txt", sep=""), sep="\t")

#df_total <- fread(paste("GSEA_", GOI, ".txt", sep=""), header=T)

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


###

tgc <- summarySE(df_total, measurevar="NES", groupvars="cell_type")

tgc_padj <- summarySE(df_total, measurevar="padj", groupvars="cell_type")

tgc <- merge(tgc,tgc_padj,by="cell_type")

tgc2 <- tgc
tgc2$cell_type <- factor(tgc2$cell_type)



bar <- ggplot(tgc2, aes(reorder(cell_type, NES),NES)) +
    geom_bar(aes(fill = padj), stat="identity") +
    #geom_errorbar(aes(ymin=NES-se.x, ymax=NES+se.x),
    #              width=.2,                    # Width of the error bars
    #              position=position_dodge(.9),
    #              color="brown")+
    scale_fill_distiller(palette = "RdBu", direction=1, trans = "log10", breaks = c(0.01, 0.05, 0.5, 1)) +
    #scale_fill_viridis(option="viridis", direction=-1, trans = "log10", breaks = c(0, 0.01, 0.05,0.2,0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    coord_flip()+ 
    theme_classic() +
    theme(legend.position="bottom") 

ggsave(paste("BarPlotNES.Retina_", GOI, "OIRvsNORM_P17.eps", sep=""), dpi = 300, width = 4, height = 4 )




write.table(tgc2, paste("tgc2", GOI,"_VldlrKOvsCTL.celltype_NES_summary.txt", sep=""), sep="\t")


aov.dfb <- aov(NES ~ cell_type, data=df_total)

posthoc <- TukeyHSD(x=aov.dfb, "cell_type",  conf.level=0.95)


summary(aov.dfb)

print(model.tables(aov.dfb,"means"),digits=3)


posthoc.save <- as.data.frame(posthoc$cell_type)

write.table(posthoc.save, paste("Anova_Tukey_Retina_", GOI,".txt", sep=""), sep="\t")




q("no")
