# ---------------------------------------------------------------------------
# 06_senescence - 06_cellphonedb_oir_p17_kras.py
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# the senescence-merged compartment, run on cells taken from the annotated
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
#
# LANGUAGE - this file was named .py on the cluster and is not python
#     Measured: 82 R statements, 9 shell lines, 0 python statements. It is an R
#     script with the CellPhoneDB command line embedded in it as a shell block
#     (the `cd` / `COUNTS=` / `METADATA=` / `cellphonedb method statistical_analysis`
#     lines). Renamed to .R here; the shell block is left in place because it is
#     the only record of how CellPhoneDB was invoked, and it is the step that
#     runs between writing the input and reading the output back.
#
#     or write into directories that no longer exist. See ../../EDIT_POLICY.md.

# ===========================================================================
# 0. parameters
# ===========================================================================
PROJ_DIR <- "/project/def-jsjoyal/gaelcge"
TC_DIR   <- file.path(PROJ_DIR, "Retina/OIR_NORM_TimeCourse")
B_DIR    <- file.path(TC_DIR, "Aligned/Cond_Sorting")
# the python 3.7 virtualenv CellPhoneDB ran in; it no longer exists
PY37_ENV <- "~/ENV_python3.7.0/bin/activate"

##In R

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
library(doubletFinder)
library(magrittr)
library(harmony)
library(SingleR)
library(ktplots)


###Make input data from Seurat object for cellphonedb

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/Mapping"))

retina <- readRDS("retina_subclustered_renamed.rds")

##Subset data (optional)

retina <- SetAllIdent(retina, id = "KRAS_cell_subtype")

retina <- SubsetData(retina, cells.use = NULL, subset.name = NULL, ident.use = NULL, ident.remove = c("Immune cells 1", "Immune cells 2", "Immune cells 4", "Immune cells 5"))

retina <- SetAllIdent(retina, id = "Cond_TimePoint")

retina <- SubsetData(retina, cells.use = NULL, subset.name = NULL, ident.use = "OIR_P17", ident.remove = NULL)

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/input/OIR_P17_KRAS"))

#From cellphoneDB v2 preprint
count_raw <- retina@raw.data[,retina@cell.names]
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
write.table(count_norm, "count_OIR_P17_KRAS.txt", sep="\t", quote=F)

#Using Seurat norm log data (not recommanded)
#count <- retina@data[,retina@cell.names]
#count <- apply(count, 2, function(x) (x*1))
#write.table(count, "log_norm_count.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote=F)

## generating meta file
meta_data <- cbind(rownames(retina@meta.data), retina@meta.data[,"KRAS_cell_subtype", drop=F]) # cluster is the user’s corresponding cluster column
write.table(meta_data, "metadata_OIR_P17_KRAS.txt", sep = "\t", col.names = TRUE, row.names = FALSE)


   ##########################
   ##########################


###Run cellPhoneDb in python environement

source ENV_python3.7.4/bin/activate

cd B_DIRSubclustering/MergingSenescence/CellphoneDB/input/OIR_P17_KRAS

DATE=$(date +"%d.%m.%y_%H_%M_%S")

COUNTS=count_OIR_P17_KRAS.txt
METADATA=metadata_OIR_P17_KRAS.txt

OUT=B_DIRSubclustering/MergingSenescence/CellphoneDB/output/

CPUS=$SLURM_CPUS_PER_TASK

echo "Date:" $DATE
echo "Counts:" ${COUNTS[${SLURM_ARRAY_TASK_ID}]}
echo "Meta data:" ${METADATA[${SLURM_ARRAY_TASK_ID}]}
echo "Output:" ${OUT}
echo "nCpus:" ${CPUS}

cellphonedb method statistical_analysis ${METADATA} ${COUNTS} --project-name=OIR_P17_KRAS --output-path ${OUT} --threads ${CPUS} --counts-data=gene_name




cellphonedb plot dot_plot --means-path=${OUT}/OIR_P17/means.txt --pvalues-path=${OUT}/OIR_P17/pvalues.txt --verbose


cellphonedb plot heatmap_plot ${METADATA} --pvalues-path=${OUT}/OIR_P17/pvalues.txt

###

####################################################


                    #### make plot in R####

####################################################


##    Making dotplot using Joel function

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/output/OIR_P17_KRAS"))

library(ggplot2)

dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){

  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, family = 'Arial'),
        axis.text.y = element_text(size=12, colour = "black", family = 'Arial'),
        axis.title=element_blank(),
        text = element_text('Arial'),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if (output_extension == '.pdf') {
      ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
      ggsave(filename, width = width, height = height, limitsize=F, dpi = 300)
  }
}


pvals = read.table(paste0("pvalues.txt"), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
means = read.table(paste0("means.txt"), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)

cell_interactions <- c("Immune cells 3|KRAS negative Endothelial cells"
                      , "KRAS negative Endothelial cells|Immune cells 3"
                      , "Immune cells 3|KRAS positive Endothelial cells"
                      , "KRAS positive Endothelial cells|Immune cells 3"
#                     , "Immune cells 3|KRAS positive Pericytes"
#                     , "KRAS positive Pericytes|Immune cells 3"
#                     , "Immune cells 3|KRAS positive Astrocytes"
#                     , "KRAS positive Astrocytes|Immune cells 3"
#                     , "Immune cells 3|KRAS positive Muller glia"
#                     , "KRAS positive Muller glia|Immune cells 3"
                      )
cell_interactions <- c(cell_interactions, "interacting_pair")

selected_pvals <- pvals[,cell_interactions]

sig_interactions <- NULL
for(i in 1:nrow(selected_pvals)){
  if(any(selected_pvals[i,]<0.01)){
    sig_interactions <- c(sig_interactions, pvals$interacting_pair[i])
  }
}


dir <- file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/result_plot/OIR_P17_KRAS/")

dot_plot(filename = paste0(dir,"dotplot_OIR_P17.png"), means_path=paste0("means.txt"), pvalues_path=paste0("pvalues.txt"), 
  output_extension=".png", selected_columns = cell_interactions[-length(cell_interactions)], selected_rows = sig_interactions, 
  height = 0.4*length(sig_interactions))





#############   Making dotplot Using kplots


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
library(doubletFinder)
library(magrittr)
library(harmony)
library(SingleR)
library(ktplots)


setwd(file.path(B_DIR, "Subclustering/MergingSenescence/Mapping"))

retina <- readRDS("retina_subclustered_renamed.rds")

##Subset data (optional)

retina <- SetAllIdent(retina, id = "Cond_TimePoint")

retina <- SubsetData(retina, cells.use = NULL, subset.name = NULL, ident.use = "OIR_P17", ident.remove = NULL)

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/output/OIR_P17"))


pvals <- read.delim("pvalues.txt", check.names = FALSE) # pvalues.txt and means.txt are output from cpdb
means <- read.delim("means.txt", check.names = FALSE) 


setwd(file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/result_plot/OIR_P17"))


png("plot_cpdb.png")
plot_cpdb(cell_type1 = "Bcell", # cell_type1 and cell_type2 will call grep, so this will accept regex arguments. Additional options for grep can be specified, such as fixed = TRUE to stop grep from misinterpreting/converting symbols.
  cell_type2 = "Tcell",
  retina,
  "cell_subtype", # column name where the cell ids are located in the metadata
  means,
  pvals,
  split.by = "group", # column name where the grouping column is
  genes = c("CXCR3", "CD274", "CXCR5"))
dev.off()



#################
#################


##Making heatmap plot  Using cellphonedb plot

dir <- file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/input/OIR_P17_KRAS/")

meta_file <- paste0(dir,"metadata_OIR_P17_KRAS.txt")  #required

dir <-file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/output/OIR_P17_KRAS/")

pvalues_file <- paste0(dir,"pvalues.txt")  #The pvalues output file [./out/pvalues.txt]

count_filename <- "heatmap_count.png" #Filename of the output plot [heatmap_count.pdf]

log_filename <- "heatmap_log_count.png" #Filename of the output plot using log-count of interactions [heatmap_log_count.pdf]

count_network_filename <- "count_network.csv" #Filename of the output network file [count_network.txt]

interaction_count_filename <- "interaction_count.txt" #Filename of the output interactions-count file [interactions_count.txt]

setwd(file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/result_plot/OIR_P17_KRAS"))

source(file.path(B_DIR, "Subclustering/MergingSenescence/CellphoneDB/result_plot/heatmap_cellphonedb.R"))


heatmaps_plot(meta_file, pvalues_file, 
              count_filename=count_filename, 
              log_filename =log_filename, 
              count_network_filename=count_network_filename, 
              interaction_count_filename=interaction_count_filename,
              count_network_separator=",",
              interaction_count_separator=" = ", 
              show_rownames = T, show_colnames = T,
              scale="none", cluster_cols = T, border_color='white', cluster_rows = T, fontsize_row=11,
              fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
              col1 = "slateblue", col2 = 'yellow', col3 = 'red', meta_sep='\t', pvalues_sep='\t', pvalue=0.05)



##Select cell interaction for heatmap

count_network_filename_result <- read.csv(count_network_filename)

count_network_filename_result$interaction <- paste(count_network_filename_result$SOURCE,count_network_filename_result$TARGET,sep="_")

count_network_filename_result <-count_network_filename_result[2:9,]

my_palette <- colorRampPalette(c("blue", "yellow", "red"))

png("heatmap_immune3_interaction.png", width=2500, height=2000, res = 300)
heatmap.2(cbind(count_network_filename_result$count, count_network_filename_result$count), 
          dendrogram = "row", labCol = "", labRow = count_network_filename_result$interaction, col=my_palette(50), 
          symbreak=FALSE, trace='none', cexRow=1, cexCol= 1, Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(4,10), lwid=c(4,8), 
          key.par=list(mar=c(5, 3, 3, 3)), margins=c(5,30), key.ylab=NA, srtCol=45)
dev.off()


setEPS()
postscript("heatmap_immune3_interaction.eps", width=8, height=6)
heatmap.2(cbind(count_network_filename_result$count, count_network_filename_result$count), 
          dendrogram = "row", labCol = "", labRow = count_network_filename_result$interaction, col=my_palette(50), 
          symbreak=FALSE, trace='none', cexRow=1, cexCol= 1, Colv="Rowv", revC=TRUE, key.title = NA, density.info="none", lhei=c(4,10), lwid=c(4,8), 
          key.par=list(mar=c(5, 3, 3, 3)), margins=c(5,30), key.ylab=NA, srtCol=45)
dev.off()










