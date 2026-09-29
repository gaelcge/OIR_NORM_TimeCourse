# ---------------------------------------------------------------------------
# 02_neuroglial - 03_monocle_trajectory.R
#
# Part of 2_compartments_p14_p17. Subclustering and downstream analysis of
# neurons and glia, run on cells taken from the annotated
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
#     or write into directories that no longer exist. See ../../EDIT_POLICY.md.
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
options(warn=-1)
library(plotly)
library(tidyverse)
library(ggrepel)
library(gtools)
library(data.table)
library(gplots)
library(useful)
library(locfit)
library(Matrix)
library(monocle)
library(DDRTree)
library(pheatmap)
library(reshape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(magrittr)
library(harmony)
library(RColorBrewer)
library(Seurat)

#Set directory of dataset to analyse


setwd(file.path(B_DIR, "Subclustering/NeuroGlialCell/CellIdentification"))

NeuroGlia <- readRDS("NeuroGlia_subclustered_renamed.rds")

summary(NeuroGlia@ident)
summary(NeuroGlia@meta.data)


#Start monocle analysis


setwd(file.path(B_DIR, "Subclustering/NeuroGlialCell/Monocle"))
## Load the expression matrix 



##############################

####Run monocle

# Where 'data_to_be_imported' can either be a Seurat object
# or an SCESet.

#importCDS(data_to_be_imported)

# We can set the parameter 'import_all' to TRUE if we'd like to
# import all the slots from our Seurat object or SCESet.
# (Default is FALSE or only keep minimal dataset)

HSMM <- importCDS(NeuroGlia, import_all = TRUE)

cds <- updateCDS(HSMM)

##Or convert from raw data

HSMM_expr_matrix <- NeuroGlia@raw.data
HSMM_expr_matrix <- HSMM_expr_matrix[,colnames(NeuroGlia@data)]
HSMM_sample_sheet <- NeuroGlia@meta.data
HSMM_gene_annotation <- as.data.frame(HSMM_expr_matrix[,1])
HSMM_gene_annotation[,1] <- rownames(HSMM_gene_annotation)
colnames(HSMM_gene_annotation) <- "gene_short_name"

#Once these tables are loaded, you can create the CellDataSet object like this.
#Don't accidentally convert to a dense expression matrix
#The output from a number of RNA-Seq pipelines, including CellRanger, is already in a sparseMatrix format (e.g. MTX). 
#If so, you should just pass it directly to newCellDataSet without first converting it to a dense matrix (via as.matrix(), 
#because that may exceed your available memeory.

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
cds <- newCellDataSet(as.matrix(HSMM_expr_matrix),
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.5)

#Step 1: Noramlize and pre-process the data
#We then estimate size factors for each cell and dispersion function for the genes in the cds as usual. Monocle 3 now performs these operations using the DelayedArray packages so they work on datasets with millions of cells. The dispersion calculation and several other operations rely on the DelayedArray package in Bioconductor, which splits the operation into blocks in order to avoid exhausting the computer's memory. You can control the block size and the verbosity of these operations as shown below:

# Pass TRUE if you want to see progress output on some of Monocle 3's operations
DelayedArray:::set_verbose_block_processing(TRUE)

#####Estimate size factors and dispersions-Required
#We'll  call two functions that pre-calculate some information about the data. 
#Size factors help us normalize for differences in mRNA recovered across cells, and "dispersion" values will help us perform differential expression analysis later.

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


#Next, run the preprocessCDS() function to project the data onto the top principal components:

cds <- preprocessCDS(cds, num_dim = 20)

#Step 2: Reduce the dimensionality of the data
#Then, apply a further round of (nonlinear) dimensionality reduction using UMAP:

cds <- reduceDimension(cds, reduction_method = 'UMAP')


#Now we're ready to start learning trajectories!

#Step 3: Partition the cells into supergroups
#Rather than forcing all cells into a single developmental trajectory, Monocle 3 enables you to learn a set of trajectories that describe the biological process you're studying. For example, if you're looking at a community of immune cells responding to infection, each cell type will respond to antigen (and each other) in a different way, so they should be organized into distinct trajectories. This can also be helpful when you have small groups of outlier cells that for either technical or biological reasons are very dissimilar from the rest of the cells in your experiment. They can confuse a trajectory analysis. Monocle 3's partitioning strategy circumvents this issue because such groups often wind up in their own partition. In the Paul data, there is a small outgroup the authors classified as dendritic cells that Monocle automatically partitions away from the main trajectory.

#In Monocle 3, we recognize "disjoint" trajectories by drawing on ideas from Alex Wolf and colleagues, who recently introduced the concept of abstract graph participation. Monocle 3 implements the test for cell community connectedness from Wolf et al via the partitionCells() function, which divides the cells into "supergroups".

cds <- partitionCells(cds)

#Step 4: Learn the principal graph
#Now that the cells are paritioned, we can organize each supergroup into a separate trajectory. The default method for doing this in Monocle 3 is SimplePPT, which assumes that each trajectory is a tree (albeit one that may have multiple roots). Learn these trees with the learnGraph function:

cds <- learnGraph(cds,  RGE_method = 'SimplePPT')

#After you've learned the graph, it's time to assign each cell's pseudotime value with orderCells(). But before we get to that, let's see how to plot the trajectory so we know where the beginning of the trajectory should be.

#Step 5: Visualize the trajectory
#Once the you've learned the trajectory, you can visualize the it and color the cells in different ways.

print(head(pData(cds)))

png("Plot_cell_trajectory.Cond_TimePoint.png", width=2000, height=1500, res=150)
plot_cell_trajectory(cds,
                     color_by = "Cond_TimePoint")
dev.off()

png("Plot_cell_trajectory.Cond_Cell_type_final.png", width=2000, height=1500, res=150)
plot_cell_trajectory(cds,
                     color_by = "Cell_type_final")
dev.off()

png("Plot_cell_trajectory.Cond_Condition.png", width=2000, height=1500, res=150)
plot_cell_trajectory(cds,
                     color_by = "Condition")
dev.off()

cds <- clusterCells(cds,
                        method = 'louvain',
                        res = 1e-6,
                        louvain_iter = 1,
                        verbose = T,
                        python_home="")
    

png("Plot_cell_trajectory.Cond_NDUFA4L2.png", width=2000, height=1500, res=150)
plot_cell_clusters(cds,
                   color_by = 'Tissue',
                   cell_size = 0.1,
                   show_group_id = T)  +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

#Adjusting the start of pseudotime with orderCells
#In the trajectory above, you can see that the "Multipotent progenitor" cells are located roughly in the middle of the longer trajectory segment. We know from extensive past study of hematopoeisis that these are the "root" cell state that generates all the others in the primary trajectory. We need to tell Monocle that these cells are the "beginning" of the trajectory. In Monocle 2, this wouldn't be possible, because the software required that the root be one of the leaves of the tree. Monocle 3 allows you to specify an internal part of the tree as the root. There are two ways to do this:

#You can specifiy the name of a specific cell or principal graph node as the root by passing it to orderCells().
#You can call orderCells() with no root specified and it will open a window for you to select the root(s) with your mouse cursor. This latter option is only available in interactive sessions, and doesn't work in Jupyter notebooks.
#Once you select one or more roots, orderCells() computes the shortest path from each cell's location on the princiapl graph to the nearest root node. That is, a cell's pseudotime value is the geodesic distance from it to the nearest root, traveling over the graph. Any cell that is not reachable from some root will be assigned a pseudotime value of infinity.

#You may find it helpful to automatically pick the root according to any number of biologically-driven criteria. For example, you could find the nodes at which cells expressing a certain marker gene are concentrated. Or we could select the node where cells from an early experimental timepoint land. Here, we provide one such helper function to show you how to do this kind of thing:

# a helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)

  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

#In the above function, we are accessing pr_graph_cell_proj_closest_vertex which is just a matrix with a single column that stores for each cell, the ID of the principal graph node it's closest to. This is handy for computing statistics (e.g. with dplyr) about the principal graph nodes and which cells of what type map to them.

#Now we can call this function to automatically find the node in the principal graph where our multipotent progenitors reside:

MPP_node_ids = get_correct_root_state(cds,
                                      cell_phenotype =
                                        'Cell_type_final', "Muller Glia 1")

cds <- orderCells(cds, root_pr_nodes = MPP_node_ids)



png("Plot_cell_trajectory.Cond_Cell_type_final.ordered.png", width=2000, height=1500, res=150)
plot_cell_trajectory(cds)
dev.off()


#Learning and visualizing trajectories in 3D
#Sometimes, projecting cells into three dimensions instead of two can make the biological process you're studying easier to interpret. In Monocle 3, we provide the plot_3d_cell_trajectory function to plot a dataset in 3 dimensions. The block of code below shows you how to learn a trajectory in 3D. Similarly to the 2d case, it:

#Reduces the data down into three dimensions using UMAP
#Learns and the trajectory trajectory and orders the cells
#Visualizes the 3D trajectory

cds <- reduceDimension(cds, max_components = 3,
                       reduction_method = 'UMAP',
                       metric="cosine",
                       verbose = F)


cds <- partitionCells(cds)

cds <- learnGraph(cds,
                  max_components = 3,
                  RGE_method = 'SimplePPT',
                  partition_component = T,
                  verbose = F)

cds <- orderCells(cds,
                  root_pr_nodes =
                    get_correct_root_state(cds,
                                           cell_phenotype = 'Cell_type_final', 
                                           "Muller Glia 1"))


plot_3d_cell_trajectory(cds,
                        color_by="Cell_type_final",
                        webGL_filename=
                          paste(getwd(), "/trajectory_3D.html", sep=""),
                        show_backbone=TRUE,
                        useNULL_GLdev=TRUE)


plot_3d_cell_trajectory(cds,
                        color_by="Condition",
                        webGL_filename=
                          paste(getwd(), "/trajectory_3D.Condition.html", sep=""),
                        show_backbone=TRUE,
                        useNULL_GLdev=TRUE)


#Identifying genes that vary in expression over a trajectory
#We are often interested in finding genes that are differentially expressed across a single-cell trajectory. Monocle 3 introduces a new approach for finding such genes that draws on a powerful technique in spatial correlation analysis, the Moran’s I test. Moran’s I is a measure of multi-directional and multi-dimensional spatial autocorrelation. The statistic tells you whether cells at nearby positions on a trajectory will have similar (or dissimilar) expression levels for the gene being tested. Although both Pearson correlation and Moran’s I ranges from -1 to 1, the interpretation of Moran’s I is slightly different: +1 means that nearby cells will have perfectly similar expression (as in the right panel below); 0 represents no correlation (center), and -1 means that neighboring cells will be *anti-correlated* (left).

pr_graph_test <- principalGraphTest(cds, k=3, cores=1)


#We can easily view the top differentially expressed genes as follows:

dplyr::add_rownames(pr_graph_test) %>%
    dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-qval)) %>% head(3)

#You can check out the test results for one or more genes like this:

pr_graph_test[fData(cds)$gene_short_name %in% c("NDUFA4L2"),]

#The code below reports that overall, there are more than XXX DE genes over the whole trajectory:

nrow(subset(pr_graph_test, qval < 0.01))

#Once you've identified differentially expressed genes, you'll often want to visualize their expression levels on the trajectory. The plot bleow shows each cell with detectable levels of Hbb-b1 superimposed on the trajectory. Hbb-b1 is a subunit of beta globin, a highly specific marker of erythroid cells. Indeed, it is largely restricted to the erythroid branch of the trajectory. Values are log-transformed and scaled into Z scores. Cells with no expression are not shown to avoid overplotting.

plot_3d_cell_trajectory(cds, markers = c('NDUFA4L2'),
                        webGL_filename=paste(getwd(), "/trajectory_3D.NDUFA4L2.html", sep=""),
                        show_backbone=TRUE,
                        useNULL_GLdev=TRUE)













#Filtering low-quality cells-Recommended

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
    num_cells_expressed >= 10))


##The vector expressed_genes now holds the identifiers for genes expressed in at least 50 cells of the data set. 
#We will use this list later when we put the cells in order of biological progress. 

#CellDataSet objects provide a convenient place to store per-cell scoring data: the phenoData slot. 
#Simply include scoring attributes as columns in the data frome you used to create your CellDataSet container. 
#You can then easily filter out cells that don't pass quality control. 

print(head(pData(HSMM)))
    
#If you are using RPC values to measure expression, as we are in this vignette, it's also good to look at the distribution of mRNA totals across the cells:

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
            2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
            2*sd(log10(pData(HSMM)$Total_mRNAs)))

png("qplot_retina_monocle.Total_mRNAs.lab.png")
qplot(Total_mRNAs, data = pData(HSMM), color = Labo, geom =
"density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)
dev.off()

png("qplot_retina_monocle.Total_mRNAs.Batch.png", width=1000)
qplot(Total_mRNAs, data = pData(HSMM), color = Batch, geom =
"density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)
dev.off()


png("qplot_retina_monocle.Total_mRNAs.Dataset.png")
qplot(Total_mRNAs, data = pData(HSMM), color = orig.ident, geom =
"density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)
dev.off()


png("qplot_retina_monocle.Total_mRNAs.Cond_TimePoint.png")
qplot(Total_mRNAs, data = pData(HSMM), color = Cond_TimePoint, geom =
"density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)
dev.off()

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
      pData(HSMM)$Total_mRNAs < upper_bound]

#Once you've excluded cells that do not pass your quality control filters, you should verify that 
#the expression values stored in your CellDataSet follow a distribution that is roughly lognormal:
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.

print(head(pData(HSMM)))

table(pData(HSMM)$cell_type)

table(pData(HSMM)$Cell_type_final)

png("qplot_retina_monocle.CellTypeFinal_distribution.png")
pie <- ggplot(pData(HSMM),
aes(x = factor(1), fill = factor(Cell_type_final))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

######
######Constructing Single Cell Trajectories


# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress through myogenesis. We ultimately want a set of genes that increase (or decrease) in expression as a function of progress through the process we're studying.

# Ideally, we'd like to use as little prior knowledge of the biology of the system under study as possible. We'd like to discover the important ordering genes from the data, rather than relying on literature and textbooks, because that might introduce bias in the ordering. We'll start with one of the simpler ways to do this, but we generally recommend a somewhat more sophisticated approach called "dpFeature".

# One effective way to isolate a set of ordering genes is to simply compare the cells collected at the beginning of the process to those at the end and find the differentially expressed genes, as described above. The command below will find all genes that are differentially expressed in response to the switch from growth medium to differentiation medium:


diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
              fullModelFormulaStr = "~Media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

















#Alternative2 = Selecting genesfrom singalling pathway

name_pathway <- "TipsVsStalk.TOP250"

#~ ---------------------------------------------------------------------------
#~ RELIC - commented out when this repository was assembled.
#~ This setwd() points into another project's output tree, which no longer
#~ exists on the cluster. setwd() to a missing directory aborts an R script,
#~ so the line is disabled rather than deleted; the working directory stays at
#~ this step's own output directory. See ../EDIT_POLICY.md.
#~ ---------------------------------------------------------------------------
#~ GSEA <- read.table(paste("/RQexec/gaelcge/Gene_list/",name_pathway,".csv",sep=""), row.names=1)

ordering_genes <- rownames(GSEA)

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
 
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
 
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}


ordering_genes <- convertMouseGeneList(ordering_genes)


ordering_genes <- toupper(ordering_genes)

                  
HSMM <- setOrderingFilter(HSMM, ordering_genes)


png("plot_ordering_genes.alt1.png")
plot_ordering_genes(HSMM)
dev.off()

HSMM <- reduceDimension(HSMM, max_components = 2,
    method = 'DDRTree')


HSMM <- orderCells(HSMM)

png("plot_cell_trajectory.alt1.png")
plot_cell_trajectory(HSMM, color_by = "State")
dev.off()

png("plot_cell_trajectory.alt1_Dataset.png")
plot_cell_trajectory(HSMM, color_by = "Dataset", show_tree = TRUE, show_backbone = TRUE, backbone_color = "black",
  markers = FALSE, use_color_gradient = FALSE, markers_linear = FALSE,
  show_cell_names = FALSE, show_state_number = FALSE, cell_size = 1.5,
  cell_link_size = 1, cell_name_size = 2, state_number_size = 2.9,
  show_branch_points = TRUE, theta = 0)
dev.off()


png("plot_cell_trajectory.alt1.wrap.png")
plot_cell_trajectory(HSMM, color_by = "State") +
    facet_wrap(~State, nrow = 1)
dev.off()

png("plot_spanning_tree.png")
plot_spanning_tree(HSMM) 
dev.off()



blast_genes <- row.names(subset(fData(HSMM),
gene_short_name %in% ordering_genes))

png("plot_genes_jitter.alt1_Dataset.png", res=150, height = 7000, width = 1000)
plot_genes_jitter(HSMM[blast_genes,],
    grouping = "Dataset",
  min_expr = 0.1)
dev.off() 




##Analyzing Branches in Single-Cell Trajectories

#BEAM takes as input a CellDataSet that's been ordered with orderCells and the name of a branch point in the trajectory. 
#It returns a table of significance scores for each gene. Genes that score significant are said to be branch-dependent in their expression.

BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

save(BEAM_res, file="BEAM_res.RData")

#You can visualize changes for all the genes that are significantly branch dependent using a special type of heatmap. 
#This heatmap shows changes in both lineages at the same time. It also requires that you choose a branch point to inspect. 
#Columns are points in pseudotime, rows are genes, and the beginning of pseudotime is in the middle of the heatmap. 
#As you read from the middle of the heatmap to the right, you are following one lineage through pseudotime. 
#As you read left, the other. The genes are clustered hierarchically, so you can visualize modules of genes that have similar lineage-dependent expression patterns.


png("plot_genes_branch_1_heatmap.5clusters.png", res = 300, height = 3000, width = 2000)
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                          pval < 0.01)),],
                                          branch_point = 1,
                                          num_clusters = 5,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T)
dev.off()

HSMM_genes_branch <- row.names(subset(fData(HSMM),
          gene_short_name %in% row.names(subset(BEAM_res,
                                          pval < 0.01))))
                                          
HSMM_genes_branch.fate1 <- row.names(subset(fData(HSMM),
          gene_short_name %in% c("PLTP", "PTMS", "BICD2", "BHLHB9", "PTP4A2", "TFP1", "MAP4K4", "PLAT", 
          "GM26924", "LPAR6")))

HSMM_genes_branch.fate2 <- row.names(subset(fData(HSMM),
          gene_short_name %in% c("MT-CYTB", "GM9803", "TOPBP1", "SLC2A1", "QSER1", "GSTM1", "SLC16A1", "ESAM", "IGFBP7",
          "GAS5", "CDIP1", "ANKRD11")))
          
HSMM_genes_branch.Pre_branched <- row.names(subset(fData(HSMM),
          gene_short_name %in% c("LY6A", "LY6C1", "CUL1", "CNEP1R1")))
          
          
png("plot_genes_branch_1_pseudotime.Dataset.p01.png", res = 300, height = 20000, width = 3000)
plot_genes_branched_pseudotime(HSMM[HSMM_genes_branch,],
                       branch_point = 1,
                       color_by = "Dataset",
                       ncol = 1)
dev.off()

png("plot_genes_branch_1_pseudotime.State.p01.png", res = 300, height = 20000, width = 3000)
plot_genes_branched_pseudotime(HSMM[HSMM_genes_branch,],
                       branch_point = 1,
                       color_by = "State",
                       ncol = 1)
dev.off()


png(filename="FeaturePlot_retina_CTL_OIR.TXNRD1.png", res = 300, width=2000, height=2000)
FeaturePlot(object = retina_CTL_OIR, features.plot = "TXNRD1", cols.use = c("lightgrey", 
    "blue"), pt.size = 0.5)
dev.off()


png(filename="FeaturePlot_cells_for_monocle.TXNRD1.png", res = 300, width=2000, height=2000)
FeaturePlot(object = cells_for_monocle, features.plot = "TXNRD1", cols.use = c("lightgrey", 
    "blue"), pt.size = 0.5)
dev.off()

png(filename="SplitDotPlotGG_retina_CTL_OIR.TXNRD1.png", width=2000, height=2000, bg = "transparent", res = 300)
SplitDotPlotGG(retina_CTL_OIR,
               grouping.var = "Dataset",
               genes.plot = "TXNRD1",
               plot.legend = TRUE,
               x.lab.rot = TRUE,
               dot.scale=5)
dev.off()


#END
q("no")
