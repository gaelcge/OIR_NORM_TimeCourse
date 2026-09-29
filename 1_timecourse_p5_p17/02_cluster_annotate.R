# ---------------------------------------------------------------------------
# Step 02 - Clusters, cell-type labels, and composition across the time course.
#
# ESTABLISHES
#     59 populated clusters collapsed to 21 cell types and 20 cell groups, the
#     cluster marker table, and the cell-type composition per condition and
#     timepoint. This is the annotation every later step inherits.
#
# WHY THIS WAY
#     Clusters are cut from the 2D UMAP embedding, not from the PC neighbour
#     graph: FindNeighbors(reduction = "umap", dims = 1:2). This is a deliberate
#     house choice, and it has a consequence that has to be reported with it -
#     the UMAP's own parameters and seed become clustering parameters. Seurat
#     defaults for n.neighbors, min.dist and seed.use were used, so a Seurat
#     version that changes one of those defaults would change the cluster
#     labels without any visible cause. Report this as clustering on the 2D
#     embedding; it is not clustering in expression space.
#
#     The trade-off, stated plainly: UMAP preserves local neighbourhoods but not
#     global distances, so clusters cut from the layout partly reflect the
#     embedding rather than expression alone. The reason it is preferred here is
#     that the embedding is what the analyst reads, so the clusters match what
#     is visible in the figure and can be named against retinal anatomy.
#
#     Resolution is 1. The cluster assignments live in `integrated_snn_res.1`.
#     The cluster-to-cell-type vector below has 60 entries because level 59
#     exists but is empty: cluster 59 is `Contaminant` and holds 0 cells, so it
#     contributes nothing to any composition figure.
#
#     Neuronal_progenitor_cells are grouped with Amacrine_cells downstream. That
#     is a substantive call, not a typo - see METHODS_CHOICES.md.
#
# REQUIRES
#     step 01 (reads its integrated object)
#
# INPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#     functions/seurat_to_loom.R                sourced, not a step
#
# OUTPUTS
#     <OBJ_DIR>/<project_name>.<res>.<Dim>.<perp>.Seurat_object.integrated.rds
#                                               re-saved with clusters + labels
#     <OBJ_DIR>/Seurat_object.integrated_P12_14_17.rds
#     <OBJ_DIR>/{P12,P14,P17}.Seurat_object.integrated.rds
#     Seurat_object.cluster.markers.txt         FindAllMarkers, positive only
#     Tsne/Umap/DotPlot/Vln figures             see the README step table
#     *.Seurat.loom_object/                     loom exports for python tools
#
# USAGE
#     Rscript 02_cluster_annotate.R
#     Sections run top to bottom; section 0 is the parameter block.
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

# object identity, interpolated into every output filename by the original
# scripts. res=1 is the clustering resolution actually used - confirmed from the
# object, whose clustering column is `integrated_snn_res.1`.
project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"
res  <- 1
Dim  <- 20
DIM_nb <- c(1:20)
perp <- 30

# clustering, on the 2D UMAP embedding
K_PARAM     <- 20
PRUNE_SNN   <- 1/15
CLUSTER_RES <- 1        # -> integrated_snn_res.1
RANDOM_SEED <- 0

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


Seurat_object.integrated <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))


#Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset. Commands are identical to the standard workflow, but do not run the ScaleData function after integration. You can see that after integration, cells group by their biological cell type (which has been pre-annotated), instead of by their underlying technology.

Seurat_object.integrated <- RunPCA(Seurat_object.integrated, verbose = FALSE)

Seurat_object.integrated <- RunUMAP(Seurat_object.integrated, dims = DIM_nb, reduction.name = "umap")

plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "umap", group.by = c("Cond_Sorting", "ident", "TimePoint"),
  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
    byrow = TRUE, override.aes = list(size = 3))))

png(filename=paste("UmapPlot", project_name, res, Dim, perp, "ident.cca.png", sep="."), width=1500, height=1500, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


Seurat_object.integrated <- RunTSNE(Seurat_object.integrated, reduction = "pca", cells = NULL, assay = "integrated", 
  dims = DIM_nb, features = NULL, seed.use = 1, tsne.method = "Rtsne",
  add.iter = 0, dim.embed = 2, distance.matrix = NULL,
  reduction.name = "tsne", reduction.key = "tSNE_")

plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "tsne", group.by = c("Cond_Sorting", "ident", "TimePoint"),
  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
    byrow = TRUE, override.aes = list(size = 3))))

png(filename=paste("TsnePlot", project_name, res, Dim, perp, "ident.cca.png", sep="."), width=1500, height=1500, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


##Find cluster

png(filename=paste("DimHeatmap", project_name, res, Dim, perp, "dim.1-12.png", sep="."), width=1000, height=2000, bg = "white", res = 150)
DimHeatmap(Seurat_object.integrated, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

png(filename=paste("ElbowPlot", project_name, res, Dim, perp, "dim.1-20.png", sep="."), width=1000, height=1500, bg = "white", res = 150)
ElbowPlot(Seurat_object.integrated, ndims = 30, reduction = "pca")
dev.off()


Seurat_object.integrated <- FindNeighbors(Seurat_object.integrated, reduction = "umap", dims = c(1:2),
  assay = NULL, features = NULL, k.param = 20, compute.SNN = TRUE,
  prune.SNN = 1/15, nn.eps = 0, verbose = TRUE,
  force.recalc = TRUE, do.plot = FALSE, graph.name = NULL)



Seurat_object.integrated <- FindClusters(Seurat_object.integrated, graph.name = NULL,
  modularity.fxn = 1, initial.membership = NULL, weights = NULL,
  node.sizes = NULL, resolution = res, algorithm = 1, n.start = 10,
  n.iter = 10, random.seed = 0, group.singletons = TRUE,
  temp.file.location = getwd(), edge.file.name = NULL, verbose = TRUE)


# Print cluster ident


plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "tsne", group.by = c("ident"),
  split.by = NULL, shape.by = NULL, order = NULL, label = TRUE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))+ NoLegend())

png(filename=paste("TsnePlot", project_name, res, Dim, perp, "cluster_ident.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "umap", group.by = c("ident"),
  split.by = NULL, shape.by = NULL, order = NULL, label = TRUE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))+ NoLegend())

png(filename=paste("umapPlot", project_name, res, Dim, perp, "cluster_ident.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


#Map other paramaters = Cond_TimePoint

plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "tsne", group.by = c("Cond_TimePoint"),
  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))+ theme(legend.position = "top"))

png(filename=paste("TsnePlot", project_name, res, Dim, perp, "Cond_TimePoint.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "umap", group.by = c("Cond_TimePoint"),
  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))+ theme(legend.position = "top"))

png(filename=paste("umapPlot", project_name, res, Dim, perp, "Cond_TimePoint.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


#Map other paramaters = TimePoint

plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "tsne", group.by = c("TimePoint"),
  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))+ theme(legend.position = "top"))

png(filename=paste("TsnePlot", project_name, res, Dim, perp, "TimePoint.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


plots <- DimPlot(Seurat_object.integrated, cells = NULL, cols = NULL,
  pt.size = NULL, reduction = "umap", group.by = c("TimePoint"),
  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
  label.size = 4, repel = FALSE, cells.highlight = NULL,
  cols.highlight = "#DE2D26", sizes.highlight = 1,
  na.value = "grey50", combine = FALSE, ncol = NULL)
plots <- lapply(X = plots, FUN = function(x) x + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))+ theme(legend.position = "top"))

png(filename=paste("umapPlot", project_name, res, Dim, perp, "TimePoint.png", sep="."), width=1000, height=1000, bg = "white", res = 150)
CombinePlots(plots)
dev.off()


###Check QC parameters

png(filename=paste("FeaturePlot", project_name, res, Dim, perp, "nFeature.png", sep="."), width=2000, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object.integrated, features=c("nFeature_RNA", "nCount_RNA"), dims = c(1, 2), cells = NULL,
  cols = c("yellow", "red"), pt.size = NULL, order = FALSE,
  min.cutoff = NA, max.cutoff = NA, reduction = NULL,
  split.by = NULL, shape.by = NULL, slot = "data", blend = FALSE,
  blend.threshold = 0.5, label = FALSE, label.size = 4,
  repel = FALSE, ncol = NULL, combine = TRUE, coord.fixed = FALSE,
  by.col = TRUE)
dev.off()


#
png(filename=paste("VlnPlotQC", project_name, res, Dim, perp, "Cond_Sorting.regressed.png", sep="."), width=700, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object.integrated, features=c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Cond_Sorting", ncol = 2)
dev.off()

png(filename=paste("VlnPlotQC", project_name, res, Dim, perp, "Cluster.ident.regressed.png", sep="."), width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object.integrated, features=c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = NULL, ncol = 2)
dev.off()


####Check features

DefaultAssay(Seurat_object.integrated) <- "RNA"

Seurat_object.integrated <- NormalizeData(Seurat_object.integrated, verbose = FALSE)
Seurat_object.integrated <- ScaleData(Seurat_object.integrated, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch", "S.Score", "G2M.Score"))

markers.retina.dotplot <- rev(c("RHO", "LHX1", "SLC17A6", "PAX6", "GAD1", "SLC5A7", "SLC6A9", "OPN1MW", "VSX2", "OTX2", 
                                "PRDM1", "RLBP1", "GFAP", "IGFBP5", "NDUFA4L2", "PECAM1", "KCNJ8", "CX3CR1", "VIM", "FBN1", "C1QA", "NEFL"))


png(filename=paste(project_name,res, Dim, perp, "DotPlot_markers.cca.png", sep="."), res = 150, width=2000, height=3000)
DotPlot(Seurat_object.integrated, features= markers.retina.dotplot, cols = c("blue", "red"),
  col.min = NA, col.max = NA, dot.min = 0, dot.scale = 8) + RotatedAxis()
dev.off()




### Assigning annotation
Seurat_object.integrated <- SetIdent(Seurat_object.integrated, cells = NULL, value="integrated_snn_res.1")

cell_type_assigned <- c("Bipolar_cells_1",      # cluster 0
                      "Rods",     # cluster 1
                      "Rods",     # cluster 2
                      "Rods",     # cluster 3
                      "Bipolar_cells_2",     # cluster 4
                      "Rods",     # cluster 5
                      "Glycinergic_amacrine_cells",     # cluster 6
                      "Cones",     # cluster 7
                      "Amacrine_cells",     # cluster 8
                      "Bipolar_cells_2",     # cluster 9
                      "Rods",     # cluster 10
                      "Muller_glia",     # cluster 11
                      "Bipolar_cells_1",     # cluster 12
                      "Bipolar_cells_2",     # cluster 13
                      "Muller_glia",     # cluster 14
                      "Cones",    # cluster 15
                      "Bipolar_cells_1",     # cluster 16
                      "Rods",     # cluster 17
                      "Early_rods",     # cluster 18
                      "Rods",     # cluster  19
                      "Rods",     # cluster 20
                      "Bipolar_cells_2",     # cluster 21
                      "Neuronal_progenitor_cells",     # cluster 22
                      "Glycinergic_amacrine_cells",     # cluster 23
                      "Bipolar_cells_2",     # cluster 24
                      "Bipolar_cells_1",     # cluster 25
                      "Rods",     # cluster 26
                      "Rods",     # cluster 27
                      "Early_muller_glia",     # cluster 28
                      "Rods",     # cluster 29
                      "Bipolar_cells_1",     # cluster 30
                      "Bipolar_cells_2",     # cluster 31
                      "Bipolar_cells_1",     # cluster 32
                      "Muller_glia",     # cluster 33
                      "Rods",     # cluster 34
                      "Rods",     # cluster 35
                      "Rods",     # cluster 36
                      "Rods",     # cluster 37
                      "Muller_glia",     # cluster 38
                      "Retinal_ganglion_cells",     # cluster 39
                      "Amacrine_cells",     # cluster 40
                      "Early_bipolar_cells",     # cluster 41
                      "Rods",     # cluster 42
                      "Muller_glia",     # cluster 43
                      "Immune_cells",     # cluster 44
                      "Early_amacrine_cells",     # cluster 45
                      "Bipolar_cells_1",     # cluster 46
                      "Pericytes",     # cluster 47
                      "Opticin_cells",     # cluster 48
                      "Rods",     # cluster 49
                      "Amacrine_cells",     # cluster 50
                      "Endothelial_cells",     # cluster 51
                      "Activated_muller_glia",     # cluster 52
                      "Early_muller_glia",     # cluster 53
                      "Rods",     # cluster 54
                      "Horizontal_cells",     # cluster 55
                      "Bipolar_cells_1",     # cluster 56
                      "Cholinergic_amacrine_cells",     # cluster 57
                      "Astrocytes",     # cluster 58
                      "Contaminant"     # cluster 59
                                   )    



names(cell_type_assigned) <- levels(Seurat_object.integrated)
Seurat_object.integrated <- RenameIdents(Seurat_object.integrated, cell_type_assigned)

Seurat_object.integrated[["Cell_Type"]] <- Idents(object = Seurat_object.integrated)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

png(filename="tsne.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


Seurat_object.integrated <- RenameIdents(Seurat_object.integrated, "Activated_muller_glia" = "Muller_glia")


Seurat_object.integrated[["Cell_Group"]] <- Idents(object = Seurat_object.integrated)


###FInd markers for non-annotated clusters

Seurat_object.markers <- FindAllMarkers(Seurat_object.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Seurat_object.markers, "Seurat_object.cluster.markers.txt")

Seurat_object.markers <- read.table("Seurat_object.cluster.markers.txt")

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

png(filename="VlnPlot_top2Markers.sct.png", width=3500, height=4500, bg = "white", res = 50)
VlnPlot(Seurat_object.integrated, features = unique(top2$gene))
dev.off()

png(filename="FeaturePlot_top2Markers.sct.png", width=3000, height=4000, bg = "white", res = 150)
FeaturePlot(Seurat_object.integrated, features = c(unique(top2$gene)))
dev.off()

png(filename="DoHeatmap_top10Markers.sct.png", width=2000, height=3000, bg = "white", res = 150)
DoHeatmap(Seurat_object.integrated, features = top10$gene) + NoLegend()
dev.off()

png(filename="DoHeatmap_top2Markers.sct.png", width=900, height=700, bg = "white", res = 150)
DoHeatmap(Seurat_object.integrated, features = top2$gene) + NoLegend()
dev.off()


png(filename="DotPlot_markers.top2.sct.png", res = 150, width=2000, height=2000)
DotPlot(
  Seurat_object.integrated,
  assay = NULL,
  unique(top2$gene),
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
dev.off()


##REmove contaminants

Seurat_object.integrated <- SubsetData(Seurat_object.integrated, ident.remove = "Contaminant")

png(filename="umap.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

png(filename="tsne.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


colorsplot <- DiscretePalette(21, palette = "polychrome")

png(filename="umap.sct_filtered.annotated.polychrome.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, cols = colorsplot) + NoLegend()
dev.off()

png(filename="tsne.sct_filtered.annotated.polychrome.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = TRUE, pt.size = 0.5, cols = colorsplot) + NoLegend()
dev.off()


png(filename="umap.sct_filtered.annotated.polychrome.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()

png(filename="tsne.sct_filtered.annotated.polychrome.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()

colorsplot <- DiscretePalette(21, palette = "stepped")

png(filename="umap.sct_filtered.annotated.stepped.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()

png(filename="tsne.sct_filtered.annotated.stepped.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()


colorsplot <- DiscretePalette(21, palette = "glasbey")

png(filename="umap.sct_filtered.annotated.glasbey.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()

png(filename="tsne.sct_filtered.annotated.glasbey.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()


setEPS(bg = "white", family = "Times", width=8, height=6)
postscript("umap.sct_filtered.annotated.alphabet.legend.eps")
DimPlot(Seurat_object.integrated, reduction = "umap", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()




colorsplot <- DiscretePalette(21, palette = "alphabet2")

png(filename="umap.sct_filtered.annotated.alphabet2.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()

png(filename="tsne.sct_filtered.annotated.alphabet2.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()


colorsplot <- DiscretePalette(21, palette = "alphabet")

png(filename="umap.sct_filtered.annotated.alphabet.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "umap", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()

png(filename="tsne.sct_filtered.annotated.alphabet.legend.png", width=2000, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object.integrated, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colorsplot)
dev.off()



###Make proportion graph

phase_counts <- Seurat_object.integrated@meta.data %>%
  group_by(Cond_TimePoint) %>%
  count(Cell_Type) %>%
  mutate(freq = n / sum(n)*100)  


png(filename=paste("Cell_Type_proportions_Cond_TimePoint.png", sep="."), width=4000, height=3000, res = 300) 
ggplot(phase_counts,aes(x=reorder(Cell_Type,-freq),y=freq,fill=Cell_Type)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cond_TimePoint, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

ggplot(phase_counts,aes(x=Cond_TimePoint,y=freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_Type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette="YlGnBu", direction=1)
ggsave(filename=paste("Cell_Type_proportions_Condition_perNew_Cell_SubType.YlGnBu.eps", sep="."), width=8, height=8) 



phase_counts <- phase_counts %>% group_by(Cond_TimePoint) %>% mutate(freq_sum = sum(freq))

png(filename=paste("Cell_Type_proportions_Cond_TimePoint.stacked.png", sep="."), width=2000, height=1500, res=300)
ggplot(phase_counts,aes(x=Cond_TimePoint,y=freq,fill=Cell_Type)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cond_TimePoint') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



###save annotated Seurat object

Seurat_object.integrated[["Cell_Type"]]$Cell_Type <- gsub("_", "", Seurat_object.integrated[["Cell_Type"]]$Cell_Type)

saveRDS(Seurat_object.integrated, paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))

Seurat_object.integrated <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))

##Convert to loon object

source("SeuratToLoom.R")

seuratToLoom(obj=paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."),dir= paste("NORM_OIR_time_course", "Seurat.loom_object", sep="."))



##Subet p12 to P17 

Seurat_object.integrated <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))

Seurat_object.integrated_P12_14_17 <- SetIdent(Seurat_object.integrated, value = "TimePoint")

Seurat_object.integrated_P12_14_17 <- SubsetData(Seurat_object.integrated_P12_14_17, ident.use = c("P12", "P14", "P17"))


phase_counts <- Seurat_object.integrated_P12_14_17@meta.data %>%
  group_by(Cond_TimePoint) %>%
  count(Cell_Type) %>%
  mutate(freq = n / sum(n)*100)  


png(filename=paste("Cell_Type_proportions_P12_14_17_Cond_TimePoint.png", sep="."), width=4000, height=3000, res = 300) 
ggplot(phase_counts,aes(x=reorder(Cell_Type,-freq),y=freq,fill=Cell_Type)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cond_TimePoint, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

ggplot(phase_counts,aes(x=Cond_TimePoint,y=freq,fill=Cond_TimePoint)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_Type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette="YlGnBu", direction=1)
ggsave(filename=paste("Cell_Type_proportions_P12_14_17_Condition_perNew_Cell_SubType.YlGnBu.eps", sep="."), width=8, height=8) 



phase_counts <- phase_counts %>% group_by(Cond_TimePoint) %>% mutate(freq_sum = sum(freq))

png(filename=paste("Cell_Type_proportions_P12_14_17_Cond_TimePoint.stacked.png", sep="."), width=2000, height=1500, res=300)
ggplot(phase_counts,aes(x=Cond_TimePoint,y=freq,fill=Cell_Type)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Cond_TimePoint') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




saveRDS(Seurat_object.integrated_P12_14_17, paste("Seurat_object.integrated_P12_14_17.rds", sep="."))

##Convert to loon object

source("SeuratToLoom.R")

seuratToLoom(obj=paste("Seurat_object.integrated_P12_14_17.rds", sep="."),dir= paste("NORM_OIR_time_course_P12_14_17", "Seurat.loom_object", sep="."))


##Subet p12 P14 P17 separately 

Seurat_object.integrated <- readRDS(paste(project_name, res, Dim, perp, "Seurat_object.integrated.rds", sep="."))

for (i in c("P12", "P14", "P17")) {

Seurat_object.integrated_Psub <- SetIdent(Seurat_object.integrated, value = "TimePoint")

Seurat_object.integrated_Psub <- SubsetData(Seurat_object.integrated_Psub, ident.use = i)


saveRDS(Seurat_object.integrated_Psub, paste(i,"Seurat_object.integrated.rds", sep="."))

##Convert to loon object

source("SeuratToLoom.R")

seuratToLoom(obj=paste(i,"Seurat_object.integrated.rds", sep="."),dir= paste("NORM_OIR",i, "Seurat.loom_object", sep="."))

}




q("no")
