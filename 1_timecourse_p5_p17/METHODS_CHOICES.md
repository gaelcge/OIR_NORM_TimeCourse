# Methods choices — time course, P5 to P17

Choices in this folder that a reader could reasonably mistake for an error, and
the reason each one is not. Numbered so they can be cited from a reviewer
response. Where a choice departs from what its source does, the source is
cited here; where it follows the source, the citation is in the README beside
the step instead, so no citation appears twice.

---

## 1. Clusters are cut from the 2D UMAP embedding, not from the PC graph

`02_cluster_annotate.R` builds the neighbour graph with
`FindNeighbors(reduction = "umap", dims = 1:2, k.param = 20, prune.SNN = 1/15)`
and then runs `FindClusters(resolution = 1, algorithm = 1, random.seed = 0)`.
The standard Seurat workflow builds that graph on the first 20–50 principal
components; this does not.

**Why.** The embedding is what the analyst reads. Clusters cut from the 2D
layout correspond to the groups visible in the figure, which makes them
possible to name against retinal anatomy — 59 clusters resolving into rods,
cone subsets, bipolar subsets, amacrine subsets, Müller glia, astrocytes,
pericytes, endothelium and immune cells.

**What it costs, stated so it is not discovered by a reviewer.** UMAP preserves
local neighbourhoods but not global distances or densities, so clusters cut
from the layout partly reflect the embedding's parameters rather than
expression alone, and the same data re-embedded with a different seed can split
or merge groups. Concretely: `RunUMAP` was called with Seurat's defaults for
`n.neighbors`, `min.dist` and `seed.use`, and those defaults are therefore
clustering parameters here. A Seurat version that changed one of them would
change the cluster labels with no visible cause.

**How to report it.** "Clusters were defined by Louvain community detection on a
shared-nearest-neighbour graph of the 2D UMAP embedding (k = 20, prune.SNN =
1/15, resolution = 1, seed = 0)." Not "clustering on expression space", and not
"clustering on principal components".

**What is not available.** The PCA-graph clustering was not computed alongside,
so the object carries only `integrated_snn_res.1` and there is no
same-object comparison to show a reviewer who asks for one.

Departure from: Butler A et al. Integrating single-cell transcriptomic data
across different conditions, technologies, and species. Nat Biotechnol
2018;36:411–420. doi:10.1038/nbt.4096 — whose workflow clusters on principal
components.

---

## 2. The clustering resolution is 1, and the script used to say otherwise

The parameter block sets `res = 1`; the cluster original additionally assigned
`res = 0.8` on the line immediately before `FindClusters`. That second
assignment did not take effect in the run that produced the object, established
two independent ways:

- the object's clustering column is named `integrated_snn_res.1`
- every output filename on the cluster carries `1.20.30`, and a run at 0.8
  would have written `0.8.20.30`

The line is deleted in this repository, because a published script containing it
would tell a reader the analysis used 0.8. This is one of the five permitted
edits in `EDIT_POLICY.md`.

---

## 3. `percent.crystal < 0.025` — an unusual QC filter

Standard scRNA-seq QC filters on library size, gene count and mitochondrial
fraction. `01_integrate_by_sorting.R` adds a crystallin fraction, summing CRYA
and CRYB genes, and removes cells above 2.5%.

**Why.** Dissected retina carries lens fragments. Crystallin-high droplets are
lens contamination, not a retinal population, and left in they form clusters
that would be annotated as a cell type.

**Known dead line.** `CRYG.genes` is computed on the line above and never used;
only CRYA and CRYB enter the filter. Left in place — the script produced
published numbers and the line changes nothing.

---

## 4. Cell-cycle scores are regressed inside SCTransform, not corrected after

`vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch", "S.Score",
"G2M.Score")` is passed to `SCTransform` per sorted fraction, before
integration anchors are found.

**Why.** The retina between P5 and P17 is still proliferating, so cell cycle is
real biological signal that would otherwise dominate the early timepoints and
drive the anchors. Regressing inside the normalisation means the anchors are
computed on residuals already free of cycle structure. Correcting afterwards
would leave the integration itself cycle-driven.

**The cost.** Proliferation is deliberately removed as an axis of variation, so
this object cannot be used to study proliferation across the time course. A
question about progenitor cycling needs the un-regressed object.

---

## 5. `Neuronal_progenitor_cells` are grouped with `Amacrine_cells`

`annotation/celltype_grouping.csv` maps `Neuronal_progenitor_cells` →
`Amacrine_cells` for the differential-expression steps. This is a scientific
call, not a typo, and it is the entry in this file most likely to be
challenged: a progenitor population and a differentiated interneuron population
are being pooled.

The grouping is applied only in steps 04–06, which need cell-type groups with
enough cells per condition to test. The fine label survives untouched in the
object's `Cell_Type` column (982 cells), so any analysis that needs the
progenitors separate can recover them.

**No independent precedent is cited for this grouping.** It was chosen for this
dataset. A reader who disagrees can re-run steps 04–06 against `Cell_Type`
instead of `General_CellType` and the ten-row CSV is there to make the
difference explicit.

---

## 6. Differential-expression tables are filtered on the unadjusted p-value

Every table `04_de_genes.R` writes is produced by
`subset(DGE_test, p_val < 0.05)` and named `..._p0.05.txt`. The volcano plots in
the same script colour points on `p_val_adj < 0.05`. So the tables and the
figures in that step do not use the same threshold, and the filename does not
say which one it means.

**Correction family and direction.** Seurat's `FindMarkers` returns `p_val_adj`
as Bonferroni across all genes in the assay, which is conservative; the
unadjusted `p_val` is what the tables were cut on. Both columns are present in
every written table, so a reader can re-filter.

**Not reconciled.** These tables accompany published results. Re-filtering them
would change what the repository reports relative to the paper. Stated here
instead, so that `p0.05` in a filename is not read as FDR 0.05.

---

## 7. Human gene sets are applied to mouse data by symbol identity

MSigDB v7.1 human symbol collections are used directly against the mouse
expression matrix. The merged count matrix carries uppercased mouse symbols
(verified: `MT-ATP6`, `0610005C13RIK`), so the join succeeds on string
identity rather than on any orthology mapping.

**The limitation.** Mouse genes whose symbol differs from the human orthologue
do not match and are silently absent from every gene set. No orthology mapping
(biomaRt, MGI homology) was applied. Gene-set results should therefore be read
as enrichment over the intersection of each set with the uppercased mouse
symbol space, not over the full human set.

---

## 8. GSVA is not reproducible from this repository

`03_gsva_senescence.R` combines the MSigDB collections with three curated
`.gmx` files — `MALLETTE_SENESCENCE_UP`,
`GLOBAL_SENESCENCE_LITERATURE_CURATED`, `SAWCHYN_UNBIASED_UP` — that live in a
collaborator's project directory and are not readable by this account
(permission denied, checked 2026-09-29). They are not committed and their
contents are recorded nowhere here.

Consequence: step 03 cannot be re-run as published, and neither can
`05_de_functions.R`, which loads step 03's output as a `"GO"` assay. The MSigDB
half of the collection is pinned by checksum in `reference/PROVENANCE.txt`; the
three curated sets are not specified at all. A re-run with MSigDB alone will
produce a different score matrix.

`kcdf = "Gaussian"` is correct for the input and is not a departure — GSVA is
given log-normalised expression rather than counts.

---

## 9. `GO_HEXOSE_CATABOLIC_PROCESS.gmt` is empty

The file `04_de_genes.R` reads for that term is 0 bytes of content on the
cluster (sha256 in `reference/PROVENANCE.txt`). The script later indexes
`gsc_final[["GO_HEXOSE_CATABOLIC_PROCESS"]]` while assembling a glycolysis /
pentose-phosphate gene list, so that term contributed no genes.

Whether the affected section ran at all cannot be determined from the files, and
it was not re-run to find out. Left in the script and recorded here.

---

## 10. 244 lines were deleted from the fgsea step

Source lines 353–596 of `fgsea_working.R` were a second analysis contrasting
`ident.1 = "VldlrKO"` against `ident.2 = "CTL"` per cell type. Removed, on this
evidence:

- those identities do not exist in this object's metadata
- no output of that section is present anywhere in the fgsea directory — no
  `*Vldlr*` file, no `barplotNESPvalue.PR.*`

So it cannot have run on these cells, and a published script containing it would
invite the question of which results came from it.

**Two related things were deliberately NOT changed.** The gene set named
`JOYAL_FAO_VldlrKO_Retina` (around line 193) is a fatty-acid-oxidation
signature derived from a Vldlr-KO retina study, used here as one of the pathways
tested in OIR — the name refers to the signature's provenance, not to a
contrast, and it stays. And one `write.table` inside the retained OIR section
names its output `tgc2<GOI>_VldlrKOvsCTL.celltype_NES_summary.txt`, mislabelling
an OIR-vs-NORM result; no file of that name exists on the cluster, so no
published result carries it, and the line is left as it ran.

---

## 11. The fgsea pathway was swept by editing the script

`GOI` is assigned a single pathway name in the body of `06_fgsea.R`, yet four
pathway output directories exist on the cluster. The step was run four times
with that one line edited between runs. The four values are recorded as
`GOI_SWEPT` in the parameter block so the repository does not imply a single
pathway was examined; the body is unchanged and still uses its own single
assignment.

**What survives of that step.** Two `.eps` barplots per pathway, eight files in
total. The `GSEA_<GOI>.txt`, `tgc2*` and `Anova_Tukey_*` tables the script
writes are not present in any pathway directory.

---

## 12. Step 03b exists and feeds nothing

`03b_gsva_oir_only.R` sat in a directory named `Cpu`, which suggests a
performance variant. It is not one — it covers OIR cells only, uses MSigDB
without the curated senescence sets, and writes `gsva.exprs.MsigDB.csv`, which
no other script reads. Kept and numbered so that a reader comparing this
repository against the cluster directory does not find an unexplained output.
