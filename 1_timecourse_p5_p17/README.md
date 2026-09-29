# 1_timecourse_p5_p17 — the P5→P17 time course

**What this layer establishes:** how retinal cell-type composition changes
across postnatal development and what OIR does to it, from a single integrated
object of **51,151 cells × 21,705 features** spanning P5, P7, P10, P12, P14 and
P17; then which genes and which gene-set programmes differ between OIR and
normoxia within each cell type.

**The object:** `RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting.1.20.30.Seurat_object.integrated.rds`
— integrated across the two sorted fractions by SCT anchors, clustered at
resolution 1 on the 2D UMAP embedding into 59 populated clusters, annotated to
21 cell types and 20 cell groups. Seurat object version 3.1.5.

> This is **not** the object used by `2_compartments_p14_p17`. Cluster numbers
> and label sets do not transfer between the two. Read item 1 of
> `METHODS_CHOICES.md` before combining any output from the two folders.

Remember the design constraint from the root README: **OIR cells exist only at
P12, P14 and P17.** Every condition contrast below is restricted accordingly.

## Steps

| Step | Script | Consumes | Produces | Establishes |
|---|---|---|---|---|
| 01 | `01_integrate_by_sorting.R` | merged matrix; Regev cell-cycle list | `…Seurat_object.integrated.rds`, `DGE.data.sparse.….rds`, initial QC violins | the integrated object; that the two sorted fractions can be reconciled without regressing out condition |
| 02 | `02_cluster_annotate.R` | step 01 | the object re-saved with clusters and labels; `Seurat_object.cluster.markers.txt`; t-SNE/UMAP/DotPlot figures; P12–P17 and per-timepoint subsets; loom exports | 59 clusters → 21 cell types, and the composition shift across condition and timepoint |
| 03 | `03_gsva_senescence.R` | step 02; MSigDB v7.1; three unavailable curated sets | `gsva.exprs.MsigDB_h_c2_c5.senescence.all.v7.1.symbols.csv` | per-cell gene-set scores — **not reproducible, see below** |
| 03b | `03b_gsva_oir_only.R` | step 02; MSigDB v7.1 | `gsva.exprs.MsigDB.csv` | nothing downstream; kept so the cluster directory has no unexplained output |
| 04 | `04_de_genes.R` | step 02 | per-cell-type DE tables, volcano plots | the OIR-vs-normoxia gene-level response at P14, at P17, and pooled |
| 05 | `05_de_functions.R` | **step 03**, step 02 | gene-set DE tables, Nebulosa density plots, clustered heatmaps | which gene-set programmes, senescence among them, differ by condition within cell type |
| 06 | `06_fgsea.R` | step 02; MSigDB v7.1; curated pathway lists | 8 NES barplots across 4 pathways | pathway-level enrichment from unfiltered rankings |
| — | `functions/seurat_to_loom.R` | — | — | a helper sourced by step 02; not a step |

**Run order note.** GSVA is numbered **03**, before the differential-expression
steps, because `05_de_functions.R` reads its output file and loads it onto the
object as a `"GO"` assay. That dependency is the reason for the numbering and is
easy to miss from the original directory layout, where GSVA sat in a sibling
folder to `DifferentialExpression`.

`01` and `02` are one logical step split across two files by the 8.5 h
walltime — `01` builds and saves the object, `02` reads it back. The Slurm
wrapper originally ran both in sequence.

---

## 01 — `01_integrate_by_sorting.R`

QC filter (`nFeature_RNA` 100–6000, `nCount_RNA` < 10000, `percent.mito` < 0.10,
`percent.crystal` < 0.025) → cell-cycle scoring → per-fraction `SCTransform`
regressing `nFeature_RNA`, `percent.mito`, `Batch`, `S.Score`, `G2M.Score` →
`SelectIntegrationFeatures(nfeatures = 3000)` → `PrepSCTIntegration` →
`FindIntegrationAnchors(normalization.method = "SCT", dims = 1:20)` →
`IntegrateData`.

Runtime ~7 h on 16 cores / 180 GB.

**Citations.**
- **Anchor-based integration**, the method this step rests on — Stuart T et al.
  Comprehensive integration of single-cell data. *Cell* 2019;177:1888–1902.
  doi:[10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)
- **SCTransform**, the normalisation the anchors are computed on, and the
  mechanism by which the confounders above are regressed out — Hafemeister C,
  Satija R. *Genome Biol* 2019;20:296.
  doi:[10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1)
- **The cell-cycle gene lists** — Tirosh I et al. *Science* 2016;352:189–196.
  doi:[10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501).
  The cached copy is identified by sha256 in `../reference/PROVENANCE.txt`; the
  script splits it positionally at lines 43/44, which is correct for that copy
  and only for that copy.
- The **crystallin QC filter** has `no precedent; chosen because dissected
  retina carries lens fragments and crystallin-high droplets cluster as a
  spurious cell type`. See `METHODS_CHOICES.md` item 3.

## 02 — `02_cluster_annotate.R`

`RunPCA` → `RunUMAP(dims = 1:20)` → `RunTSNE(perplexity = 30)` →
`FindNeighbors(reduction = "umap", dims = 1:2, k.param = 20, prune.SNN = 1/15)`
→ `FindClusters(resolution = 1, algorithm = 1, random.seed = 0)` →
`FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)` → hand
annotation → composition figures → P12/P14/P17 subsets and loom exports.

The cluster→cell-type and cell-type→group maps are committed beside the script
as `annotation/cluster_annotation.csv` (60 rows) and
`annotation/celltype_grouping.csv` (10 rows). The second is identical in all
three scripts that apply it — checked pair by pair — so it is safe to read as
the single definition.

**Citations.**
- **Louvain community detection**, the clustering algorithm
  (`algorithm = 1`) — Blondel VD et al. *J Stat Mech* 2008;P10008.
  doi:[10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)
- **UMAP**, which here supplies both the embedding and the space the neighbour
  graph is built in — McInnes L, Healy J, Melville J. arXiv:1802.03426.
  *Preprint; no peer-reviewed version.*
- **t-SNE** — van der Maaten L, Hinton G. *J Mach Learn Res* 2008;9:2579–2605.
- **SingleR** is loaded and informs the labelling, but the committed labels are
  assigned by hand from marker expression — Aran D et al. *Nat Immunol*
  2019;20:163–172.
  doi:[10.1038/s41590-018-0276-y](https://doi.org/10.1038/s41590-018-0276-y)
- **Clustering on the 2D UMAP embedding rather than on principal components has
  no precedent cited**; it is a deliberate departure from the Seurat workflow
  above, chosen because the clusters then correspond to what is visible in the
  figure and can be named against retinal anatomy. Its cost is stated in
  `METHODS_CHOICES.md` item 1, which is where the departure is documented.

## 03 / 03b — GSVA

`gsva(method = "gsva", kcdf = "Gaussian")` over MSigDB v7.1 h + c2 + c5,
combined in step 03 with three curated senescence sets.

> **Step 03 is not reproducible from this repository.** The three curated
> `.gmx` files are held by a collaborator and are not readable (checked
> 2026-09-29). Their contents are recorded nowhere here. A re-run with MSigDB
> alone yields a different score matrix, and step 05 inherits the problem. Full
> statement in `../reference/PROVENANCE.txt` §3.

**Citations.**
- **GSVA** — Hänzelmann S, Castelo R, Guinney J. *BMC Bioinformatics* 2013;14:7.
  doi:[10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)
  `kcdf = "Gaussian"` follows this paper's guidance for continuous,
  log-normalised input; it is not a departure.
- **MSigDB v7.1** — Liberzon A et al. *Cell Syst* 2015;1:417–425.
  doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004);
  Subramanian A et al. *PNAS* 2005;102:15545–15550.
  doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).
  Release and checksums in `../reference/PROVENANCE.txt`.

## 04 — `04_de_genes.R`

`FindMarkers` per cell type, OIR vs normoxia at P14, at P17, and pooled.

> **The written tables are filtered on the unadjusted p-value** while the volcano
> plots colour on the adjusted one, and the filenames read `p0.05` either way.
> The adjusted column is in every table. `METHODS_CHOICES.md` item 6.

**Citations.**
- **Wilcoxon rank-sum test** as implemented in `FindMarkers` — Seurat, Butler A
  et al. doi:[10.1038/nbt.4096](https://doi.org/10.1038/nbt.4096)
- **Why a rank-based test rather than a parametric or zero-inflated one**, from
  an independent benchmark rather than the tool's own paper — Soneson C,
  Robinson MD. Bias, robustness and scalability in single-cell differential
  expression analysis. *Nat Methods* 2018;15:255–261.
  doi:[10.1038/nmeth.4612](https://doi.org/10.1038/nmeth.4612)
- **Multiple-testing correction**: Seurat's `p_val_adj` is Bonferroni across
  assay genes. Where FDR is discussed instead — Benjamini Y, Hochberg Y. *J R
  Stat Soc B* 1995;57:289–300.
  doi:[10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- Testing **within cell type rather than across the whole retina** has `no
  precedent cited; chosen because OIR shifts cell-type composition (step 02), so
  a whole-retina contrast cannot separate a change in proportion from a change
  in expression`.

## 05 — `05_de_functions.R`

Loads the step 03 GSVA matrix onto the object as a `"GO"` assay, CLR-normalises
it, and runs the same `FindMarkers` contrasts on gene sets that step 04 runs on
genes.

**Citations.**
- **GSVA** as above, for the scores being tested.
- **CLR normalisation** for non-count feature assays — Seurat, Stuart T et al.
  doi:[10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)
- **Nebulosa**, for the kernel-density feature plots — Alquicira-Hernandez J,
  Powell JE. *Bioinformatics* 2021;37:2485–2487.
  doi:[10.1093/bioinformatics/btab003](https://doi.org/10.1093/bioinformatics/btab003)
- Putting gene-set scores in an assay **so that the same test serves genes and
  gene sets** has `no precedent cited; chosen so the two analyses are
  comparable rather than being two different tests`.

## 06 — `06_fgsea.R`

Rankings from `FindMarkers(logfc.threshold = 0, min.pct = 0, test.use = "wilcox")`
— unfiltered on purpose — then `fgsea` per cell type.

The pathway was **swept by editing one line**: `GOI_SWEPT` in the parameter block
records the four values run (`GO_AEROBIC_RESPIRATION`,
`REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS`,
`SIRT3_TARGET_GENES`, `WP_NAD_METABOLISM_SIRTUINS_AND_AGING`). Surviving output
is two `.eps` NES barplots per pathway, eight files; the tables the script also
writes are not present on the cluster. 244 lines of a `VldlrKO`-vs-`CTL`
analysis belonging to another experiment were removed —
`METHODS_CHOICES.md` items 10 and 11.

**Citations.**
- **fgsea** — Korotkevich G et al. Fast gene set enrichment analysis. bioRxiv
  doi:[10.1101/060012](https://doi.org/10.1101/060012). Checked 2026-09-29:
  still a preprint, no published version. Cited as a preprint knowingly.
- **GSEA**, the method fgsea implements — Subramanian A et al. *PNAS*
  2005;102:15545–15550.
  doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
- Passing an **unfiltered ranking** follows the GSEA formulation above, in which
  the statistic is computed over the complete ordered list; it is not a
  departure.
