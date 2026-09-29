# Software versions

One project, ten years, several major versions of R and Seurat. There is no
single environment that runs this repository, and a `renv.lock` or
`environment.yml` would misrepresent it: the two analysis folders ran under
different Seurat major versions, and the scripts were never re-run under a
common one. What follows is what each part actually ran under, and how that is
known.

**Confidence is marked on every line.** `measured` means read off an artifact
that exists — an object's recorded version, a module load line in a committed
wrapper. `from the API` means inferred from functions the script calls that
exist only in a particular Seurat major version. `not recorded` means it was
not captured at the time and cannot now be recovered; those entries are the
gaps, and they are stated rather than guessed.

---

## R and Seurat, by folder

| folder | R | Seurat | how known |
|---|---|---|---|
| `0_matrix` | 3.5.x | v2 | **from the API** — `raw.data`, `min.genes=`, `MakeSparse()`, `features.plot=`, `x.lab.rot=`; none exist after Seurat 2 |
| `2_compartments_p14_p17` | **3.5.0** | v2 | **measured** — `module load r/3.5.0` in the committed Slurm wrappers; API matches (`SubsetData`, `CalcVarExpRatio`, `AlignSubspace`) |
| `2_compartments_p14_p17` GSVA steps | **3.6** | — | **measured** — `module load nixpkgs/16.09 gcc/7.3.0; module load r/3.6` in the GSVA wrappers |
| `1_timecourse_p5_p17` | not recorded | **3.1.5** | **measured** — the integrated object records `@version 3.1.5`; API matches (`min.features=`, `FindIntegrationAnchors`, `PrepSCTIntegration`) |

The Slurm wrappers are committed with their original `module load` lines
precisely because those lines are the only surviving record of the R version for
two of these steps.

### Reading the objects today

**measured 2026-09-29:** the `…AlignedBySorting.1.20.30` object (Seurat object
version 3.1.5, 12.2 GB) loads without error under **Seurat 5.5.0 /
SeuratObject 5.4.0** on R 4.4.0, and its `@meta.data`, `Assays()` and
`Reductions()` are all readable. So the objects survive; it is the *scripts*
that do not run under Seurat 5, because the functions they call were removed.

Functions used here that no longer exist in Seurat 5: `SubsetData`,
`MakeSparse`, `CalcVarExpRatio`, `AlignSubspace`, `SetIdent(value=)` in its v2
form, and the `features.plot=` / `x.lab.rot=` / `do.return=` / `pc.genes=`
argument names.

---

## Packages the analysis depends on

Versions were not pinned in any script, and `sessionInfo()` was not captured in
any surviving log, so **no package version below is recorded**. Each entry gives
the citation for the method, which is what a reader needs in order to know what
was done; where an independent benchmark exists that justifies the choice over
alternatives, it is given on its own line.

### Integration and normalisation
- **Seurat** — Butler A et al. Integrating single-cell transcriptomic data
  across different conditions, technologies, and species. *Nat Biotechnol*
  2018;36:411–420. doi:[10.1038/nbt.4096](https://doi.org/10.1038/nbt.4096)
  — the CCA alignment used by `2_compartments_p14_p17`.
- **Seurat v3 anchors** — Stuart T et al. Comprehensive integration of
  single-cell data. *Cell* 2019;177:1888–1902.
  doi:[10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)
  — `FindIntegrationAnchors` / `IntegrateData` in `1_timecourse_p5_p17`.
- **SCTransform** — Hafemeister C, Satija R. Normalization and variance
  stabilization of single-cell RNA-seq data using regularized negative binomial
  regression. *Genome Biol* 2019;20:296.
  doi:[10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1)
- **harmony** — Korsunsky I et al. Fast, sensitive and accurate integration of
  single-cell data with Harmony. *Nat Methods* 2019;16:1289–1296.
  doi:[10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)
  — loaded by most scripts; the committed integration paths use Seurat anchors
  rather than Harmony, so it is a dependency of the environment rather than of
  the published result.
- **DoubletFinder** — McGinnis CS, Murrow LM, Gartner ZJ. DoubletFinder:
  doublet detection in single-cell RNA sequencing data using artificial nearest
  neighbors. *Cell Syst* 2019;8:329–337.
  doi:[10.1016/j.cels.2019.03.003](https://doi.org/10.1016/j.cels.2019.03.003)
  — loaded, and commented out in the committed P5–P17 integration. Note the
  library is spelled `doubletFinder` in one script and `DoubletFinder` in
  another; only the latter is the package name.

### Cell cycle
- Tirosh I et al. Dissecting the multicellular ecosystem of metastatic melanoma
  by single-cell RNA-seq. *Science* 2016;352:189–196.
  doi:[10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)
  — the S and G2M gene lists, applied via `CellCycleScoring`. The cached copy is
  checksummed in `reference/PROVENANCE.txt`.

### Dimensionality reduction
- **UMAP** — McInnes L, Healy J, Melville J. UMAP: Uniform Manifold
  Approximation and Projection for dimension reduction. arXiv:1802.03426 (2018).
  *Preprint; no peer-reviewed version.*
- **t-SNE** — van der Maaten L, Hinton G. Visualizing data using t-SNE. *J Mach
  Learn Res* 2008;9:2579–2605.
- **PHATE** (`phateR`) — Moon KR et al. Visualizing structure and transitions in
  high-dimensional biological data. *Nat Biotechnol* 2019;37:1482–1492.
  doi:[10.1038/s41587-019-0336-3](https://doi.org/10.1038/s41587-019-0336-3)
  — loaded widely; no committed step produces a PHATE embedding that survives.
- **Louvain** — Blondel VD et al. Fast unfolding of communities in large
  networks. *J Stat Mech* 2008;P10008.
  doi:[10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)
  — `FindClusters(algorithm = 1)`.

### Imputation (tried, not used)
- **MAGIC** (`Rmagic`) — van Dijk D et al. Recovering gene interactions from
  single-cell data using data diffusion. *Cell* 2018;174:716–729.
  doi:[10.1016/j.cell.2018.05.061](https://doi.org/10.1016/j.cell.2018.05.061)
  — present and commented out; see item 5 of
  `2_compartments_p14_p17/METHODS_CHOICES.md` for why it was rejected.

### Differential expression
- **Wilcoxon rank-sum** as implemented in `FindMarkers(test.use = "wilcox")`.
- Soneson C, Robinson MD. Bias, robustness and scalability in single-cell
  differential expression analysis. *Nat Methods* 2018;15:255–261.
  doi:[10.1038/nmeth.4612](https://doi.org/10.1038/nmeth.4612)
  — the independent benchmark supporting a rank-based test for this purpose,
  rather than the tool's own paper.
- **limma** — Ritchie ME et al. limma powers differential expression analyses
  for RNA-sequencing and microarray studies. *Nucleic Acids Res* 2015;43:e47.
  doi:[10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)
- **Benjamini–Hochberg** — Benjamini Y, Hochberg Y. Controlling the false
  discovery rate. *J R Stat Soc B* 1995;57:289–300.
  doi:[10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
  — note that Seurat's `p_val_adj` is **Bonferroni** across assay genes, not BH;
  the distinction matters for how the tables are read.

### Gene sets
- **GSVA** — Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation
  analysis for microarray and RNA-Seq data. *BMC Bioinformatics* 2013;14:7.
  doi:[10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)
- **fgsea** — Korotkevich G et al. Fast gene set enrichment analysis. bioRxiv
  doi:[10.1101/060012](https://doi.org/10.1101/060012)
  — **checked 2026-09-29: still a preprint, no peer-reviewed version exists.**
  Cited as a preprint deliberately, not for lack of looking.
- **MSigDB** — Liberzon A et al. The Molecular Signatures Database hallmark gene
  set collection. *Cell Syst* 2015;1:417–425.
  doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004);
  Subramanian A et al. *PNAS* 2005;102:15545–15550.
  doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
  — releases **v7.1** and **v6.2** are both used; see
  `reference/PROVENANCE.txt`.
- **Nebulosa** — Alquicira-Hernandez J, Powell JE. Nebulosa recovers single-cell
  gene expression signals by kernel density estimation. *Bioinformatics*
  2021;37:2485–2487.
  doi:[10.1093/bioinformatics/btab003](https://doi.org/10.1093/bioinformatics/btab003)
- **SingleR** — Aran D et al. Reference-based analysis of lung single-cell
  sequencing reveals a transitional profibrotic macrophage. *Nat Immunol*
  2019;20:163–172.
  doi:[10.1038/s41590-018-0276-y](https://doi.org/10.1038/s41590-018-0276-y)
  — loaded in the annotation step; the committed labels are assigned by hand
  from markers, so SingleR informs rather than produces them.

### Trajectory
- **Monocle 2 / DDRTree** — Qiu X et al. Reversed graph embedding resolves
  complex single-cell trajectories. *Nat Methods* 2017;14:979–982.
  doi:[10.1038/nmeth.4402](https://doi.org/10.1038/nmeth.4402)
- Trapnell C et al. The dynamics and regulators of cell fate decisions.
  *Nat Biotechnol* 2014;32:381–386.
  doi:[10.1038/nbt.2859](https://doi.org/10.1038/nbt.2859)
  — the original pseudotime formulation.
- **dyno / dynverse** — Saelens W et al. A comparison of single-cell trajectory
  inference methods. *Nat Biotechnol* 2019;37:547–554.
  doi:[10.1038/s41587-019-0071-9](https://doi.org/10.1038/s41587-019-0071-9)
  — both the tool and the independent benchmark that justifies comparing
  trajectory methods rather than trusting one.

### Cell–cell interaction
- **NicheNet** — Browaeys R, Saelens W, Saeys Y. NicheNet: modeling
  intercellular communication by linking ligands to target genes. *Nat Methods*
  2020;17:159–162.
  doi:[10.1038/s41592-019-0667-5](https://doi.org/10.1038/s41592-019-0667-5)
  — **the reference networks used cannot be identified.** They were read from a
  collaborator's copy on a path that is no longer readable, and no release was
  recorded. The method is specified; the exact networks are not.
- **CellPhoneDB** — Efremova M et al. CellPhoneDB: inferring cell–cell
  communication from combined expression of multi-subunit ligand–receptor
  complexes. *Nat Protoc* 2020;15:1484–1506.
  doi:[10.1038/s41596-020-0292-x](https://doi.org/10.1038/s41596-020-0292-x)
  — **unversioned.** Run on the command line inside a python 3.7 virtualenv
  between two committed scripts; neither the CellPhoneDB version nor its
  database release was recorded, and the virtualenv no longer exists. The
  interaction results cannot be reproduced to the version.
- **ktplots** — used for the CellPhoneDB heatmaps. Repository, no archived
  release DOI; the commit used was not recorded.

### Python
- A **python 3.7 virtualenv** (`ENV_python3.7.0`) was used for CellPhoneDB and
  for the python GSVA steps. Its package set was not captured and the
  environment is gone. The python GSVA scripts import `pysnow`, `flask` and
  `flask_snow`, which are unrelated template leftovers but are required for the
  files to import as written.

---

## What is missing, plainly

- No `sessionInfo()` from any run.
- No package versions for anything above.
- No Drop-seq tools version for the upstream alignment.
- No CellPhoneDB version or database release.
- No NicheNet network release.
- No record of which MSigDB release goes with which figure, beyond the
  folder-level split (v7.1 / v6.2).

Nothing here can be recovered by inspecting the files, so none of it is
guessed. For future work in this project, capturing `sessionInfo()` into the
step's output directory costs one line and closes most of this list.
