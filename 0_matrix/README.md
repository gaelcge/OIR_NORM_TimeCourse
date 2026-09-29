# 0_matrix — the deposited count matrix

**What this layer establishes:** the raw UMI count matrix for the P14 and P17
samples that was deposited at GEO as **GSE150703**, and the first QC view of the
unfiltered data.

It is numbered first because it is where a reader starting from the public
deposit enters the project — **not** because the analyses depend on it. Neither
`1_timecourse_p5_p17` nor `2_compartments_p14_p17` reads its output; both read
the merged Drop-seq matrix directly.

| Step | Script | Consumes | Produces | Establishes |
|---|---|---|---|---|
| 01 | `01_make_dge_matrix.R` | merged matrix (4.67 GB, `ctb-jsjoyal`) | `retina_NORM_OIR_P14_P17_C57_WR_CD73FT_RawUMI_Count_DGEmatrix.txt`, `VlnPlotQC.Cond_Sorting.Init.png` | the GEO deposit, and that per-batch gene counts are comparable enough to integrate |

## 01 — `01_make_dge_matrix.R`

Selects the P14/P17 normoxic and OIR columns out of the merged matrix, writes
them as a tab-delimited count matrix, builds a Seurat object with
`min.cells = 3` / `min.genes = 100`, computes mitochondrial and crystallin
fractions, parses the sample metadata out of the column names, and renders one
violin plot of gene counts per batch.

Runtime: minutes. Memory is the constraint, not CPU — the input is a 4.67 GB
text matrix read with `data.table::fread` into a data frame.

**This step is not reproducible, and the reason is worth stating precisely.**
The McCarroll-lab columns are reduced with `sample(names(...), 6000)` and no
seed is set, so the exact 6,000 cells in the deposited matrix cannot be
recovered from this code. They are recoverable from the deposit itself, which is
why **GSE150703, not this script, is the authoritative definition of the P14/P17
cell set.** The two analysis folders are unaffected: they read all McCarroll
cells from the merged matrix.

**Seurat v2 API.** `raw.data`, `min.genes=`, `MakeSparse()`, `features.plot=`,
`x.lab.rot=`. This script will not run under Seurat 3 or later without being
rewritten. See `../SOFTWARE_VERSIONS.md`.

### Citations for this step

- **Drop-seq**, the assay that produced the counts — Macosko EZ et al. *Cell*
  2015;161:1202–1214.
  doi:[10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002)
- **Seurat**, for the object construction and the QC metric conventions
  (`percent.mito` as a fraction of UMIs on mitochondrial genes) — Butler A et
  al. *Nat Biotechnol* 2018;36:411–420.
  doi:[10.1038/nbt.4096](https://doi.org/10.1038/nbt.4096)
- **The data, by accession, together with the paper that generated it** — GEO
  **GSE150703**. The accession identifies the release; it does not credit the
  work, and the paper does not identify the release, so both belong in any
  citation of these data.

There is **no precedent cited for the crystallin QC filter** applied in the
analysis folders; it was chosen for this dataset because dissected retina
carries lens fragments. See item 3 of
`../1_timecourse_p5_p17/METHODS_CHOICES.md`.
