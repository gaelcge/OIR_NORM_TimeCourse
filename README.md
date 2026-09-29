# OIR_NORM_TimeCourse

Single-cell RNA-seq (Drop-seq) of the mouse retina across a developmental time
course in oxygen-induced retinopathy (OIR) and in normoxia, analysed to identify
the cellular mechanisms of proliferative retinopathy. The repository holds the
analysis code: integration and clustering of the merged count matrix, cell-type
annotation, differential expression by cell type and by condition, gene-set
activity, subclustering of individual compartments, trajectory inference, and
cell–cell interaction analysis.

The work spans roughly a decade of cluster analysis, over several major versions
of R and Seurat. That history is visible in the repository and is documented
rather than smoothed over — see `SOFTWARE_VERSIONS.md` and `EDIT_POLICY.md`.

---

## Data statement

**No data are committed to this repository.** Every script takes its input
locations from a parameter block at the top of the file.

| what | where | size |
|---|---|---|
| Count matrices | NCBI GEO, accession **GSE150703** | — |
| Merged Drop-seq matrix (all samples, all timepoints) | cluster, `ctb-jsjoyal` allocation, `Sequencing/Merging/TimeCourseOIR/Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt` | 4.67 GB |
| Seurat objects | cluster, beside the scripts that wrote them | 1.4–12.2 GB each |
| Reference gene sets and databases | cluster; specified by checksum in `reference/PROVENANCE.txt` | — |

To run anything here you need the merged matrix or the GEO deposit, and a copy
of the reference gene sets. `reference/PROVENANCE.txt` records the source,
release and sha256 of every reference input, and states plainly which ones are
**not obtainable** — the senescence gene-set analysis depends on three curated
files held by a collaborator and is not reproducible from this repository.

---

## Experimental design, and the constraint it places on every result

The main object carries **51,151 cells × 21,705 features** over six timepoints
and two conditions. The design is **not balanced**, and the imbalance bounds
what can be compared:

| cells | P5 | P7 | P10 | P12 | P14 | P17 |
|---|---|---|---|---|---|---|
| **Normoxia** | 3,022 | 1,557 | 1,863 | 2,097 | 17,959 | 3,487 |
| **OIR** | — | — | — | 3,179 | 9,895 | 8,092 |

**OIR cells exist only at P12, P14 and P17.** P5, P7 and P10 are normoxia only.
An OIR-versus-normoxia contrast is therefore available at three timepoints and
nowhere else; the earlier timepoints describe normal retinal development, not a
paired comparison. This is why the analyses subset to P12–P17 or to P14/P17
wherever a condition contrast is made.

The two sorted fractions are also incomplete:

| cells | P5 | P7 | P10 | P12 | P14 | P17 |
|---|---|---|---|---|---|---|
| `WR` | 709 | 1,557 | 1,863 | 2,295 | 21,949 | 6,979 |
| `Cd73ft` | 2,313 | — | — | 2,981 | 5,905 | 4,600 |

**P7 and P10 have no `Cd73ft` fraction.** Since integration anchors are computed
across `Sorting`, those two timepoints enter the integration through a single
fraction.

Other design facts carried forward in code: cells come from two laboratories
(37,766 Joyal, 13,385 McCarroll), across 22 sequencing batches and replicates
labelled `r1, r2, r3, r5, r8` — `r4`, `r6` and `r7` are absent from the object.
Batch is regressed out during normalisation.

### What the two fractions are

**`WR`** is whole retina — the full dissociated retina, unsorted.

**`Cd73ft`** is the CD73 flow-through: the CD73-negative fraction collected
after sorting on CD73 (Nt5e). CD73 marks photoreceptor precursors, so the
flow-through is depleted of photoreceptors and correspondingly enriched for the
rarer populations — vascular endothelium, pericytes, glia and immune cells.

This is why the two fractions must be integrated rather than pooled, and why
integration is anchored on `Sorting`: the fractions differ in composition **by
design**, so their difference is a sampling artefact to reconcile, not biology
to preserve. It also explains why the rare cell types in this dataset are
recoverable at all — 453 pericytes, 350 endothelial cells, 507 immune cells and
223 astrocytes out of 51,151 would be far scarcer without the enrichment.

Read the cell counts with this in mind: a cell-type proportion computed over all
cells is a proportion over a deliberately non-representative mixture, not over
the retina. Proportions in the composition figures are therefore comparable
*between conditions at matched timepoints*, where the fraction mix is similar,
and are not estimates of true retinal abundance.

### The OIR protocol

Pups and their nursing dam are exposed to **75% O₂ from P7 to P12**, then
returned to room air; neovascularisation peaks around **P17**. Normoxic
littermates remain in room air throughout.

This is what makes the design table above the shape it is. There is no OIR
before P12 because there is no retinopathy before P12 — the hyperoxic phase is
vessel loss, and the proliferative phase only begins after the return to room
air. The normoxia-only timepoints P5, P7 and P10 therefore describe normal
retinal development, and they are in the dataset as the developmental baseline
against which the P12–P17 OIR response is read, not as the missing half of a
paired comparison.

- Smith LE et al. *Invest Ophthalmol Vis Sci* 1994;35:101–111. PMID 7507904
- Connor KM et al. *Nat Protoc* 2009;4:1565–1573.
  doi:[10.1038/nprot.2009.187](https://doi.org/10.1038/nprot.2009.187)

---

## Layout

Folders are numbered by the order in which data enter the analysis, not by
figure number — figure numbers change during revision and data layers do not.

```
0_matrix/                 the deposited count matrix and its first QC view
1_timecourse_p5_p17/      P5-P17 integration, clustering, annotation, DE, gene sets
2_compartments_p14_p17/   P14/P17 integration, and all subclustering,
  compartments/           trajectory and cell-cell interaction work
reference/                PROVENANCE.txt only; the databases are not committed
figures/                  assembled multi-panel figures
```

### Two analyses, two objects — read this before combining anything

`1_timecourse_p5_p17` and `2_compartments_p14_p17` are **independent analyses of
the same merged matrix**, built years apart, and their outputs are not
interchangeable:

| | `1_timecourse_p5_p17` | `2_compartments_p14_p17` |
|---|---|---|
| object | `…AlignedBySorting.1.20.30` | `…AlignedByCond_Sorting.3.17.30` |
| timepoints | P5 → P17 | P14, P17 |
| integration anchored on | `Sorting` | `Cond_Sorting` |
| clustering | resolution 1 on the 2D UMAP, 59 clusters → 21 cell types | resolution 3 on 17 CCA dimensions, labelled by hand |
| Seurat | 3.1.5 | v2, under R 3.5.0 |
| establishes | the time course, DE and gene-set activity | the compartments: immune, neuroglial, vascular, pericyte, senescence |

**Cluster numbers and cell-type label sets do not transfer between them.** A
figure that puts a cluster identity from one beside a cluster identity from the
other is wrong. Each folder's `METHODS_CHOICES.md` opens by saying which object
it runs on.

| folder | what it establishes | n |
|---|---|---|
| `0_matrix` | the GSE150703 count matrix; not read by either analysis | P14/P17 subset |
| `1_timecourse_p5_p17` | cell-type composition across development and OIR; DE and gene-set activity per cell type | 51,151 cells |
| `2_compartments_p14_p17` | subtypes within each compartment; trajectories; ligand–receptor interactions | P14/P17 subset |

---

## Upstream: how the count matrix was made

This repository is the analysis stage. Everything below happened before it and
is described here rather than run from here.

**Library preparation and sequencing.** Drop-seq on dissociated mouse retina,
following the published protocol; cells were taken as whole retina (`WR`) and as
the CD73 flow-through (`Cd73ft`, the CD73-negative fraction, depleted of
photoreceptor precursors) at each timepoint — see the design section above. Sample identity is
encoded in the cell barcode column names of the merged matrix as
`Condition_TimePoint_Sorting_Labo_Replicate_CellBarCode`, which is what every
downstream script parses to build its metadata.

- Macosko EZ et al. Highly parallel genome-wide expression profiling of
  individual cells using nanoliter droplets. *Cell* 2015;161:1202–1214.
  doi:[10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002)
  — the Drop-seq method.

**Alignment and digital expression.** Reads were aligned and per-sample digital
gene expression (DGE) matrices generated with the Drop-seq tools pipeline. The
tool versions used were not recorded at the time and cannot be recovered from
the files; the archived sequencing data
(`TimeCourseOIR_sequencingfiles.tar.gz`, 39.8 GB) sits beside the merged matrix
on the cluster.

**Merging.** The per-sample DGEs were combined into the single 4.67 GB matrix by
`MergeDGE.R`, which lives beside its output and its Slurm wrapper in
`Sequencing/Merging/TimeCourseOIR/` on the `ctb-jsjoyal` allocation. It is not
committed here: it produced the input to this repository rather than being part
of it. Its run log (`MergeDGE.Rout`) is archived in the same directory.

**Gene symbols.** The merged matrix carries **uppercased mouse** gene symbols
(`MT-ATP6`, `0610005C13RIK`). This matters downstream: the MSigDB collections
used for gene-set analysis are human symbol sets, and they are joined to these
data by string identity with no orthology mapping. See item 7 of
`1_timecourse_p5_p17/METHODS_CHOICES.md`.

**The OIR model.** Oxygen-induced retinopathy is the standard mouse model of
proliferative retinopathy.

- Smith LE et al. Oxygen-induced retinopathy in the mouse. *Invest Ophthalmol
  Vis Sci* 1994;35:101–111. PMID 7507904 — the original model.
- Connor KM et al. Quantification of oxygen-induced retinopathy in the mouse: a
  model of vessel loss, vessel regrowth and pathological angiogenesis. *Nat
  Protoc* 2009;4:1565–1573.
  doi:[10.1038/nprot.2009.187](https://doi.org/10.1038/nprot.2009.187)

---

## Known caveats carried forward in code

These are the things a reader needs to know before using any output. Each is
documented in full in the `METHODS_CHOICES.md` of the folder it affects.

1. **Clusters in `1_timecourse_p5_p17` are cut from the 2D UMAP embedding**, not
   from principal components. The UMAP's parameters and seed are therefore
   clustering parameters. Report it as such.
2. **The senescence GSVA is not reproducible.** Three curated gene-set files are
   held by a collaborator and are unreadable; the step that scores them and the
   step that consumes those scores both depend on them.
3. **Differential-expression tables are filtered on the unadjusted p-value**
   while the accompanying volcano plots colour on the adjusted one. Filenames
   reading `p0.05` do not mean FDR 0.05. The adjusted column is present in every
   table.
4. **`0_matrix/01_make_dge_matrix.R` draws 6,000 cells without a seed**, so the
   exact cell set of the GEO deposit is not recoverable from the code — only
   from the deposit. Neither analysis folder depends on that subsample.
5. **Two MSigDB releases are in use** — v7.1 in `1_timecourse_p5_p17`, v6.2 in
   `2_compartments_p14_p17`. Enrichment results are not comparable term by term
   across the two folders.
6. **373 lines across 10 files of `2_compartments_p14_p17` are commented out as
   relics** of other projects, with a `#~ ` prefix and a banner. Nothing was
   deleted; if an expected result seems missing, search the `#~` blocks.
7. **CellPhoneDB is unversioned.** Neither the tool version nor its reference
   database release was recorded, so the interaction results cannot be
   reproduced to the version.
8. **Cell cycle is deliberately regressed out** of the P5–P17 integration, so
   that object cannot be used to study proliferation across development.
9. **The cell mixture is not representative of the retina by design.** Half the
   cells come from a CD73-depleted fraction, so a proportion computed over all
   cells is a proportion over an enriched mixture. Composition figures are
   comparable between conditions at matched timepoints; they are not estimates
   of true retinal abundance. The fraction mix also differs by timepoint — P7
   and P10 are `WR` only — so proportions are not comparable along the time
   axis either.

---

## Reproducing anything here

1. Obtain the count matrix (GSE150703) or the merged matrix, and the reference
   gene sets listed in `reference/PROVENANCE.txt`.
2. Read `SOFTWARE_VERSIONS.md`. These scripts were written against Seurat v2 and
   Seurat 3.1.5; several use functions (`SubsetData`, `MakeSparse`,
   `features.plot`, `min.genes`) that no longer exist in Seurat 5. Objects saved
   by them still load under Seurat 5.5.0.
3. Set the two path variables at the top of the step you want to run.
4. Run steps in numbered order within a folder. `EDIT_POLICY.md` states exactly
   what was changed relative to the code that ran on the cluster, and
   `PUBLISHED_TREE.md` maps every file here to its original path.

## Citing this repository

Cite the accession together with the paper that generated it; the accession
alone does not credit the work and the paper alone does not identify the release.
Code: cite this repository by URL and the commit you used — a bare repository
URL is not a citation, because the code cloned next year may not be the code
that produced the figure.
