# Methods choices — compartments, P14 and P17

This folder is the earlier of the two analyses and the messier one. It is also
where the subclustering, trajectory and cell–cell interaction work lives, which
exists nowhere else in the project.

Read item 1 before anything else in this folder.

---

## 1. This folder runs on a different object from `1_timecourse_p5_p17`

Everything here reads
`RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30`:

| | this folder | `1_timecourse_p5_p17` |
|---|---|---|
| timepoints | P14, P17 | P5, P7, P10, P12, P14, P17 |
| integration grouped by | `Cond_Sorting` | `Sorting` |
| resolution / dims | 3 / 17 | 1 / 20 |
| Seurat | v2, under R 3.5.0 | 3.1.5 |
| cells | — | 51,151 |

Both are integrations of the same merged Drop-seq matrix, so the two are not
independent data — but they are independent analyses. **Cluster numbers do not
transfer between them, and neither do cell-type label sets.** A figure combining
a cluster identity from one with a cluster identity from the other would be
wrong.

---

## 2. Integration is anchored on condition crossed with sorted fraction

`01_integrate_by_cond_sorting.R` splits on `Cond_Sorting` — the four
condition × fraction groups — and anchors CCA across them.

**Why this is worth flagging.** Anchoring on a grouping that includes condition
aligns the OIR and normoxic subspaces onto each other. That makes cell types
directly comparable across condition, which is what the subclustering work
needs. It also shrinks the condition effect, which is the quantity the
differential-expression steps then measure — the effect is partly removed by
the same step that makes the comparison possible.

The later analysis anchors on `Sorting` alone for exactly this reason. Both
choices are defensible for their purpose; they should not be described as the
same method.

---

## 3. Resolution 3 over-partitions on purpose

`02_cluster.R` clusters at resolution 3 on 17 CCA dimensions, and
`03_annotate_celltypes.R` then collapses the result by hand from marker
expression.

**Why.** Over-clustering and merging puts the biological judgement in the
labelling step, where it is visible in a named mapping, rather than in a
resolution parameter where it is invisible. The cost is that a cluster count at
this resolution carries no meaning on its own and should not be reported as a
number of cell populations.

---

## 4. Cells poorly explained by CCA were inspected, then kept

Two integrations exist: `01_integrate_by_cond_sorting.R` (the pipeline) and
`01b_integrate_with_cca_discard.R` (a variant that removes cells whose
`var.ratio.pca` is below 0.5, i.e. whose profile is better explained by PCA than
by the shared CCA subspace).

`02b_cluster_discarded.R` embeds the discarded cells and plots PECAM1, CLDN5
and MKI67 on them, and `04_cluster_markers.R` computes their markers. The
question was whether the filter preferentially removes a real population —
endothelium and proliferating cells being the ones at risk — rather than
low-quality droplets.

**The pipeline reads the unfiltered object.** The discarding variant is in the
repository because the evidence for not discarding is what it produced. Do not
mistake `01b`/`02b` for the main line.

---

## 5. MAGIC imputation was tried and rejected

`07_gsva_prepare.R` contains a MAGIC imputation block, commented out in the
original. Every downstream path is named `NotImputed` for this reason.

**Why rejected.** Imputation raises the within-cell correlation between genes
that share a pathway — which is precisely the quantity GSVA scores. Gene-set
enrichment on imputed values is therefore inflated by construction. The
un-imputed matrix is the conservative input.

The commented block is left in place as the record that the alternative was
considered.

---

## 6. 1000 cells per identity

`07_gsva_prepare.R` also writes a subsample of 1000 cells per cell-type
identity, which the `Subset1000CellperIdent` outputs are built from.

**Why.** Cell-type sizes in this dataset differ by roughly two orders of
magnitude. Without equalising, a gene-set score distribution over all cells
describes rods. Equalising costs power in the abundant types and is the reason
two GSVA variants exist rather than one.

---

## 7. MSigDB v6.2 here, v7.1 in the sibling folder

This folder's gene-set steps read
`Gene_list/Gmt.file/{h.all,c3.all,sup_h.c2.c5.c6.c7.all}.v6.2.symbols.gmt`
and one file dated `270519`. The P5–P17 analysis reads v7.1.

**Do not substitute one for the other.** Between releases MSigDB adds, removes
and renames sets, so a term present in one analysis may be absent or redefined
in the other, and enrichment results are not comparable across the two folders
term by term. Both releases and their checksums are in
`reference/PROVENANCE.txt`.

---

## 8. Two GSVA implementations, both of which ran

`08`/`08b` use the Bioconductor GSVA package on a senescence-focused collection
at P17. `09`/`09b` use the python GSVA package on the broad MSigDB collections
across all cells. Both produced output that is on the cluster — from 2.5 MB to
2.9 GB.

**A known oddity in the python scripts, left as it was.** They import `pysnow`,
`flask` and `flask_snow`, which have nothing to do with gene-set scoring. Those
are template leftovers and were present in the scripts as run; the scripts will
not import without those packages installed.

---

## 9. Compartments are subclustered independently

Each of the six compartments under `compartments/` is re-embedded and
re-clustered on its own cells, and `11_subcluster_convergence.R` reconciles the
label sets onto one object afterwards.

**Why not one clustering.** No single resolution resolves microglial activation
states and rod subtypes at the same time: the resolution that splits the former
shatters the latter. Independent subclustering is the standard way out, and the
convergence step is where its cost is paid — the sub-labels come from six
separate embeddings, so a distance between two cells in different compartments
has no meaning.

**A gap in the chain.** Several compartment scripts read objects from
`CellIdentification/` and `CellIdentification_5/` directories. The `_5` suffix
indicates a fifth iteration, and the earlier iterations sit in
`ImmuneCells/Deprecated/Mapping_1` through `Mapping_5` on the cluster, which
this repository does not publish. So for those compartments the committed
script reads an object whose producing script is not in the repository. Stated
rather than papered over.

---

## 10. 373 lines are commented out as relics from other projects

60 blocks in 10 files are disabled with a `#~ ` prefix. The two qualifying
criteria and the full reasoning are in `../EDIT_POLICY.md` item 6. In summary:
`setwd()` calls into directories that no longer exist, and blocks operating on
objects the script never loads — several of which filter on
`Dataset %in% c("OIR.P17.CD31.Sirt3KO.S129", ...)`, a column absent from this
object, identifying them as code from a CD31-sorted Sirt3 knockout experiment.

Worst affected: `compartments/02_neuroglial/02_de_by_condition.R`, where 151
lines operate on vascular-endothelium objects from that other experiment.

Nothing was deleted. If a result you expected is missing, the code for it may be
behind a `#~` — search for it before assuming it was never written.

---

## 11. The NicheNet step depends on files that are not reachable

`compartments/06_senescence/04_nichenet.Rmd` reads all three NicheNet reference
networks from a heart-maturation project belonging to another user on a
different allocation, and writes one figure back into that user's directory.
Those paths are **not readable** by this account (checked 2026-09-29) and are
commented out as relics.

The networks themselves are the standard published NicheNet references, so the
analysis is specified even though the paths are dead — see the layer README for
the citation and the official deposit. What cannot be recovered is which
release of those networks was used.

---

## 12. CellPhoneDB was run on two cell groupings

`05_cellphonedb_oir_p17.py` and `06_cellphonedb_oir_p17_kras.py` differ in the
cell grouping they export as CellPhoneDB input — the second adds a KRAS-defined
grouping. Both write into `CellphoneDB/input/` and read back from
`CellphoneDB/output/`, and `07_cellphonedb_heatmap.R` renders the result.

The CellPhoneDB run itself is not in this repository: it happens between the
two, on the command line, inside a python 3.7 virtualenv. The wrapper that
activates that environment is committed; the environment is not, and no version
of CellPhoneDB or of its reference database was recorded at the time. That
means the interaction results cannot be reproduced to the version — a real gap,
and the reason `SOFTWARE_VERSIONS.md` marks CellPhoneDB as unversioned.
