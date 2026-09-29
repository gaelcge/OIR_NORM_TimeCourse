# 2_compartments_p14_p17 — compartments at P14 and P17

**What this layer establishes:** the cellular detail of the OIR response —
subtypes within each retinal compartment, the trajectories that connect them,
and the ligand–receptor interactions between them. This is the only place in the
project where subclustering, pseudotime and cell–cell interaction analysis
exist.

**The object:** `RETINA_RYTVELA_OIRTimeCourse_AlignedByCond_Sorting.3.17.30.Seurat_object.rds`
→ annotated as `…Seurat_object.annotated.rds` by step 03. P14 and P17 only,
CCA-integrated across the four condition × fraction groups, clustered at
resolution 3 on 17 dimensions, Seurat v2 under R 3.5.0.

> **This is a different object from the one in `1_timecourse_p5_p17`,** built
> years earlier with a different integration, a different resolution and a
> different Seurat major version. Cluster numbers and label sets do not
> transfer. `METHODS_CHOICES.md` item 1 has the comparison table; read it first.

This folder is the older and rougher of the two. **373 lines across 10 files are
commented out with a `#~ ` prefix** — relics carried in from other projects,
disabled rather than deleted. If a result you expect seems absent, search the
`#~` blocks before concluding it was never written. `../EDIT_POLICY.md` item 6.

## Run order

Steps `01`–`10` at folder level, then everything under `compartments/`, then
**`11`** — which is numbered last because it consumes the per-compartment
subclustered objects and merges their labels back onto one object.

`01b` and `02b` are a **variant, not the main line.** They apply Seurat v2's
CCA-versus-PCA variance-ratio filter and set the removed cells aside; step 04
then characterises those cells. The pipeline object is the one from `02`. The
variant is in the repository because the evidence for *not* discarding is what
it produced — see `METHODS_CHOICES.md` item 4.

| Step | Script | Consumes | Produces | Establishes |
|---|---|---|---|---|
| 01 | `01_integrate_by_cond_sorting.R` | merged matrix; cell-cycle list | `Clustering_noDiscarding/…Seurat_object.rds` | the P14/P17 object, all cells retained |
| 01b | `01b_integrate_with_cca_discard.R` | merged matrix; cell-cycle list | `Clustering/…Seurat_object.rds` | the same integration with `var.ratio.pca < 0.5` cells removed |
| 02 | `02_cluster.R` | step 01 | clustered object, cluster figures | the resolution-3 partition that step 03 annotates |
| 02b | `02b_cluster_discarded.R` | step 01b | `Clustering/discarded/Seurat_object.integrated.discard_17CC.rds` | that the discarded cells are not a coherent population — PECAM1/CLDN5/MKI67 checked |
| 03 | `03_annotate_celltypes.R` | step 02 | `CellTypeAnnotation/…Seurat_object.annotated.rds` | the cell-type labels every later step reads |
| 04 | `04_cluster_markers.R` | steps 02 **and 02b** | marker tables and figures | the marker evidence behind the step 03 labels, and what the CCA filter would have cost |
| 05 | `05_de_by_celltype.R` | step 03 | DE tables and figures per cell type | the cell-type axis, pooled across condition |
| 06 | `06_de_by_condition.R` | step 03 | DE tables and figures per condition | the OIR response within each cell type |
| 07 | `07_gsva_prepare.R` | step 03 | expression matrices, and a 1000-cells-per-identity subsample | the GSVA inputs; that MAGIC imputation was tried and rejected |
| 08 | `08_gsva_p17_senescence.R` | step 07 | `retina_p17_data_gsva_out_senescence_May2020.csv` | per-cell senescence scores at P17 |
| 08b | `08b_gsva_1000cells_senescence.R` | step 07 | `…senescence_feb2020.csv` | the size-equalised version of 08 (an earlier attempt, not a sweep) |
| 09 | `09_gsva_msigdb.py` | step 07 | `gsva.exprs.Full_eset.NonImputed.sup_h.c2.c5.c6.c7.csv` | broad MSigDB scores via the python implementation |
| 09b | `09b_gsva_msigdb_1000cells.py` | step 07 | `…scaled.c3.all.v6.2.symbols.csv` | the same on the subsample |
| 10 | `10_fgsea.R` | step 03; step 06 contrasts | fgsea tables and NES figures | pathway enrichment per cell type |
| 11 | `11_subcluster_convergence.R` | step 03 **and every compartment** | the converged object and its figures | one object carrying every compartment's fine labels |

### compartments/

Each compartment is re-embedded and re-clustered **on its own cells**, because
no single resolution resolves microglial activation states and rod subtypes at
the same time. The cost is that sub-labels come from six separate embeddings, so
a distance between cells in different compartments is meaningless —
`METHODS_CHOICES.md` item 9.

| folder | scripts | establishes |
|---|---|---|
| `01_immune_cells` | `02_de_by_subtype.R`, `03_de_by_condition.R`, `04_gsva_prepare.R` | microglial and immune subtypes and their OIR response |
| `02_neuroglial` | `01_subcluster.R`, `02_de_by_condition.R`, `03_monocle_trajectory.R`, `04_dyno_trajectory.Rmd` | neuronal and glial subtypes, and a pseudotime ordering through them |
| `03_vascular_endothelium` | `01_subcluster.R` | endothelial subtypes including tip and proliferating states |
| `04_pericytes` | `01_subcluster.R` | pericyte heterogeneity |
| `05_glial_ecs` | `01_subcluster.R`, `02_monocle_trajectory.R` | glia and endothelium analysed jointly, and their trajectory |
| `06_senescence` | `01_subcluster.R`, `02_convergence.R`, `03_gsva.R`, `04_nichenet.Rmd`, `05_cellphonedb_oir_p17.py`, `06_cellphonedb_oir_p17_kras.py`, `07_cellphonedb_heatmap.R` | the senescence-associated populations and their signalling to other compartments |

**A gap in the chain, stated rather than hidden.** Several compartment scripts
read objects from `CellIdentification/` and `CellIdentification_5/`. The `_5`
suffix is a fifth iteration; the earlier iterations sit in
`ImmuneCells/Deprecated/Mapping_1` … `Mapping_5` on the cluster and are **not**
published here. For those compartments the committed script therefore reads an
object whose producing script is absent from this repository.
`METHODS_CHOICES.md` item 9.

**`01_immune_cells` has no `01_`.** Its subclustering script is one of the
deprecated-mapping iterations above, so the numbering starts at `02`. The number
is left unused rather than renumbering the steps that follow.

---

## Citations by step

### 01, 01b, 02, 02b — integration and clustering
- **CCA alignment**, the integration method of this lineage — Butler A et al.
  Integrating single-cell transcriptomic data across different conditions,
  technologies, and species. *Nat Biotechnol* 2018;36:411–420.
  doi:[10.1038/nbt.4096](https://doi.org/10.1038/nbt.4096).
  `CalcVarExpRatio` and the `var.ratio.pca < 0.5` threshold of step 01b are that
  paper's own alignment diagnostic and its suggested cut.
- **Louvain community detection** — Blondel VD et al. *J Stat Mech* 2008;P10008.
  doi:[10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)
- **The cell-cycle gene lists** — Tirosh I et al. *Science* 2016;352:189–196.
  doi:[10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)
- Anchoring on **`Cond_Sorting` rather than on fraction alone** is a departure
  from the later analysis in this same project and has `no external precedent
  cited; chosen because the compartment work needs cell types comparable across
  condition`. Its cost — partial removal of the condition effect the DE steps
  then measure — is in `METHODS_CHOICES.md` item 2.
- **Resolution 3 as deliberate over-partitioning** has `no precedent cited;
  chosen so the biological judgement sits in the named label mapping of step 03
  rather than in a resolution parameter`.

### 03, 04 — annotation and markers
- **Seurat** for `FindAllMarkers` — Butler A et al., as above.
- **SingleR**, loaded to inform the labelling — Aran D et al. *Nat Immunol*
  2019;20:163–172.
  doi:[10.1038/s41590-018-0276-y](https://doi.org/10.1038/s41590-018-0276-y)

### 05, 06 — differential expression
- **Wilcoxon rank-sum** via `FindMarkers` — Butler A et al., as above.
- **Why a rank-based test**, from an independent benchmark — Soneson C, Robinson
  MD. *Nat Methods* 2018;15:255–261.
  doi:[10.1038/nmeth.4612](https://doi.org/10.1038/nmeth.4612)
- **limma**, used in parts of these steps — Ritchie ME et al. *Nucleic Acids
  Res* 2015;43:e47.
  doi:[10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

### 07, 08, 08b, 09, 09b — gene-set activity
- **GSVA** — Hänzelmann S, Castelo R, Guinney J. *BMC Bioinformatics* 2013;14:7.
  doi:[10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)
- **MSigDB v6.2** — Liberzon A et al. *Cell Syst* 2015;1:417–425.
  doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004);
  Subramanian A et al. *PNAS* 2005;102:15545–15550.
  doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).
  **v6.2 here, v7.1 in the sibling folder** — see `../reference/PROVENANCE.txt`.
- **MAGIC**, present and commented out — van Dijk D et al. *Cell*
  2018;174:716–729.
  doi:[10.1016/j.cell.2018.05.061](https://doi.org/10.1016/j.cell.2018.05.061).
  Rejected here; the reason is a departure and is in `METHODS_CHOICES.md` item 5.
- **Subsampling to 1000 cells per identity** has `no precedent cited; chosen
  because cell-type sizes differ by two orders of magnitude and an unequalised
  score distribution describes rods`.

### 10 — fgsea
- **fgsea** — Korotkevich G et al. bioRxiv
  doi:[10.1101/060012](https://doi.org/10.1101/060012). Still a preprint as of
  2026-09-29.
- **GSEA**, the method — Subramanian A et al.
  doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

### compartments — trajectory
- **Monocle 2 / DDRTree**, used by `03_monocle_trajectory.R` and
  `05_glial_ecs/02_monocle_trajectory.R` — Qiu X et al. Reversed graph embedding
  resolves complex single-cell trajectories. *Nat Methods* 2017;14:979–982.
  doi:[10.1038/nmeth.4402](https://doi.org/10.1038/nmeth.4402)
- **The original pseudotime formulation** — Trapnell C et al. *Nat Biotechnol*
  2014;32:381–386. doi:[10.1038/nbt.2859](https://doi.org/10.1038/nbt.2859)
- **dyno**, used by `04_dyno_trajectory.Rmd`, and simultaneously the independent
  benchmark that justifies checking a trajectory against more than one
  method — Saelens W et al. A comparison of single-cell trajectory inference
  methods. *Nat Biotechnol* 2019;37:547–554.
  doi:[10.1038/s41587-019-0071-9](https://doi.org/10.1038/s41587-019-0071-9)
- **PHATE** (`phateR`) is loaded across these scripts — Moon KR et al. *Nat
  Biotechnol* 2019;37:1482–1492.
  doi:[10.1038/s41587-019-0336-3](https://doi.org/10.1038/s41587-019-0336-3)
  — no committed step produces a surviving PHATE embedding.

### compartments/06_senescence — cell–cell interaction
- **NicheNet** — Browaeys R, Saelens W, Saeys Y. *Nat Methods* 2020;17:159–162.
  doi:[10.1038/s41592-019-0667-5](https://doi.org/10.1038/s41592-019-0667-5).
  The published reference networks are deposited by the authors under their own
  DOI; **the copy used here cannot be identified** — it was read from a
  collaborator's unreadable path and no release was recorded.
  `METHODS_CHOICES.md` item 11.
- **CellPhoneDB** — Efremova M et al. *Nat Protoc* 2020;15:1484–1506.
  doi:[10.1038/s41596-020-0292-x](https://doi.org/10.1038/s41596-020-0292-x).
  **Unversioned**: the run happened on the command line between the two
  committed scripts, in a python 3.7 environment that no longer exists, and
  neither the tool version nor its database release was recorded. The
  interaction results cannot be reproduced to the version.
  `METHODS_CHOICES.md` item 12.
- **ktplots**, for the heatmaps — public repository, no archived release DOI,
  and the commit used was not recorded. A repository URL alone is not a
  citation, so this is recorded as a gap rather than presented as one.
