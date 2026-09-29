# Published tree — provenance of every file

Every path in this repository, mapped to the absolute path of the file it came
from on Narval, so any file here can be diffed against its original.

`dropped` counts lines deleted (all documented in the relevant
`METHODS_CHOICES.md`); `relic` counts lines commented out with `#~` as code
carried in from other projects (`EDIT_POLICY.md` item 6). Where both are 0, the
only changes are the header, the parameter block and the path substitutions.

| repository path | original path on Narval | dropped | relic |
|---|---|---|---|
| `0_matrix/01_make_dge_matrix.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/MakeDGEmatrix.R` |  |  |
| `1_timecourse_p5_p17/01_integrate_by_sorting.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering/MakeSeurat.R` |  |  |
| `1_timecourse_p5_p17/01_integrate_by_sorting.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering/ScriptR1.sh` |  |  |
| `1_timecourse_p5_p17/02_cluster_annotate.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering/MakeSeurat.2.R` | 4 |  |
| `1_timecourse_p5_p17/03_gsva_senescence.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR/run_gsva.R` |  |  |
| `1_timecourse_p5_p17/03_gsva_senescence.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR/ScriptR.sh` |  |  |
| `1_timecourse_p5_p17/03b_gsva_oir_only.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR/Cpu/run_gsva.R` |  |  |
| `1_timecourse_p5_p17/03b_gsva_oir_only.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/GSVA/P14_P17_OIR/Cpu/ScriptR.sh` |  |  |
| `1_timecourse_p5_p17/04_de_genes.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/DifferentialExpression/Genes/DifferentialAnalysis.R` |  |  |
| `1_timecourse_p5_p17/05_de_functions.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/DifferentialExpression/Functions/DifferentialAnalysis.R` |  |  |
| `1_timecourse_p5_p17/06_fgsea.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/DifferentialExpression/fgsea/fgsea_working.R` |  |  |
| `1_timecourse_p5_p17/annotation/celltype_grouping.csv` | `extracted from the RenameIdents blocks of DE/Genes, DE/Functions and fgsea (verified identical)` |  |  |
| `1_timecourse_p5_p17/annotation/cluster_annotation.csv` | `extracted from Clustering/MakeSeurat.2.R (cell_type_assigned vector)` |  |  |
| `1_timecourse_p5_p17/functions/seurat_to_loom.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering/SeuratToLoom.R` |  |  |
| `2_compartments_p14_p17/01_integrate_by_cond_sorting.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Clustering_noDiscarding/MakeSeurat.R` |  |  |
| `2_compartments_p14_p17/01_integrate_by_cond_sorting.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Clustering_noDiscarding/ScriptR1.sh` |  |  |
| `2_compartments_p14_p17/01b_integrate_with_cca_discard.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Clustering/MakeSeurat.R` |  |  |
| `2_compartments_p14_p17/01b_integrate_with_cca_discard.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Clustering/ScriptR1.sh` |  |  |
| `2_compartments_p14_p17/02_cluster.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Clustering_noDiscarding/MakeSeurat.2.R` |  |  |
| `2_compartments_p14_p17/02b_cluster_discarded.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Clustering/MakeSeurat.2.R` |  |  |
| `2_compartments_p14_p17/03_annotate_celltypes.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/CellTypeAnnotation/CellIDs.R` |  |  |
| `2_compartments_p14_p17/04_cluster_markers.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/MarkerIdentification/FindClustermarkers.R` |  |  |
| `2_compartments_p14_p17/04_cluster_markers.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/MarkerIdentification/ScriptR1.sh` |  |  |
| `2_compartments_p14_p17/05_de_by_celltype.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/DifferentialAnalysis/CellType/Differential_analysis.R` |  | 12 |
| `2_compartments_p14_p17/06_de_by_condition.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/DifferentialAnalysis/Condition/Differential_analysis.R` |  | 46 |
| `2_compartments_p14_p17/07_gsva_prepare.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/CellGSVA.part1.R` |  |  |
| `2_compartments_p14_p17/08_gsva_p17_senescence.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/run_gsva.R` |  |  |
| `2_compartments_p14_p17/08_gsva_p17_senescence.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/ScriptR.sh` |  |  |
| `2_compartments_p14_p17/08b_gsva_1000cells_senescence.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/Subset1000CellperIdent/run_gsva.R` |  |  |
| `2_compartments_p14_p17/08b_gsva_1000cells_senescence.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/Subset1000CellperIdent/ScriptR.sh` |  |  |
| `2_compartments_p14_p17/09_gsva_msigdb.py` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/run_GSVA.py` |  |  |
| `2_compartments_p14_p17/09_gsva_msigdb.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/ScriptPython.sh` |  |  |
| `2_compartments_p14_p17/09b_gsva_msigdb_1000cells.py` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/GSVA/NotImputed/Subset1000CellperIdent/run_GSVA.py` |  |  |
| `2_compartments_p14_p17/10_fgsea.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/fgsea/fgsea_working.R` |  |  |
| `2_compartments_p14_p17/11_subcluster_convergence.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/Retina_Subcluster_convergence.R` |  |  |
| `2_compartments_p14_p17/compartments/01_immune_cells/02_de_by_subtype.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/ImmuneCells/DifferentialAnalysis/CellType/DifferentialAnalysis.R` |  | 14 |
| `2_compartments_p14_p17/compartments/01_immune_cells/03_de_by_condition.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/ImmuneCells/DifferentialAnalysis/Conditions/DifferentialAnalysis.R` |  | 6 |
| `2_compartments_p14_p17/compartments/01_immune_cells/04_gsva_prepare.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/ImmuneCells/GSVA/CellGSVA.part1.R` |  |  |
| `2_compartments_p14_p17/compartments/01_immune_cells/04_gsva_prepare.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/ImmuneCells/GSVA/NotImputed/ScriptR1.sh` |  |  |
| `2_compartments_p14_p17/compartments/02_neuroglial/01_subcluster.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/NeuroGlialCell/Retina_analysis.R` |  |  |
| `2_compartments_p14_p17/compartments/02_neuroglial/02_de_by_condition.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/NeuroGlialCell/DifferentialAnalysis/Condition/DifferentialAnalysis.R` |  | 151 |
| `2_compartments_p14_p17/compartments/02_neuroglial/03_monocle_trajectory.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/NeuroGlialCell/Monocle/Monocle_analysis.R` |  | 1 |
| `2_compartments_p14_p17/compartments/02_neuroglial/04_dyno_trajectory.Rmd` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/NeuroGlialCell/Dyno/Dyno_analysis.test.Rmd` |  |  |
| `2_compartments_p14_p17/compartments/03_vascular_endothelium/01_subcluster.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/VascularEndothelium/Mapping/Retina_analysis.R` |  |  |
| `2_compartments_p14_p17/compartments/04_pericytes/01_subcluster.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/Pericytes/Retina_analysis.R` |  |  |
| `2_compartments_p14_p17/compartments/05_glial_ecs/01_subcluster.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/Glial_ECs/Retina_analysis.R` |  | 75 |
| `2_compartments_p14_p17/compartments/05_glial_ecs/02_monocle_trajectory.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/Glial_ECs/Monocle/Monocle_analysis.R` |  | 2 |
| `2_compartments_p14_p17/compartments/06_senescence/01_subcluster.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/Mapping/Retina_analysis.R` |  | 62 |
| `2_compartments_p14_p17/compartments/06_senescence/02_convergence.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/Retina_Subcluster_convergence.R` |  |  |
| `2_compartments_p14_p17/compartments/06_senescence/03_gsva.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/GSVA/run_gsva.R` |  |  |
| `2_compartments_p14_p17/compartments/06_senescence/03_gsva.slurm.sh` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/GSVA/ScriptR.sh` |  |  |
| `2_compartments_p14_p17/compartments/06_senescence/04_nichenet.Rmd` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/run_nichenet_all_celltypes.Rmd` |  | 4 |
| `2_compartments_p14_p17/compartments/06_senescence/05_cellphonedb_oir_p17.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/CellphoneDB/result_plot/OIR_P17/CellPhoneDB.py` |  |  |
| `2_compartments_p14_p17/compartments/06_senescence/06_cellphonedb_oir_p17_kras.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/CellphoneDB/result_plot/OIR_P17_KRAS/CellPhoneDB.py` |  |  |
| `2_compartments_p14_p17/compartments/06_senescence/07_cellphonedb_heatmap.R` | `/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/Cond_Sorting/Subclustering/MergingSenescence/CellphoneDB/result_plot/heatmap_cellphonedb.R` |  |  |

## Renamed on the way in

Two files were named `.py` on the cluster and are not python. Measured: 82 R
statements, 9 shell lines and 0 python statements in each. They are R scripts
with the CellPhoneDB command line embedded as a shell block. Renamed to `.R`;
the shell block is left in place, since it is the only record of how CellPhoneDB
was invoked.

| repository path | original name |
|---|---|
| `…/06_senescence/05_cellphonedb_oir_p17.R` | `CellPhoneDB.py` |
| `…/06_senescence/06_cellphonedb_oir_p17_kras.R` | `CellPhoneDB.py` |

## Files written for this repository (no cluster original)

| path | what it is |
|---|---|
| `README.md` | the study, the design table, the data statement, the caveats |
| `EDIT_POLICY.md` | the six permitted classes of change, and what was not changed |
| `SOFTWARE_VERSIONS.md` | R/Seurat versions per folder, with how each is known, and the gaps |
| `PUBLISHED_TREE.md` | this file |
| `.gitignore` | excludes every data extension plus `data/`, `raw/`, `reference/` |
| `reference/PROVENANCE.txt` | source, release and sha256 of every reference input; states which are unobtainable |
| `0_matrix/README.md` | step table and citations |
| `1_timecourse_p5_p17/README.md` | step table and per-step citations |
| `1_timecourse_p5_p17/METHODS_CHOICES.md` | 12 numbered choices a reader could mistake for errors |
| `2_compartments_p14_p17/README.md` | step table, compartment table and per-step citations |
| `2_compartments_p14_p17/METHODS_CHOICES.md` | 12 numbered choices, opening with the two-object warning |

## Not published, and why

| on Narval | why it is not here |
|---|---|
| `Aligned/Cond_Sorting/Clustering_noDiscarding/MakeSeurat.2.R` variants under `Deprecated/` | superseded lineages; 17 scripts under `Deprecated/`, `Aligned/Deprecated/` and `SeuratV3/Deprecated/` |
| `Aligned/SeuratV3/Deprecated/Aligned_Cond_Sorting/` | an abandoned third integration |
| `Deprecated/Non_Aligned/` | the pre-alignment era |
| `Subclustering/ImmuneCells/Deprecated/Mapping_1` … `Mapping_5` | iterations of the immune subclustering; the fifth produced the object the committed compartment scripts read, so that input has no committed producer (see `2_compartments_p14_p17/METHODS_CHOICES.md` item 9) |
| `Sequencing/Merging/TimeCourseOIR/MergeDGE.R` | upstream of this repository; described in prose in the root README |
| every `.rds`, `.txt` matrix, `.loom`, `.csv` output | data, never committed |
