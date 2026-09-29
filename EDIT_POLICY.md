# Edit policy

This repository is a reorganisation of analysis code that ran on the Digital
Research Alliance of Canada cluster Narval between 2018 and 2026. **Nothing was
re-run to produce it.** The results, objects and figures the manuscript reports
were produced by the code as it stood on the cluster, and this repository exists
to make that code readable — not to improve it.

That distinction sets a hard limit on what was allowed to change. Renaming a
variable, tidying a loop or "fixing" a threshold would break the correspondence
between the published numbers and the published code, and there is no way to
verify such a change without re-running an analysis that takes 8 h on 180 GB.
So the transformation was restricted to five classes of change, listed below,
and every edited file was diffed against its cluster original to confirm no
sixth class crept in.

## The five permitted changes

**1. Paths lifted into a section-0 parameter block.**
Every script opened with `setwd()` and read its inputs through absolute
`/home/gaelcge/projects/...` paths written inline, sometimes several times in
one file. Those are replaced by a `DATA_DIR` assignment in the parameter block
at the top, with paths built from it. The values are unchanged — only their
location in the file is. This is the one change that touches many lines, and it
is the change that makes the code runnable by anyone other than its author.

**2. The stale input path corrected.**
`MakeDGEmatrix.R` and `SeuratV3/.../MakeSeurat.R` both read the merged matrix
from `.../projects/def-jsjoyal/gaelcge/Sequencing/Merging/TimeCourseOIR/`. That
directory does not exist. The file is on the `ctb-jsjoyal` allocation:
`/project/ctb-jsjoyal/gaelcge/Sequencing/Merging/TimeCourseOIR/Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt`
(4.67 GB). The scripts as they stood on the cluster could not run; the data
evidently moved allocations at some point after they last ran. Corrected via
`DATA_DIR`.

**3. A personal email address removed.**
Every Slurm wrapper carried `#SBATCH --mail-user=<author's personal gmail>`.
Replaced with a placeholder. This changes no computation.

**4. The Slurm account string corrected.**
The wrappers request `--account=def-jsjoyal`. Alliance schedulers now require
the resource-suffixed form and reject the bare name when an account has
multiple allocations, so this is `--account=def-jsjoyal_cpu`. The wrappers are
kept as a record of the resources each step needed (`--mem=180G`,
`--time=08:30:00`, 16 cores), not as scripts expected to run unmodified on a
cluster that has changed twice since.

**5. One inert line deleted.**
`SeuratV3/.../MakeSeurat.2.R` assigned the clustering resolution twice: `res = 1`
in the parameter block, then `res = 0.8` on the line immediately before
`FindClusters`. The second assignment did not execute in the run that produced
the object — verified by reading the object's metadata, which carries the
clustering in a column named `integrated_snn_res.1`, and by the `1.20.30` in
every output filename. The line is removed because leaving it in a published
script would tell a reader the analysis used 0.8. Recorded in
`1_timecourse_p5_p17/METHODS_CHOICES.md`.

**6. Relic code commented out, never deleted.**
A decade of work in one directory means blocks pasted in from other projects.
373 lines across 60 blocks in 10 files of `2_compartments_p14_p17` are disabled
with a `#~ ` prefix and a banner saying why. Two things qualify, and only two —
each is a property the file itself proves, not a judgement about the science:

- **A `setwd()` into another project's tree.** Seven such calls point into
  `Retina/Retina_Mix/Test13_CD31_NORM_OIR_Sirt3/...` (a CD31-sorted Sirt3
  knockout experiment), `Retina/Retina_Rytvela/...`, a decommissioned
  `/RQexec/` path, or another user's home on a different allocation. Those
  directories no longer exist, and `setwd()` to a missing directory aborts an R
  script — so the line makes everything after it unreachable. Commenting it out
  leaves the working directory at the step's own output directory.
- **A statement operating on an object the file never loads.** Blocks
  referencing `VascularEndothelium`, `VascularProlEndothelium`,
  `VascularTipEndothelium`, `Seurat_object.P12_P14_P17_ECs`, `NeuroGlia` and
  similar in scripts that assign none of them. These were pasted from a live
  session that had the object in memory; run on their own they fail with
  "object not found". Several of them filter on
  `Dataset %in% c("OIR.P17.CD31.Sirt3KO.S129", ...)` — a metadata column this
  object does not have — which is what identifies the dataset they came from.

The heaviest case is
`compartments/02_neuroglial/02_de_by_condition.R`: 151 of its lines operate on
vascular-endothelium objects from the Sirt3 experiment.

Commenting rather than deleting was the author's instruction and is the right
call: the text stays readable and restorable, and a reader can see that the
disabled code was never part of this analysis rather than wondering what was
removed.

## What was deliberately NOT changed

- **Thresholds, filters and test choices**, including the differential-expression
  tables being filtered on the unadjusted p-value while the volcano plots colour
  on the adjusted one. This is documented in the relevant `METHODS_CHOICES.md`
  rather than silently reconciled.
- **The cell-type collapse repeated by hand in three scripts.** All ten pairs
  were checked and are identical in all three, so replacing them with a single
  committed CSV would have been behaviour-preserving. It was still not done: the
  CSV is committed alongside as `annotation/celltype_grouping.csv` for reading
  and checking, and the scripts keep their literals. The refactor is available
  but it is a change to code that produced published results, and it was not
  needed to make the repository readable.
- **`CRYG.genes` computed and never used** in the QC step. A harmless dead line
  in a script that produced published numbers; left in place.
- **Flat script structure.** These are not refactored into functions. A step is
  a flat script whose sections run in order, which is both the house convention
  and what these files already were.

## Provenance of every file

`PUBLISHED_TREE.md` maps every path in this repository to the absolute path of
the file it came from on Narval, so any file here can be diffed against the
original.
