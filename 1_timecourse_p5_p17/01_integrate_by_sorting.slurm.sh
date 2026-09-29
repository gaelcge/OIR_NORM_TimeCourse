#!/bin/bash
# Wrapper for step 01. Reproduced as submitted, with the author's email
# removed and the account given its required _cpu suffix. The original
# `cd` pointed at Aligned/Cond_Sorting/SeuratV3, a path that does not
# exist - the scripts are under Aligned/SeuratV3 - and MakeSeurat.2.R was
# commented out, so the wrapper as it stood ran only half the step.
# Run this from the directory holding the step scripts.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:30:00
#SBATCH --mem=180G
#SBATCH --job-name=Seurat_TimeCourse
#SBATCH --account=def-jsjoyal_cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH -o
#SBATCH --no-requeue


###EdgeR-DGE requires -lmem=46gb for 3h
###MakeSeurat requires -lmem=80gb for 7h
#

# cd <directory holding this script>


R CMD BATCH MakeSeurat.R MakeSeurat.Rout

#R CMD BATCH MakeSeurat.2.R MakeSeurat.2.Rout

	
exit
