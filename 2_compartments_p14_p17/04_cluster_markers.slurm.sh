#!/bin/bash
# Reproduced as submitted, with the author's email removed and the
# account given its required _cpu suffix. R and module versions are left
# as they were: they record the software this step actually ran under.
#SBATCH -N 1
#SBATCH -n 16
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


module load r/3.5.0 

R CMD BATCH FindClustermarkers.R FindClustermarkers.Rout

	
exit
