#!/bin/bash
# Reproduced as submitted, with the author's email removed and the
# account given its required _cpu suffix. R and module versions are left
# as they were: they record the software this step actually ran under.
#SBATCH -n 16
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=80G
#SBATCH --job-name=GSVA_R
#SBATCH --account=def-jsjoyal_cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH -o Run.report
#SBATCH --no-requeue

####
###EdgeR-DGE requires -lmem=46gb for 3h
###MakeSeurat requires -lmem=80gb for 7h
#

# cd <directory holding this script>

R CMD BATCH run_gsva.R run_gsva.Rout

	
exit
