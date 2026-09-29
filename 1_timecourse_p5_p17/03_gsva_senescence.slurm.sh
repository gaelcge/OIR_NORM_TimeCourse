#!/bin/bash
# Wrapper for step 03. Reproduced as submitted: --ntasks=26 with
# --mem-per-cpu=150G requests ~3.9 TB in aggregate, which is very unlikely
# to be what was meant. Kept as the record of what was submitted; size it
# yourself before re-running.
#SBATCH --ntasks=26
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=150G
#SBATCH --job-name=GSVA_R
#SBATCH --account=def-jsjoyal_cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH -o Run.report
#SBATCH --no-requeue


# cd <directory holding this script>

module load nixpkgs/16.09  gcc/7.3.0
module load r/3.6

R CMD BATCH run_gsva.R run_gsva.Rout
	
exit
