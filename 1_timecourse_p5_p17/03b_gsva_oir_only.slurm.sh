#!/bin/bash
# Wrapper for step 03b. Reproduced as submitted (--ntasks=40,
# --mem-per-cpu=200G); the same caveat as step 03 applies.
#SBATCH --ntasks=40
#SBATCH --time=25:00:00
#SBATCH --mem-per-cpu=200G
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
