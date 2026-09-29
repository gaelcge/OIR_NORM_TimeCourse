#!/bin/bash
# Reproduced as submitted, with the author's email removed and the
# account given its required _cpu suffix. R and module versions are left
# as they were: they record the software this step actually ran under.
#SBATCH -n 16
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=80G
#SBATCH --job-name=GSVA_py
#SBATCH --account=def-jsjoyal_cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH -o Seurat.report
#SBATCH --no-requeue

####
###EdgeR-DGE requires -lmem=46gb for 3h
###MakeSeurat requires -lmem=80gb for 7h
#

# cd <directory holding this script>

#R CMD BATCH CellGSVA.part1.R CellGSVA.part1.Rout

#R CMD BATCH CellGSVA.part2.R CellGSVA.part2.Rout

module load python/3.7.0

source ~/ENV_python3.7.0/bin/activate

module load scipy-stack

python run_GSVA.py

deactivate

	
exit
