#!/bin/bash
# Reproduced as submitted, with the author's email removed and the
# account given its required _cpu suffix. R and module versions are left
# as they were: they record the software this step actually ran under.
#SBATCH --ntasks=16
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=40G
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


#R CMD BATCH CellGSVA.part1.R CellGSVA.part1.Rout

#R CMD BATCH CellGSVA.part2.R CellGSVA.part2.Rout

module load scipy-stack

# the python 3.7 virtualenv this step ran under; point PY37_ENV at your own
source "${PY37_ENV:?set PY37_ENV to the activate script of a python 3.7 env}"

# cd <directory holding this script>

python run_GSVA.py



	
exit
