#!/bin/bash
#SBATCH --job-name=haika-bounds
#SBATCH --account=project_2004996
#SBATCH --partition=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -eo pipefail
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load r-env/461
set -u
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /scratch/project_2004996/haika-wings
Rscript --vanilla scripts/09-bound-sensitivity.R
