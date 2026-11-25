#!/bin/bash
#SBATCH --job-name=haika-bm
#SBATCH --account=project_2004996
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -eo pipefail
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load r-env/461
set -u
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /scratch/project_2004996/haika-wings
Rscript --vanilla scripts/04-fit-bm.R "${SLURM_ARRAY_TASK_ID}"
