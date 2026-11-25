#!/bin/bash
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -eo pipefail
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load "${NEFF_R_MODULE:-r-env/461}"
set -u
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "${NEFF_DIR:-$SLURM_SUBMIT_DIR}"
Rscript --vanilla scripts/check-completion.R 1
Rscript --vanilla scripts/04-audit.R 1
