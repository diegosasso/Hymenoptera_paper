#!/bin/bash
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -eo pipefail
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load "${PHYLO_SAMPLING_R_MODULE:-r-env/461}"
set -u
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "${PHYLO_SAMPLING_DIR:-$SLURM_SUBMIT_DIR}"
export R_PARALLEL_PORT=$(python3 scripts/select-worker-port.py)
# Each of 100 array elements runs the 15 conditions for one replicate number.
# This stays below the partition's 200 submitted-job limit.
for condition in $(seq 1 15); do
  task=$(( (condition - 1) * 100 + SLURM_ARRAY_TASK_ID ))
  Rscript --vanilla scripts/03-run-trial.R "$task"
done

