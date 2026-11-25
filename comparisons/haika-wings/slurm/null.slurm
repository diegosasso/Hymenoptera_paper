#!/bin/bash
#SBATCH --job-name=haika-null
#SBATCH --account=project_2004996
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
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
first_replicate=$(( (SLURM_ARRAY_TASK_ID - 1) * 10 + 1 ))
last_replicate=$(( SLURM_ARRAY_TASK_ID * 10 ))
if [ "$last_replicate" -gt 199 ]; then last_replicate=199; fi
failed=0
for replicate in $(seq "$first_replicate" "$last_replicate"); do
  if ! Rscript --vanilla scripts/06-null-comparison.R "$replicate"; then
    failed=1
  fi
done
exit "$failed"
