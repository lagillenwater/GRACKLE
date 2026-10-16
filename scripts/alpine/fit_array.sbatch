#!/bin/bash
#SBATCH --job-name=grackle-fit
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.out
#
# One cross-fitting task per array index. ARM and KIND (main | sensitivity | transfer) choose the
# task table; GRACKLE_RUN_ID names the run. The array size is set at submission: --array=0-215%50
# for main, 0-11 for sensitivity, 0-9 for transfer. A completed task exits at once, so a failed
# index is resubmitted alone with --array=<index>.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${ARM:?ARM not set}" "${KIND:?KIND not set}" "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} arm=$ARM kind=$KIND run=$GRACKLE_RUN_ID"
if [[ "$KIND" == "main" ]]; then
  TASKS="analysis/hcmi_persistence/tasks/main_${ARM}.csv"
else
  TASKS="$GRACKLE_SCRATCH/$GRACKLE_RUN_ID/tasks/${KIND}_${ARM}.csv"
fi
Rscript analysis/hcmi_persistence/05_fit.R --tasks "$TASKS" --index "$((SLURM_ARRAY_TASK_ID + 1))" ${ITERATIONS:+--iterations "$ITERATIONS"}
