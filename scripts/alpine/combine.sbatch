#!/bin/bash
#SBATCH --job-name=grackle-combine
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Runs 06_combine.R with ARM and STAGE (main | sensitivity | compare). Submitted with
# --dependency=afterany:<array>, so it runs even if some fits failed; those settings are excluded
# and listed in selection.csv.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${ARM:?ARM not set}" "${STAGE:?STAGE not set}" "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID arm=$ARM stage=$STAGE run=$GRACKLE_RUN_ID"
Rscript analysis/hcmi_persistence/06_combine.R --arm "$ARM" --stage "$STAGE" --cores "$SLURM_CPUS_PER_TASK"
