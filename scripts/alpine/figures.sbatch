#!/bin/bash
#SBATCH --job-name=grackle-figures
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=7G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID run=$GRACKLE_RUN_ID"
Rscript analysis/hcmi_persistence/07_figures.R
