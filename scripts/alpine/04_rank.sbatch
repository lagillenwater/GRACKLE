#!/bin/bash
#SBATCH --job-name=grackle-04-rank
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/04_rank.R
