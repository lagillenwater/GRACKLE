#!/bin/bash
#SBATCH --job-name=grackle-02-samples
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3840M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
Rscript analysis/hcmi_persistence/02_sample_table.R
