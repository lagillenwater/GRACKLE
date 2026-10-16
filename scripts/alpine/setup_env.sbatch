#!/bin/bash
#SBATCH --job-name=grackle-setup-env
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Builds /projects/$USER/software/envs/grackle from environment.yml with a standalone micromamba,
# then writes environment.lock.linux-64.txt into this checkout. Pull it back and commit it.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; mkdir -p logs
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
export XDG_CACHE_HOME="/scratch/alpine/$USER/cache"; mkdir -p "$XDG_CACHE_HOME"
scripts/dev/setup_env.sh "/projects/$USER/software/envs/grackle" "/projects/$USER/software/bin"
