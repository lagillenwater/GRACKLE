#!/bin/bash
#SBATCH --job-name=grackle-align
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Runs 08_align.R on the pair-alignment fits of GRACKLE_RUN_ID, then the cluster-program heatmaps
# (09_cluster_programs.R). Submitted with
# --dependency=afterany:<align array>, so it runs even if some fits failed; their settings are
# marked ineligible in settings.csv.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
: "${GRACKLE_RUN_ID:?GRACKLE_RUN_ID not set}"
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID run=$GRACKLE_RUN_ID"
Rscript analysis/hcmi_persistence/08_align.R --cores "$SLURM_CPUS_PER_TASK"
Rscript analysis/hcmi_persistence/09_cluster_programs.R
