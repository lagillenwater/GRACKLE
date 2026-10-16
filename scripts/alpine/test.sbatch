#!/bin/bash
#SBATCH --job-name=grackle-test
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=7G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Runs the testthat suite on an Alpine compute node in the pinned environment. Set E2E=true to
# include the end-to-end fixture test.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"; cd "$REPO"; source scripts/alpine/env.sh
echo "commit=$(git rev-parse --short HEAD) job=$SLURM_JOB_ID"
if [[ "${E2E:-false}" == "true" ]]; then make e2e GRACKLE_ENV="$GRACKLE_ENV"; fi
make test GRACKLE_ENV="$GRACKLE_ENV"
