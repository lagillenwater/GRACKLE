# Sourced by every GRACKLE Alpine job. Puts the micromamba environment on PATH: `module load
# anaconda` fails inside batch jobs (see fm-pm-evaluator's ralpine notes), and putting an
# environment's bin first on PATH is the pattern proven to work. Also sets the data and scratch
# roots read by analysis/hcmi_persistence/config.R.
export GRACKLE_ENV="/projects/$USER/software/envs/grackle"
export PATH="$GRACKLE_ENV/bin:$PATH"
export RETICULATE_PYTHON="$GRACKLE_ENV/bin/python"
export R_LIBS_USER="$GRACKLE_ENV/lib/R/library"
export GRACKLE_DATA="/projects/$USER/grackle_data/hcmi"
export GRACKLE_SCRATCH="/scratch/alpine/$USER/grackle/hcmi_persistence"
export XDG_CACHE_HOME="/scratch/alpine/$USER/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export TF_NUM_INTRAOP_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export TF_NUM_INTEROP_THREADS=1
export TF_CPP_MIN_LOG_LEVEL=3
mkdir -p logs "$XDG_CACHE_HOME"
