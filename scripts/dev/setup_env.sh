#!/usr/bin/env bash
# Builds the GRACKLE conda-forge environment with a standalone micromamba binary.
# The same script runs on the laptop and on Alpine (scripts/alpine/setup_env.sbatch).
#
#   scripts/dev/setup_env.sh <env-prefix> <micromamba-dir>
#
# Laptop:  scripts/dev/setup_env.sh "$HOME/.local/envs/grackle" "$HOME/.local/bin"
# Alpine:  scripts/dev/setup_env.sh "/projects/$USER/software/envs/grackle" "/projects/$USER/software/bin"
#
# Writes environment.lock.<platform>.txt (an explicit package list) into the repository root.
set -euo pipefail
PREFIX="${1:?env prefix required}"
MAMBA_DIR="${2:?micromamba directory required}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64) PLATFORM=osx-arm64 ;;
  Linux-x86_64) PLATFORM=linux-64 ;;
  *) echo "unsupported platform $(uname -s)-$(uname -m)" >&2; exit 1 ;;
esac
MAMBA="$MAMBA_DIR/micromamba"
if [[ ! -x "$MAMBA" ]]; then
  mkdir -p "$MAMBA_DIR"
  curl -Ls "https://micro.mamba.pm/api/micromamba/$PLATFORM/latest" \
    | tar -xj -C "$MAMBA_DIR" --strip-components=1 bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$(dirname "$PREFIX")/.mamba"
if [[ -d "$PREFIX/conda-meta" ]]; then
  "$MAMBA" install -y -p "$PREFIX" -f "$REPO/environment.yml"
else
  "$MAMBA" create -y -p "$PREFIX" -f "$REPO/environment.yml"
fi
# tidyestimate is CRAN-only; pin it.
"$PREFIX/bin/Rscript" -e 'if (!requireNamespace("tidyestimate", quietly = TRUE) || as.character(packageVersion("tidyestimate")) != "1.1.1") remotes::install_version("tidyestimate", version = "1.1.1", repos = "https://cloud.r-project.org", upgrade = "never"); stopifnot(requireNamespace("tidyestimate", quietly = TRUE))'
"$MAMBA" env export -p "$PREFIX" --explicit > "$REPO/environment.lock.$PLATFORM.txt"
RETICULATE_PYTHON="$PREFIX/bin/python" "$PREFIX/bin/Rscript" -e 'library(tensorflow); x <- tf$constant(1); stopifnot(as.numeric(x$numpy()) == 1); cat("tensorflow", tf$version$VERSION, "ok\n")'
echo "environment ready: $PREFIX ($PLATFORM)"
