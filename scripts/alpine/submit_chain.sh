#!/usr/bin/env bash
# Submits the HCMI persistence pipeline as a Slurm dependency chain through ralpine.
# Run it from the repository root on the laptop. Each dependency is read back with
# `ralpine jobinfo`, because a chain that silently drops a dependency runs out of order.
#
#   scripts/alpine/submit_chain.sh --from prep       # 02 -> 03 -> 03b -> 04 (data already downloaded)
#   scripts/alpine/submit_chain.sh --from fit        # fit arrays, combine, sensitivity, transfer, compare, figures
#   RUN_ID=<id> scripts/alpine/submit_chain.sh --from compare   # compare and figures for an existing run
set -euo pipefail
RALPINE="$(dirname "${BASH_SOURCE[0]}")/ralpine"
[[ "${1:-}" == "--from" ]] || { echo "usage: $0 --from prep|fit|compare" >&2; exit 1; }
FROM="${2:?--from needs prep, fit or compare}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%d)_$(git rev-parse --short HEAD)}"
EXP="--export=ALL,GRACKLE_RUN_ID=$RUN_ID"
ARMS=(coess_2d coess_3d string_mech string_all)

job_id() { grep -oE '[0-9]+$' <<<"$1" | tail -1; }
submit() { local out; out="$("$RALPINE" submit "$@")"; echo "$out" >&2; job_id "$out"; }
check_dependency() {
  local id="$1" want="$2" dep
  dep="$("$RALPINE" jobinfo "$id" | grep -oE 'Dependency=[^ ]+' || true)"
  if [[ "$dep" != *"$want"* ]]; then
    echo "job $id was submitted without its dependency ($dep, wanted $want); cancel it with: $RALPINE cancel $id" >&2
    exit 1
  fi
}
echo "run_id: $RUN_ID"

if [[ "$FROM" == "prep" ]]; then
  J2=$(submit scripts/alpine/02_sample_table.sbatch)
  J3=$(submit scripts/alpine/03_expression.sbatch "--dependency=afterok:$J2"); check_dependency "$J3" "afterok:$J2"
  J3b=$(submit scripts/alpine/03b_gene_graphs.sbatch "--dependency=afterok:$J3"); check_dependency "$J3b" "afterok:$J3"
  J4=$(submit scripts/alpine/04_rank.sbatch "--dependency=afterok:$J3"); check_dependency "$J4" "afterok:$J3"
  echo "prep chain: 02=$J2 03=$J3 03b=$J3b 04=$J4"
  exit 0
fi

if [[ "$FROM" == "fit" ]]; then
  WAIT=()
  for ARM in "${ARMS[@]}"; do   # four arms, each its own chain
    F=$(submit scripts/alpine/fit_array.sbatch --array=0-215%50 "$EXP,ARM=$ARM,KIND=main")
    C=$(submit scripts/alpine/combine.sbatch "--dependency=afterany:$F" "$EXP,ARM=$ARM,STAGE=main"); check_dependency "$C" "afterany:$F"
    S=$(submit scripts/alpine/fit_array.sbatch --array=0-11 "--dependency=afterok:$C" "$EXP,ARM=$ARM,KIND=sensitivity"); check_dependency "$S" "afterok:$C"
    T=$(submit scripts/alpine/fit_array.sbatch --array=0-9 "--dependency=afterok:$C" "$EXP,ARM=$ARM,KIND=transfer"); check_dependency "$T" "afterok:$C"
    CS=$(submit scripts/alpine/combine.sbatch "--dependency=afterany:$S" "$EXP,ARM=$ARM,STAGE=sensitivity"); check_dependency "$CS" "afterany:$S"
    echo "$ARM: fits=$F combine=$C sensitivity=$S transfer=$T combine_sensitivity=$CS"
    WAIT+=("$C" "$T" "$CS")
  done
  DEP="afterany:$(IFS=:; echo "${WAIT[*]}")"
fi

if [[ "$FROM" == "fit" || "$FROM" == "compare" ]]; then
  CMP=$(submit scripts/alpine/combine.sbatch ${DEP:+"--dependency=$DEP"} "$EXP,ARM=all,STAGE=compare")
  [[ -n "${DEP:-}" ]] && check_dependency "$CMP" "afterany:"
  FIG=$(submit scripts/alpine/figures.sbatch "--dependency=afterok:$CMP" "$EXP"); check_dependency "$FIG" "afterok:$CMP"
  echo "compare=$CMP figures=$FIG"
fi
echo "watch with: $RALPINE sq"
