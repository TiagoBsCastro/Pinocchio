#!/usr/bin/env bash
# Automated benchmark for MASS_MAPS solver variants (Brent vs Bisection)
# Usage (from src directory):
#   ../scripts/benchmark_mass_maps.sh <paramfile> "mpirun -np 4" [--decide] [--threshold 0.03]
#
# ENV / Flags:
#   MASS_MAPS_NUM_RUNS   Number of repetitions (default 3)
#   MASS_MAPS_SOLVERS    Space separated list (default "BRENT BISECTION")
#   --decide             Emit recommendation & optional patch suggestion
#   --threshold <frac>   Relative total-time advantage needed to keep slower path (default 0.03 = 3%)
#
# Output:
#   profiling_<SOLVER>.txt  concatenated profiling sections per run
#   profiling_summary.tsv   tab-separated summary of averages
#   Recommendation printed if --decide specified.

set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <paramfile> <mpirun command (quoted)> [--decide] [--threshold 0.03]" >&2
  exit 1
fi

PARAMFILE=$1; shift
MPICMD=$1; shift
DECIDE=0
THRESH=0.03
while [ $# -gt 0 ]; do
  case "$1" in
    --decide) DECIDE=1; shift ;;
    --threshold) THRESH=$2; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

REPS=${MASS_MAPS_NUM_RUNS:-3}
SOLVERS=${MASS_MAPS_SOLVERS:-"BRENT BISECTION"}

if [ ! -f "$PARAMFILE" ]; then
  echo "Parameter file not found: $PARAMFILE" >&2
  exit 3
fi

run_case() { # solver label
  local solver=$1
  echo "==== Solver: ${solver} ===="
  make clean >/dev/null 2>&1 || true
  if [ "$solver" = "BISECTION" ]; then
    make pinocchio MASS_MAPS_SOLVER=BISECTION >/dev/null
  else
    make pinocchio >/dev/null
  fi
  local total=0
  local outfile="profiling_${solver}.txt"
  : > "$outfile"
  for r in $(seq 1 $REPS); do
    echo "-- Run $r/$REPS" >&2
    $MPICMD ./pinocchio.x $PARAMFILE 2>&1 | tee run_${solver}_$r.log | awk '/Mass maps:/{flag=1} flag{print}' >> "$outfile" || true
  done
  echo "Captured profiling in $outfile"
  echo
}

for s in $SOLVERS; do
  run_case "$s"
done

echo "Summary (average of Total and key buckets)" | tee profiling_summary.tsv
awk '/Total \(avg\)/{t=$3} /Scan condition/ {scan=$3} /Solver condition/ {solc=$3} /Solver overhead/ {solo=$3} /Bracket:/ {br=$3} /Counts:/{print FILENAME, t, scan, solc, solo, br}' profiling_*.txt | \
  awk -v out=profiling_summary.tsv 'BEGIN{hdr="solver\tTotal\tScan\tSolCond\tSolOvh\tBracket";print hdr >> out;printf("%-18s %8s %8s %8s %8s %8s\n","solver","Total","Scan","SolCond","SolOvh","Bracket");}
  {cnt[$1]++;T[$1]+=$2;SC[$1]+=$3;SO[$1]+=$4;SH[$1]+=$5;BR[$1]+=$6;}
  END{for(f in T){at=T[f]/cnt[f];asc=SC[f]/cnt[f];aso=SO[f]/cnt[f];ash=SH[f]/cnt[f];abr=BR[f]/cnt[f];printf("%-18s %8.4f %8.4f %8.4f %8.4f %8.4f\n",f,at,asc,aso,ash,abr);
      gsub("profiling_","",f); gsub(".txt","",f);printf("%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",f,at,asc,aso,ash,abr)>>out;}}'

if [ $DECIDE -eq 1 ]; then
  # Extract totals from TSV (expect exactly two rows, else just skip decision)
  mapfile -t rows < <(tail -n +2 profiling_summary.tsv)
  if [ ${#rows[@]} -ge 2 ]; then
    # Pick the smallest total
    bestSolver=""; bestTotal=1e99
    declare -A totals
    for r in "${rows[@]}"; do
      s=$(echo "$r" | awk '{print $1}')
      tot=$(echo "$r" | awk '{print $2}')
      totals[$s]=$tot
      awk -v s=$s -v tot=$tot -v best=$bestTotal 'END{}' >/dev/null
      if awk -v a=$tot -v b=$bestTotal 'BEGIN{exit !(a<b)}'; then
        bestSolver=$s; bestTotal=$tot
      fi
    done
    echo "Decision threshold (relative): $THRESH"
    for s in "${!totals[@]}"; do
      echo "  Total[$s] = ${totals[$s]}"
    done
    # Compare best against others
    keepBoth=0
    for s in "${!totals[@]}"; do
      [ "$s" = "$bestSolver" ] && continue
      rel=$(awk -v a=${totals[$s]} -v b=$bestTotal 'BEGIN{print (a/b)-1.0}')
      echo "  Relative overhead of $s vs $bestSolver: $(printf '%.2f' $(echo $rel))"
      if awk -v r=$rel -v th=$THRESH 'BEGIN{exit !(r<th)}'; then
        : # fine
      else
        keepBoth=1
      fi
    done
    if [ $keepBoth -eq 0 ]; then
      echo "Recommendation: keep ONLY solver '$bestSolver' (others slower by > threshold)." | tee -a profiling_summary.tsv
      echo "To remove the alternative from code: delete its #ifdef path and simplify build (manual step)." | tee -a profiling_summary.tsv
    else
      echo "Recommendation: performance differences < threshold; you may drop either for simplicity (default: keep BRENT)." | tee -a profiling_summary.tsv
    fi
  else
    echo "Decision: insufficient data (need >=2 solvers)." | tee -a profiling_summary.tsv
  fi
fi

if [ $DECIDE -eq 0 ]; then
  cat <<EOF
Next steps:
- Inspect profiling_*.txt for detailed bucket breakdowns.
- (Optional) re-run with --decide to print a recommendation.
EOF
fi
