#!/usr/bin/env bash
# Measure what the full all-vs-all will cost, before committing to it.
#
# The cost of one Mash comparison is linear in the sketch size, so s=100000 does
# a hundred times the work of s=1000 for the same 4.6e8 pairs. That factor is
# what decides whether the largest sketch is run over the whole collection or
# over a stratified subsample, and it is measured here rather than guessed.
#
# A random subset of N genomes (seed 42) is sketched at the three sizes and
# compared all-vs-all. Since the number of pairs grows with N squared, the full
# run is extrapolated by (30209/N)^2. The same run also measures the fraction of
# pairs that survive the screening cutoff, which is what fixes the size of the
# retained-pair file.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

N=${N:-1500}
CUTOFF=${CUTOFF:-0.28}
BENCH=$SKETCH/bench
mkdir -p "$BENCH"

manifest=$DATA/genome_manifest.tsv
[[ -s $manifest ]] || { echo "ABORT: no manifest — run 05_build_sketches.sh" >&2; exit 1; }

python3 - "$manifest" "$N" "$BENCH/genome_list.txt" <<'PY'
import random, sys
src, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
paths = [l.split("\t")[1].rstrip("\n") for l in open(src).readlines()[1:]]
random.seed(42)
sample = sorted(random.sample(paths, n))
open(out, "w").write("\n".join(sample) + "\n")
print("subset: %d of %d genomes (seed 42)" % (len(sample), len(paths)))
PY

TOTAL_PAIRS=$(python3 -c "n=30209; print(n*(n-1)//2)")
SUB_PAIRS=$(python3 -c "n=$N; print(n*(n-1)//2)")
SCALE=$(python3 -c "print($TOTAL_PAIRS/$SUB_PAIRS)")
printf 'subset pairs %d, full pairs %d, scale factor %.1f\n\n' \
       "$SUB_PAIRS" "$TOTAL_PAIRS" "$SCALE"

printf '%-8s %10s %10s %12s %10s %14s %14s\n' \
       sketch sketch_s triangle_s edges frac_kept full_est_h full_edges_est

for s in 1000 10000 100000; do
    msh=$BENCH/sub_s${s}.msh
    if [[ ! -s $msh ]]; then
        t0=$(date +%s.%N)
        "$MASH" sketch -p "$THREADS" -k 21 -s "$s" -o "${msh%.msh}" \
                -l "$BENCH/genome_list.txt" > "$LOGS/bench_sketch_s${s}.log" 2>&1
        t_sketch=$(python3 -c "print('%.1f' % ($(date +%s.%N) - $t0))")
    else
        t_sketch=cached
    fi

    t0=$(date +%s.%N)
    edges=$("$MASH" triangle -E -d "$CUTOFF" -p "$THREADS" "$msh" 2>/dev/null | wc -l)
    t_tri=$(python3 -c "print('%.1f' % ($(date +%s.%N) - $t0))")

    python3 - "$s" "$t_sketch" "$t_tri" "$edges" "$SUB_PAIRS" "$SCALE" <<'PY'
import sys
s, t_sketch, t_tri, edges, sub_pairs, scale = sys.argv[1:]
t_tri, edges, sub_pairs, scale = float(t_tri), int(edges), int(sub_pairs), float(scale)
print('%-8s %10s %10.1f %12d %9.4f%% %14.2f %14.3e'
      % ('s=' + s, t_sketch, t_tri, edges, 100.0 * edges / sub_pairs,
         t_tri * scale / 3600.0, edges * scale))
PY
done

echo
echo "full_est_h assumes the subset scales as N^2 at constant threads."
echo "full_edges_est is what the accumulator will read, and sets the size of"
echo "the retained-pair file at roughly 40 bytes per edge before compression."
