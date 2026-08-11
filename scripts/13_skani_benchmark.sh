#!/usr/bin/env bash
# Measure what an all-vs-all with skani costs, and how large its index gets.
#
# skani screens candidate pairs with its own k-mer filter before aligning, so the
# ANI axis does not need the Mash-based prefilter that D6 originally proposed.
# That removes the selection-bias problem at its root: the two metrics screen
# independently, and neither one decides which pairs the other gets to see.
#
# Nothing below ~80 % ANI is reported, which is the tool's default screen and
# also the regime where alignment-based ANI stops being meaningful. That is a
# property of ANI rather than a limitation of the design, and it is why the far
# tail of the Mash curves has no ANI counterpart.
#
# The same subset of 1 500 genomes as step 06 is used, so the two benchmarks are
# directly comparable, and the extrapolation to the full collection is again by
# N squared.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

N_LIST=$SKETCH/bench/genome_list.txt
BENCH=$SKETCH/bench_skani
SCREEN=${SCREEN:-80}
MINAF=${MINAF:-15}
T=${SKANI_THREADS:-$THREADS}

[[ -s $N_LIST ]] || { echo "ABORT: run 06_triangle_benchmark.sh first" >&2; exit 1; }
N=$(wc -l < "$N_LIST")
rm -rf "$BENCH"; mkdir -p "$BENCH"

echo "subset: $N genomes, screen -s $SCREEN, --min-af $MINAF, $T threads"
[[ $(awk '{print int($1)}' /proc/loadavg) -gt 2 ]] && \
    echo "NOTE: load average is $(cut -d' ' -f1 /proc/loadavg); timings are an upper bound"

t0=$(date +%s.%N)
"$SKANI" sketch -t "$T" -l "$N_LIST" -o "$BENCH/db" > "$LOGS/skani_sketch_bench.log" 2>&1
t_sketch=$(python3 -c "print('%.1f' % ($(date +%s.%N) - $t0))")
sz=$(du -sm "$BENCH/db" | cut -f1)

find "$BENCH/db" -name '*.sketch' | sort > "$BENCH/sketch_list.txt"
n_sk=$(wc -l < "$BENCH/sketch_list.txt")

t0=$(date +%s.%N)
"$SKANI" triangle -t "$T" -E -s "$SCREEN" --min-af "$MINAF" \
    -l "$BENCH/sketch_list.txt" -o "$BENCH/ani.tsv" > "$LOGS/skani_tri_bench.log" 2>&1
t_tri=$(python3 -c "print('%.1f' % ($(date +%s.%N) - $t0))")
edges=$(($(wc -l < "$BENCH/ani.tsv") - 1))

python3 - "$N" "$t_sketch" "$sz" "$n_sk" "$t_tri" "$edges" <<'PY'
import sys
N, t_sketch, sz, n_sk, t_tri, edges = sys.argv[1:]
N, sz, n_sk, edges = int(N), int(sz), int(n_sk), int(edges)
t_sketch, t_tri = float(t_sketch), float(t_tri)
FULL = 30209
scale = (FULL * (FULL - 1) / 2) / (N * (N - 1) / 2)
print()
print("sketching : %.1f s for %d genomes, index %d MB" % (t_sketch, n_sk, sz))
print("            full: ~%.0f min, ~%.1f GB" % (t_sketch * FULL / N / 60, sz * FULL / N / 1024))
print("triangle  : %.1f s, %d pairs reported of %d (%.3f %%)"
      % (t_tri, edges, N * (N - 1) // 2, 100 * edges / (N * (N - 1) / 2)))
print("            full: ~%.2f h, ~%.2e pairs" % (t_tri * scale / 3600, edges * scale))
PY

echo
echo "--- sample of the output:"
head -3 "$BENCH/ani.tsv" | cut -f1-5
