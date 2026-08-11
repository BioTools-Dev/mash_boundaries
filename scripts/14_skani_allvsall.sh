#!/usr/bin/env bash
# The ANI axis: alignment-based identity for every pair skani can resolve.
#
# skani screens candidates with its own k-mer filter and aligns what survives, so
# no Mash-based prefilter is involved and neither metric decides what the other
# gets to see. The two-stage design of D6, with its selection-bias correction, is
# therefore unnecessary — measured on the benchmark subset, skani reports 100 %
# of the pairs Mash places below distance 0.15, and there is not one pair with an
# ANI that Mash failed to report.
#
# --min-af is set to 5 rather than the default 15. At 15 the tool withholds pairs
# whose alignment covers little of either genome, which is exactly the distant
# congeners the genus boundary is about: coverage of the 0.15-0.20 Mash band fell
# to 42 %. At 5 the band the study cares about is complete.
#
# Coverage ends near 80 % ANI regardless of settings. That is a property of
# alignment-based ANI, not of this design, and it is why the far tail of the Mash
# curves has no ANI counterpart and must not be reported as if it had one.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

SCREEN=${SCREEN:-80}
MINAF=${MINAF:-5}
T=${SKANI_THREADS:-$THREADS}
DB=$SKETCH/skani_db
OUT=$DIST/skani_ani.tsv

manifest=$DATA/genome_manifest.tsv
[[ -s $manifest ]] || { echo "ABORT: no manifest — run 05_build_sketches.sh" >&2; exit 1; }
tail -n +2 "$manifest" | cut -f2 > "$SKETCH/genome_list.txt"

if [[ ! -d $DB ]]; then
    echo "=== sketching with skani  ($(date '+%F %T'))"
    "$SKANI" sketch -t "$T" -l "$SKETCH/genome_list.txt" -o "$DB" \
        > "$LOGS/skani_sketch.log" 2>&1
    echo "    index $(du -sh "$DB" | cut -f1)"
fi

find "$DB" -name '*.sketch' | sort > "$SKETCH/skani_sketch_list.txt"
n=$(wc -l < "$SKETCH/skani_sketch_list.txt")
n_exp=$(($(wc -l < "$manifest") - 1))
echo "sketches: $n of $n_exp genomes"
[[ $n -eq $n_exp ]] || { echo "ABORT: skani indexed $n genomes, expected $n_exp" >&2; exit 1; }

echo "=== skani triangle -s $SCREEN --min-af $MINAF -t $T  ($(date '+%F %T'))"
"$SKANI" triangle -t "$T" -E -s "$SCREEN" --min-af "$MINAF" \
    -l "$SKETCH/skani_sketch_list.txt" -o "$OUT" \
    > "$LOGS/skani_triangle.log" 2>&1

echo "pairs with ANI: $(($(wc -l < "$OUT") - 1))"
echo "output: $OUT ($(du -h "$OUT" | cut -f1))"
zstd -q -3 -T4 -f "$OUT" -o "$OUT.zst" && rm -f "$OUT"
echo "compressed: $OUT.zst ($(du -h "$OUT.zst" | cut -f1))"
echo "done ($(date '+%F %T'))"
