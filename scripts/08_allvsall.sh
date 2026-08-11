#!/usr/bin/env bash
# Run the all-vs-all and accumulate it in one streaming pass per sketch size.
#
# `mash triangle` walks the lower triangle, so each of the 4.6e8 pairs is
# compared once, and -d filters the edge list at the source: only pairs at or
# below the screening cutoff are ever written, piped straight into the
# accumulator. No intermediate file holds the raw comparison.
#
# The cutoff defaults to 0.28 because that is just under the largest distance a
# sketch of 1000 hashes can represent — one shared hash out of a thousand gives
# d = 0.296, and anything less related shares no hashes at all and is reported
# as 1. Running the three sizes over the same window is what makes them
# comparable. Sizes above s=1000 are additionally run at a wider cutoff, since
# what they see beyond the small sketch's ceiling is itself a result.
#
# Usage:  bash scripts/08_allvsall.sh [sketch_size ...]     (default: 1000)
#         MSH=<sketch> TAG=<name> CUTOFF=<d> bash scripts/08_allvsall.sh 1
#         for an axis whose sketches are not named by sketch size.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

CUTOFF=${CUTOFF:-0.28}
KEEP=${KEEP:-1}                  # write the retained-pair file
sizes=("${@:-1000}")

# The protein axis of D12 runs through this same step: MSH names a sketch built
# by 20_build_protein_sketches.sh and TAG names its outputs, and everything
# downstream is unchanged because the proteome files are keyed by the GenBank
# accession of their genome.
for s in "${sizes[@]}"; do
    msh=${MSH:-$SKETCH/type_s${s}.msh}
    [[ -s $msh ]] || { echo "ABORT: $msh missing — run 05_build_sketches.sh" >&2; exit 1; }

    tag=${TAG:-s${s}}_d${CUTOFF}
    prefix=$DIST/$tag
    keep_arg=()
    [[ $KEEP == 1 ]] && keep_arg=(--keep-pairs "$prefix"_pairs.tsv.zst)

    echo "=== all-vs-all s=$s cutoff=$CUTOFF  ($(date '+%F %T'))"
    set +e
    "$MASH" triangle -E -d "$CUTOFF" -p "$THREADS" "$msh" 2> "$LOGS/triangle_$tag.log" \
      | python3 scripts/07_accumulate.py \
            --labels "$DATA/genome_labels.tsv" \
            --prefix "$prefix" "${keep_arg[@]}" 2>&1 | tee "$LOGS/accumulate_$tag.log"
    rc=${PIPESTATUS[0]}
    set -e
    [[ $rc -eq 0 ]] || { echo "ABORT: mash triangle failed (rc=$rc)"; tail -5 "$LOGS/triangle_$tag.log"; exit 1; }

    echo "--- outputs:"
    ls -la "$prefix"_* 2>/dev/null | sed 's/^/    /'
    echo "done s=$s ($(date '+%F %T'))"
    echo
done
