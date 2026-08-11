#!/usr/bin/env bash
# Sketch the proteomes for the AAI axis of D7.
#
# The k of a Mash sketch is not transferable between alphabets. With four bases
# and k = 21 the DNA axis has a measurable ceiling of 0.296 at s = 1000 (§3.2);
# with twenty amino acids the same formula, d_max = -ln(2j/(1+j))/k with
# j = 1/s, gives 0.888 at k = 7 and 0.691 at k = 9. That difference is not
# cosmetic: measured on a pilot of 100 proteomes, the share of inter-genus pairs
# that saturate — no shared hash at all, reported as distance 1 — is 1.6 % at
# k = 7, 16 % at k = 9 and 39 % at k = 11, against 91 % on the DNA axis. The
# protein axis is built at k = 7 for that reason, with k = 9 kept as a
# sensitivity stratum rather than as the default it usually is (D12).
#
# Everything else matches the DNA sketches so that the two axes stay comparable:
# same seed, same sketch sizes, one entry per accession, built from a manifest
# rather than by walking a directory (D10).
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

KS=${KS:-"7 9"}
SIZES=${SIZES:-"1000"}
manifest=$DATA/proteome_manifest.tsv
[[ -s $manifest ]] || { echo "ABORT: no proteome manifest — run 19_fetch_proteomes.sh" >&2; exit 1; }

tail -n +2 "$manifest" | cut -f2 > "$SKETCH/proteome_list.txt"
n=$(wc -l < "$SKETCH/proteome_list.txt")
echo "proteomes in the manifest: $n"

for k in $KS; do
    for s in $SIZES; do
        out=$SKETCH/prot_k${k}_s${s}
        [[ -s $out.msh ]] && { echo "  $out.msh exists, skipped"; continue; }
        echo "=== mash sketch -a -k $k -s $s  ($(date '+%F %T'))"
        "$MASH" sketch -a -k "$k" -s "$s" -S 42 -p "$THREADS" \
            -o "$out" -l "$SKETCH/proteome_list.txt" \
            > "$LOGS/sketch_prot_k${k}_s${s}.log" 2>&1
        got=$("$MASH" info -t "$out.msh" | tail -n +2 | grep -c . || true)
        echo "    $out.msh — $got sketches"
        [[ $got -eq $n ]] || { echo "ABORT: $got sketches for $n proteomes" >&2; exit 1; }
    done
done
ls -la "$SKETCH"/prot_k*.msh
