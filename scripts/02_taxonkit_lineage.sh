#!/usr/bin/env bash
# Resolve the standard-rank lineage of every assembly's taxid.
#
# Type material is indexed by STRAIN taxid, so the species and genus a pair of
# genomes must be compared on are not in the assembly report and have to be
# derived. taxonkit reformat does that against the pinned dump, and -t emits the
# taxid of each rank alongside its name. Ranks are joined by '|' rather than a
# tab so the names field and the taxids field stay one column each.
#
# Merged taxids are followed automatically; deleted ones come back empty and are
# counted here rather than discovered later.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

raw=$DATA/genomes_raw.tsv
out=$DATA/lineage.tsv
[[ -s $raw ]] || { echo "ABORT: $raw missing — run 01_extract_metadata.py first" >&2; exit 1; }

echo "taxdump: $TAXDUMP (pinned $TAXDUMP_DATE)"

# distinct taxids only; ~30 k assemblies collapse to fewer lookups
tail -n +2 "$raw" | cut -f2 | sort -u > "$DATA/taxids.txt"
echo "distinct taxids: $(wc -l < "$DATA/taxids.txt")"

{ printf 'taxid\tnames\ttaxids\n'
  $TAXONKIT reformat -I 1 -f '{d}|{p}|{c}|{o}|{f}|{g}|{s}' -t \
      --data-dir "$TAXDUMP" "$DATA/taxids.txt" 2> "$LOGS/taxonkit.log"
} > "$out"

echo "wrote $out ($(($(wc -l < "$out") - 1)) rows)"
empty=$(tail -n +2 "$out" | awk -F'\t' '$2 ~ /^\|*$/' | wc -l)
echo "taxids with no standard-rank lineage: $empty"
[[ -s $LOGS/taxonkit.log ]] && { echo "--- taxonkit warnings (first 10):"; head -10 "$LOGS/taxonkit.log"; }
exit 0
