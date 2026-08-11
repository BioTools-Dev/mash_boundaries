#!/usr/bin/env bash
# Build the Mash sketches this study compares, from a canonical genome manifest.
#
# The manifest exists because the download left more than one FASTA in two of the
# assembly directories: NCBI datasets also wrote *_cds_from_genomic.fna and
# *_rna_from_genomic.fna, which are coding-sequence and RNA extracts, not
# genomes. A sketch built by globbing the tree therefore carries four entries
# that are not genomes at all — the prebuilt sketches in the source directory do
# carry them, which is how this was found. Selecting *_genomic.fna and excluding
# those two suffixes gives exactly one genome per accession.
#
# All sizes use k=21 and Mash's default seed 42, matching the prebuilt sketches
# so that any difference between them is the sketch size and nothing else.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

manifest=$DATA/genome_manifest.tsv
labels=$DATA/genome_labels.tsv

# --- manifest ----------------------------------------------------------------
if [[ ! -s $manifest ]]; then
    echo "building manifest from $GENOME_DIR"
    find "$GENOME_DIR" -name '*_genomic.fna' \
         ! -name '*_cds_from_genomic.fna' ! -name '*_rna_from_genomic.fna' \
      | sort > "$DATA/genome_paths.txt"
    python3 - "$DATA/genome_paths.txt" "$manifest" <<'PY'
import re, sys
src, out = sys.argv[1], sys.argv[2]
seen = {}
with open(out, "w") as fh:
    fh.write("accession\tpath\n")
    for line in open(src):
        p = line.rstrip("\n")
        m = re.match(r"(GC[AF]_\d+\.\d+)_", p.rsplit("/", 1)[-1])
        if not m:
            sys.exit("ABORT: cannot parse an accession from %s" % p)
        acc = m.group(1)
        if acc in seen:
            sys.exit("ABORT: %s has more than one genomic FASTA:\n  %s\n  %s"
                     % (acc, seen[acc], p))
        seen[acc] = p
        fh.write("%s\t%s\n" % (acc, p))
print("manifest: %d accessions" % len(seen))
PY
    rm -f "$DATA/genome_paths.txt"
fi

n_man=$(($(wc -l < "$manifest") - 1))
n_lab=$(($(wc -l < "$labels") - 1))
echo "manifest $n_man genomes / labels $n_lab genomes"
[[ $n_man -eq $n_lab ]] || { echo "ABORT: manifest and labels disagree" >&2; exit 1; }

# --- sketches ----------------------------------------------------------------
tail -n +2 "$manifest" | cut -f2 > "$SKETCH/genome_list.txt"

for s in 1000 10000 100000; do
    out=$SKETCH/type_s${s}
    if [[ -s ${out}.msh ]]; then
        echo "s=$s already built, skipping"
        continue
    fi
    echo "=== sketching s=$s  ($(date '+%F %T'))"
    /usr/bin/time -v "$MASH" sketch -p "$THREADS" -k 21 -s "$s" \
        -o "$out" -l "$SKETCH/genome_list.txt" \
        > "$LOGS/sketch_s${s}.log" 2>&1 || { echo "ABORT: mash sketch s=$s failed"; tail -20 "$LOGS/sketch_s${s}.log"; exit 1; }
    echo "    $(du -h "${out}.msh" | cut -f1)  $("$MASH" info -H "${out}.msh" | awk '/Sketches:/{print $2" sketches"}')"
done

echo "done ($(date '+%F %T'))"
