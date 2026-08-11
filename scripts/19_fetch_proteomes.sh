#!/usr/bin/env bash
# Fetch the proteome of every type genome, for the AAI axis of D7.
#
# The type material was downloaded from GenBank, and GenBank assemblies are not
# annotated by NCBI, so they carry no protein set. Their RefSeq counterparts are:
# 29,990 of the 30,209 assemblies (99.28 %) have a paired GCF accession, and
# RefSeq annotates prokaryotic genomes with PGAP, so `protein.faa` exists for
# them. That paired accession is already in the label table and is what this step
# resolves against — the alternative, reconstructing FTP paths from assembly
# names, breaks whenever the RefSeq assembly name differs from the GenBank one.
#
# Each proteome is stored under the **GenBank accession of its genome**, not
# under the RefSeq one. That single decision lets every downstream step —
# accumulator, leave-one-out, macro curve, conflict catalogue — run unchanged on
# the protein axis, because the keys are the same ones the label table uses.
#
# The 219 assemblies with no RefSeq counterpart are written to a separate list
# rather than silently dropped: they are a declared gap in the AAI axis, and the
# curves have to report on how many genomes they rest.
#
# Restartable: an accession whose file is already present is skipped, so an
# interrupted run resumes where it stopped instead of starting over.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

BATCH=${BATCH:-500}
FAA=$PROTEOME/faa
mkdir -p "$FAA" "$PROTEOME/tmp"

labels=$DATA/genome_labels.tsv
[[ -s $labels ]] || { echo "ABORT: no label table — run 03_build_labels.py" >&2; exit 1; }

# --- what to ask for, and what cannot be asked for ---------------------------
python3 - "$labels" "$PROTEOME" "$FAA" <<'PY'
import csv, os, sys
labels, prot, faa = sys.argv[1:4]
todo, nogcf = [], []
for r in csv.DictReader(open(labels), delimiter="\t"):
    gcf = (r["paired_gcf"] or "").strip()
    if not gcf.startswith("GCF"):
        nogcf.append((r["accession"], r["organism"]))
        continue
    if not os.path.exists(os.path.join(faa, r["accession"] + ".faa.gz")):
        todo.append((gcf, r["accession"]))
with open(os.path.join(prot, "pending.tsv"), "w") as fh:
    for gcf, acc in todo:
        fh.write("%s\t%s\n" % (gcf, acc))
with open(os.path.join(prot, "no_refseq.tsv"), "w") as fh:
    fh.write("accession\torganism\n")
    for acc, org in nogcf:
        fh.write("%s\t%s\n" % (acc, org))
print("to fetch: %d    already present: skipped    no RefSeq counterpart: %d"
      % (len(todo), len(nogcf)))
PY

n=$(wc -l < "$PROTEOME/pending.tsv")
[[ $n -gt 0 ]] || { echo "nothing to fetch"; }

# --- fetch in batches ---------------------------------------------------------
cut -f1 "$PROTEOME/pending.tsv" > "$PROTEOME/tmp/all_gcf.txt"
split -l "$BATCH" -d -a 4 "$PROTEOME/tmp/all_gcf.txt" "$PROTEOME/tmp/batch_"

i=0
total=$(ls "$PROTEOME"/tmp/batch_* 2>/dev/null | wc -l)
for b in "$PROTEOME"/tmp/batch_*; do
    i=$((i + 1))
    zip=$PROTEOME/tmp/b.zip
    rm -rf "$PROTEOME/tmp/x" "$zip"
    ok=0
    for attempt in 1 2 3; do
        if "$DATASETS" download genome accession --inputfile "$b" --include protein \
               --filename "$zip" --no-progressbar >/dev/null 2>&1; then
            ok=1; break
        fi
        sleep $((attempt * 5))
    done
    [[ $ok -eq 1 ]] || {
        echo "  batch $i/$total: three attempts failed, left for the next pass" >&2
        continue; }
    unzip -q -o "$zip" -d "$PROTEOME/tmp/x"

    # rename to the GenBank accession and compress; mash reads gzip directly
    python3 - "$PROTEOME/tmp/x/ncbi_dataset/data" "$PROTEOME/pending.tsv" "$FAA" <<'PY'
import csv, os, shutil, subprocess, sys
data, pending, faa = sys.argv[1:4]
gca = {}
for line in open(pending):
    g, a = line.split()
    gca[g] = a
todo = []
for gcf in os.listdir(data):
    src = os.path.join(data, gcf, "protein.faa")
    if gcf in gca and os.path.exists(src):
        dst = os.path.join(faa, gca[gcf] + ".faa")
        shutil.move(src, dst)
        todo.append(dst)
if todo:
    subprocess.run(["xargs", "-P", "8", "-n", "16", "gzip", "-1", "-f"],
                   input="\n".join(todo).encode(), check=True)
PY
    printf "  batch %d/%d  (%s files so far)\n" "$i" "$total" "$(ls "$FAA" | wc -l)"
done
rm -rf "$PROTEOME/tmp"

# --- manifest and reconciliation ---------------------------------------------
python3 - "$labels" "$FAA" "$DATA/proteome_manifest.tsv" <<'PY'
import csv, os, sys
labels, faa, out = sys.argv[1:4]
rows = list(csv.DictReader(open(labels), delimiter="\t"))
have, missing = [], []
for r in rows:
    p = os.path.join(faa, r["accession"] + ".faa.gz")
    (have if os.path.exists(p) else missing).append((r["accession"], p, r["paired_gcf"]))
with open(out, "w") as fh:
    fh.write("accession\tpath\tpaired_gcf\n")
    for acc, p, gcf in have:
        fh.write("%s\t%s\t%s\n" % (acc, p, gcf))
print("\nproteomes on disk: %d of %d genomes (%.2f %%)"
      % (len(have), len(rows), 100.0 * len(have) / len(rows)))
print("missing: %d" % len(missing))
PY
echo "manifest: $DATA/proteome_manifest.tsv"
du -sh "$FAA"
