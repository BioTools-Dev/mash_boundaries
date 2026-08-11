#!/usr/bin/env python3
"""Simulated metagenome bins with controlled completeness and contamination.

The measurements up to here take a complete genome as the query. A metagenome
bin is not that, and the boundary a bin sees is a different boundary. D13 sets
out why the bins are simulated from the type material itself rather than taken
from a real assembly: a bin built here carries its source accession, so the
answer is known exactly, whereas a real bin would need its genus assigned by the
same kind of judgement the whole study is trying to measure.

The construction is deliberately simple, because the quantity that matters to a
sketch is which k-mers survive and not how they are packaged. Each genome is cut
into contigs of a fixed length, contigs are dropped at random until the retained
fraction reaches the target completeness, and contigs from a second genome are
added until the foreign fraction reaches the target contamination. Fragmentation
by itself costs k-1 k-mers per break — about 0.1 % at 25 kb contigs and k = 21 —
so it is expected to be irrelevant; `--contig-len` exists so that expectation can
be tested rather than assumed.

Two contaminant classes are generated, because they are predicted to act in
opposite directions. A contaminant from another phylum adds k-mers that match
nothing in the reference and only dilutes the query. A congeneric contaminant
adds k-mers that match the right genus, so it should hurt the species call and
may help the genus call — and it is the realistic case, since what co-occurs in
a sample and confuses a binner is usually a relative.

Every query is its own control: completeness 1.0 with no contamination is a row
of the grid, so degradation is measured against the intact genome of the same
accession rather than against a collection-wide average.

The bins are written in batches, sketched, and deleted; only the sketch and the
manifest survive. The manifest records the bases retained and their origin,
which is what the containment arm of `scripts/24_bin_curves.py` needs.
"""
import argparse
import collections
import csv
import gzip
import os
import random
import subprocess
import sys


def read_labels(path):
    with open(path) as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_manifest(path):
    with open(path) as fh:
        return {r["accession"]: r["path"] for r in csv.DictReader(fh, delimiter="\t")}


def load_contigs(path, contig_len):
    """The genome cut into fixed-length pieces, as (header-less) sequences."""
    seqs, cur = [], []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    out = []
    for s in seqs:
        for i in range(0, len(s), contig_len):
            piece = s[i:i + contig_len]
            if len(piece) >= 200:          # a fragment shorter than this carries
                out.append(piece)          # almost no k-mers and no binner keeps it
    return out


def take(contigs, target_bases, rng):
    """Contigs drawn at random until the target number of bases is reached."""
    order = list(range(len(contigs)))
    rng.shuffle(order)
    picked, total = [], 0
    for i in order:
        if total >= target_bases:
            break
        picked.append(contigs[i])
        total += len(contigs[i])
    return picked, total


def pick_queries(lab, have, n, per_genus, seed):
    """Query genomes spread over genera, so no genus dominates the grid.

    Only genera holding more than one species are drawn from. Under
    `novel_species` the query's own species leaves the database, so in a
    monospecific genus a correct genus call is impossible before any degradation
    is applied, and including such queries would measure the composition of the
    sample rather than the effect of the bin. The restriction makes the absolute
    precision here higher than the collection-wide figure of §3.4 and is not
    comparable to it; what is comparable is each bin against its own intact
    genome, which is why completeness 1.0 is a row of the grid.
    """
    rng = random.Random(seed)
    by_genus = collections.defaultdict(dict)
    for acc in have:
        r = lab[acc]
        if r["low_quality"] == "0" and r["no_genus"] == "0":
            by_genus[r["genus"]].setdefault(r["species_taxid"], []).append(acc)
    genera = sorted(g for g, sp in by_genus.items() if len(sp) > 1)
    rng.shuffle(genera)
    picked = []
    for g in genera:
        accs = sorted(a for v in by_genus[g].values() for a in v)
        rng.shuffle(accs)
        picked.extend(accs[:per_genus])
        if len(picked) >= n:
            break
    return picked[:n]


def pick_contaminants(lab, have, queries, seed):
    """For each query, a congeneric donor and a donor from another phylum."""
    rng = random.Random(seed + 1)
    by_genus, by_phylum = collections.defaultdict(list), collections.defaultdict(list)
    for acc in have:
        r = lab[acc]
        if r["low_quality"] == "0":
            by_genus[r["genus"]].append(acc)
            by_phylum[r["phylum"]].append(acc)
    phyla = sorted(p for p in by_phylum if len(by_phylum[p]) >= 20)
    out = {}
    for q in queries:
        r = lab[q]
        # congeneric: same genus, different species. Absent for a monospecific
        # genus, and that is recorded rather than substituted.
        near = [a for a in by_genus[r["genus"]]
                if lab[a]["species_taxid"] != r["species_taxid"]]
        far_phylum = rng.choice([p for p in phyla if p != r["phylum"]])
        out[q] = (rng.choice(near) if near else "",
                  rng.choice(by_phylum[far_phylum]))
    return out


def write_bin(fh, name, native, foreign):
    fh.write(">%s\n" % name)
    for i, s in enumerate(native + foreign):
        if i:
            fh.write(">%s_c%d\n" % (name, i))
        fh.write(s + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True, help="output prefix for .msh and manifest")
    ap.add_argument("--work", required=True, help="scratch directory for the batches")
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--sketch-size", type=int, default=10000)
    ap.add_argument("--kmer", type=int, default=21)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--queries", type=int, default=400)
    ap.add_argument("--per-genus", type=int, default=2,
                    help="cap per genus, which is what keeps Streptomyces from "
                         "supplying a tenth of the grid")
    ap.add_argument("--contig-len", type=int, default=25000)
    ap.add_argument("--completeness", default="1.0,0.9,0.7,0.5,0.3,0.2,0.1")
    ap.add_argument("--contamination", default="0,0.10,0.20")
    ap.add_argument("--batch", type=int, default=500,
                    help="bins written before sketching and deleting them")
    a = ap.parse_args()

    lab = read_labels(a.labels)
    man = read_manifest(a.manifest)
    have = sorted(set(lab) & set(man))
    queries = pick_queries(lab, have, a.queries, a.per_genus, a.seed)
    donors = pick_contaminants(lab, have, queries, a.seed)
    comps = [float(x) for x in a.completeness.split(",")]
    conts = [float(x) for x in a.contamination.split(",")]
    print("queries: %d over %d genera"
          % (len(queries), len({lab[q]["genus"] for q in queries})))
    print("grid: %d completeness x %d contamination x {near,far} donors"
          % (len(comps), len(conts)))

    os.makedirs(a.work, exist_ok=True)
    rows, batch, parts = [], [], []
    rng = random.Random(a.seed + 2)

    def flush():
        """Sketch the batch, then remove it: only the sketch is worth keeping.

        Mash records the path it was handed as the sketch ID, so the batch is
        sketched from inside the working directory and passed bare file names:
        otherwise every bin would carry an absolute scratch path into a file
        meant to be published.
        """
        if not batch:
            return
        part = "part_%03d.msh" % len(parts)
        subprocess.run([a.mash, "sketch", "-o", part, "-k", str(a.kmer),
                        "-s", str(a.sketch_size), "-S", str(a.seed),
                        "-p", str(a.threads)] + batch,
                       check=True, cwd=a.work, stderr=subprocess.DEVNULL)
        parts.append(os.path.join(a.work, part))
        for f in batch:
            os.remove(os.path.join(a.work, f))
        batch.clear()

    for n, q in enumerate(queries, 1):
        contigs = load_contigs(man[q], a.contig_len)
        total = sum(len(c) for c in contigs)
        if not contigs:
            print("  sin contigs utilizables: %s" % q, file=sys.stderr)
            continue
        donor_contigs = {}
        for cls, acc in zip(("near", "far"), donors[q]):
            donor_contigs[cls] = load_contigs(man[acc], a.contig_len) if acc else []
        for c in comps:
            native, kept = take(contigs, c * total, rng)
            for g in conts:
                for cls in ("far", "near"):
                    if g > 0 and not donor_contigs[cls]:
                        continue                      # monospecific genus: no near donor
                    if g == 0 and cls == "near":
                        continue                      # the clean bin is one bin, not two
                    foreign, got = ([], 0) if g == 0 else take(
                        donor_contigs[cls], kept * g / (1 - g), rng)
                    # The bin's identity travels in its file name, because that
                    # is what Mash keeps as the sketch ID: source accession,
                    # completeness, contamination and donor class, in that order.
                    name = "%s~c%.2f~g%.2f~%s" % (q, c, g, cls if g else "none")
                    with open(os.path.join(a.work, name + ".fna"), "w") as fh:
                        write_bin(fh, name, native, foreign)
                    batch.append(name + ".fna")
                    rows.append((name, q, lab[q]["genus"], lab[q]["species_taxid"],
                                 "%.2f" % c, "%.2f" % g, cls if g else "none",
                                 donors[q][0 if cls == "near" else 1] if g else "",
                                 total, kept, got))
                    if len(batch) >= a.batch:
                        flush()
        if n % 25 == 0:
            print("  %d/%d genomas procesados, %d bins" % (n, len(queries), len(rows)))
    flush()

    subprocess.run([a.mash, "paste", a.out] + parts, check=True,
                   stderr=subprocess.DEVNULL)
    for p in parts:
        os.remove(p)

    with open(a.out + "_manifest.tsv", "w") as fh:
        fh.write("bin\tsource\tgenus\tspecies_taxid\tcompleteness\tcontamination\t"
                 "donor_class\tdonor\tsource_bases\tnative_bases\tforeign_bases\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("\n%d bins en %s.msh" % (len(rows), a.out))
    print("manifiesto en %s_manifest.tsv" % a.out)


if __name__ == "__main__":
    main()
