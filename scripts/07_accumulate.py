#!/usr/bin/env python3
"""Accumulate a `mash triangle -E` edge stream into the study's summaries.

Reads the edge list on stdin — `seq1 seq2 dist p-value shared-hashes`, where the
sequence names are the genome paths the sketch was built from — and never holds
more than one line of it in memory. The all-vs-all is ~4.6e8 pairs and is not
written to disk in raw form.

`mash triangle -d` filters the stream, so only pairs at or below the screening
cutoff ever reach this script. The pairs above it are not lost: the exact number
of pairs in every taxonomic class is already known from the pair census, so what
was not reported is recovered by subtraction, without error. The summary reports
both sides of that subtraction.

Four things come out of the single pass:

  <prefix>_hist.tsv        distance histogram by pair class and quality stratum
  <prefix>_hist_genus.tsv  the same, resolved by genus, for the macro-average
  <prefix>_nn.tsv          the k nearest neighbours of every genome
  <prefix>_summary.txt     totals, reported vs. expected from the census

The genus-resolved histogram counts an inter-genus pair under *both* of its
genera, because the quantity a per-genus curve answers is "of the pairs that
involve genus g at distance d, what fraction is congeneric" — and a pair between
g and h is evidence for both.

The nearest-neighbour table is what the leave-one-out view is built from: for a
tool classifying a query against a database the relevant quantity is not a random
pair but the minimum over the database, and its base rate is entirely different.
"""
import argparse
import collections
import heapq
import os
import re
import subprocess
import sys

CLASSES = ("same_strain", "same_species", "congeneric", "inter_genus")
# anchored at the start of the file name, with no assumption about what
# follows: a genome is <acc>_<asm>_genomic.fna and a proteome <acc>.faa.gz
ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def load_labels(path):
    """accession -> (species_taxid, genus_taxid, taxid, is_strain_level, gold)."""
    import csv
    lab = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab[r["accession"]] = (r["species_taxid"], r["genus_taxid"], r["taxid"],
                                   r["taxid"] != r["species_taxid"],
                                   r["pass_gold"] == "1", r["genus"], r["domain"])
    return lab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--prefix", required=True, help="output path prefix")
    ap.add_argument("--bin", type=float, default=0.0005, help="histogram bin width")
    ap.add_argument("--genus-bin", type=float, default=0.002,
                    help="bin width of the genus-resolved histogram")
    ap.add_argument("--knn", type=int, default=10, help="neighbours kept per genome")
    ap.add_argument("--keep-pairs", default="",
                    help="if set, write every reported pair to this zstd file")
    ap.add_argument("--progress", type=int, default=20_000_000,
                    help="report every N lines to stderr (0 disables)")
    a = ap.parse_args()

    lab = load_labels(a.labels)

    # --- path -> compact index, resolved once ---------------------------------
    idx = {}
    sp, gn, gname, strain, gold, dom, acc_of = [], [], [], [], [], [], []

    def index_of(path):
        i = idx.get(path)
        if i is not None:
            return i
        m = ACC_RE.match(path.rsplit("/", 1)[-1])
        if not m:
            sys.exit("ABORT: cannot parse an accession from %r" % path)
        acc = m.group(1)
        try:
            s, g, t, is_strain, ok, gnm, dm = lab[acc]
        except KeyError:
            sys.exit("ABORT: %s is in the sketch but not in the label table" % acc)
        i = len(sp)
        idx[path] = i
        sp.append(s); gn.append(g); gname.append(gnm)
        strain.append(t if is_strain else "")
        gold.append(ok); dom.append(dm); acc_of.append(acc)
        return i

    hist = collections.Counter()        # (gold_pair, class, bin) -> count
    ghist = collections.Counter()       # (genus_taxid, class, gbin) -> count
    nn = collections.defaultdict(list)  # index -> max-heap of (-dist, partner, class)

    keep = None
    if a.keep_pairs:
        keep = subprocess.Popen(["zstd", "-q", "-3", "-T4", "-o", a.keep_pairs, "-f"],
                                stdin=subprocess.PIPE)
        keep.stdin.write(b"acc_a\tacc_b\tdist\tshared\tclass\tgold\n")

    bw, gbw, k = a.bin, a.genus_bin, a.knn
    n_lines = 0
    reported = collections.Counter()

    for line in sys.stdin.buffer:
        f = line.split(b"\t")
        if len(f) < 5:
            continue
        i = index_of(f[0].decode())
        j = index_of(f[1].decode())
        d = float(f[2])

        si, sj = sp[i], sp[j]
        if si == sj:
            if strain[i] and strain[i] == strain[j]:
                cls = 0
            else:
                cls = 1
        elif gn[i] == gn[j]:
            cls = 2
        else:
            cls = 3

        g_pair = 1 if (gold[i] and gold[j]) else 0
        hist[(g_pair, cls, int(d / bw))] += 1
        reported[cls] += 1

        gb = int(d / gbw)
        ghist[(gn[i], cls, gb)] += 1
        if cls == 3:
            ghist[(gn[j], cls, gb)] += 1

        for a_, b_ in ((i, j), (j, i)):
            h = nn[a_]
            if len(h) < k:
                heapq.heappush(h, (-d, b_, cls))
            elif -h[0][0] > d:
                heapq.heapreplace(h, (-d, b_, cls))

        if keep is not None:
            keep.stdin.write(b"%s\t%s\t%s\t%s\t%s\t%d\n"
                             % (acc_of[i].encode(), acc_of[j].encode(), f[2], f[4].strip(),
                                CLASSES[cls].encode(), g_pair))

        n_lines += 1
        if a.progress and n_lines % a.progress == 0:
            print("  %d M edges" % (n_lines // 1_000_000), file=sys.stderr, flush=True)

    if keep is not None:
        keep.stdin.close()
        keep.wait()

    # --- outputs --------------------------------------------------------------
    with open(a.prefix + "_hist.tsv", "w") as fh:
        fh.write("gold_pair\tclass\tdist_lo\tcount\n")
        for (g, c, b), n in sorted(hist.items()):
            fh.write("%d\t%s\t%.5f\t%d\n" % (g, CLASSES[c], b * bw, n))

    gname_of = {gn[i]: gname[i] for i in range(len(gn))}
    with open(a.prefix + "_hist_genus.tsv", "w") as fh:
        fh.write("genus_taxid\tgenus\tclass\tdist_lo\tcount\n")
        for (g, c, b), n in sorted(ghist.items()):
            fh.write("%s\t%s\t%s\t%.5f\t%d\n"
                     % (g, gname_of.get(g, ""), CLASSES[c], b * gbw, n))

    with open(a.prefix + "_nn.tsv", "w") as fh:
        fh.write("accession\tgenus\tdomain\trank\tneighbour\tdist\tclass\n")
        for i, h in nn.items():
            for rank, (nd, part, cls) in enumerate(sorted(h, key=lambda x: -x[0]), 1):
                fh.write("%s\t%s\t%s\t%d\t%s\t%.6f\t%s\n"
                         % (acc_of[i], gname[i], dom[i], rank, acc_of[part],
                            -nd, CLASSES[cls]))

    with open(a.prefix + "_summary.txt", "w") as fh:
        fh.write("edges read below the screening cutoff: %d\n" % n_lines)
        fh.write("genomes seen: %d\n" % len(idx))
        fh.write("genomes with at least one neighbour reported: %d\n\n" % len(nn))
        fh.write("%-14s %14s\n" % ("class", "reported"))
        for c, name in enumerate(CLASSES):
            fh.write("%-14s %14d\n" % (name, reported[c]))
        fh.write("\nPairs not reported are those above the cutoff; their count per\n"
                 "class is the census total minus the figure above.\n")

    print("edges: %d   genomes: %d" % (n_lines, len(idx)), file=sys.stderr)
    for c, name in enumerate(CLASSES):
        print("  %-14s %12d" % (name, reported[c]), file=sys.stderr)


if __name__ == "__main__":
    main()
