#!/usr/bin/env python3
"""Validate the protein-sketch axis against a real AAI, computed with diamond.

Mash over proteomes is a proxy for amino acid identity, not AAI itself, and D7
only buys a third axis if that proxy is shown to track the index the taxonomy
actually uses. This step computes AAI the classical way — reciprocal best hits
between two proteomes, mean identity over those hits (Konstantinidis and Tiedje,
2005) — on a subset of genome pairs, and puts it side by side with the protein
sketch distance, the DNA sketch distance and the alignment-based ANI of §3.8.

The subset is not a random sample of pairs. A random sample would be
inter-genus 239 times out of 240 (§3.1) and would say nothing about the band in
dispute, so genera with several sequenced species are drawn first and unrelated
genomes are added to cover the far range. What the subset is for is the shape of
the relation between two axes, not a frequency, so enrichment is the right
design; the selection is written out so the reader can see what it was.

All the proteomes go into one diamond database and one search, rather than a run
per pair: the same alignment serves every pair in the subset, and the best hit of
a protein in each target genome is read off the shared output.

The question it answers is concrete: where does the 65 % AAI genus criterion
(Konstantinidis and Tiedje, 2005; Barco et al., 2020) fall in protein sketch
distance, and does it agree with the boundary the leave-one-out view marks.
"""
import argparse
import collections
import csv
import gzip
import os
import random
import re
import statistics
import subprocess
import sys

ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def load_labels(path):
    lab = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab[r["accession"]] = r
    return lab


def choose(lab, have, n_genera, per_genus, n_random, seed):
    """Genomes from multi-species genera, plus unrelated ones for the far range."""
    rng = random.Random(seed)
    by_genus = collections.defaultdict(list)
    for acc in have:
        r = lab[acc]
        if r["domain"] == "Bacteria" and r["pass_gold"] == "1":
            by_genus[r["genus"]].append(acc)
    # one genome per species, so a genus does not enter through strain replicates
    rich = []
    for g, accs in by_genus.items():
        by_sp = {}
        for a in accs:
            by_sp.setdefault(lab[a]["species_taxid"], a)
        if len(by_sp) >= per_genus:
            rich.append((g, sorted(by_sp.values())))
    rng.shuffle(rich)
    picked, chosen_genera = [], []
    for g, accs in rich[:n_genera]:
        picked.extend(rng.sample(accs, per_genus))
        chosen_genera.append(g)
    pool = [a for a in have if a not in set(picked)]
    picked.extend(rng.sample(pool, min(n_random, len(pool))))
    return picked, chosen_genera


def write_db(accs, faa_dir, out):
    """One FASTA with every protein tagged by the genome it came from."""
    n = 0
    with open(out, "w") as fh:
        for acc in accs:
            p = os.path.join(faa_dir, acc + ".faa.gz")
            with gzip.open(p, "rt") as g:
                i = 0
                for line in g:
                    if line.startswith(">"):
                        i += 1
                        fh.write(">%s|%d\n" % (acc, i))
                    else:
                        fh.write(line)
                n += i
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--faa", required=True, help="directory of <accession>.faa.gz")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--diamond", default="diamond")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--genera", type=int, default=12, help="multi-species genera to draw")
    ap.add_argument("--per-genus", type=int, default=4, help="genomes per drawn genus")
    ap.add_argument("--random", type=int, default=12, help="unrelated genomes added")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-id", type=float, default=30.0,
                    help="identity floor of a hit, in percent, as in the classical "
                         "definition")
    ap.add_argument("--min-cov", type=float, default=70.0,
                    help="coverage floor of the shorter sequence, in percent")
    ap.add_argument("--mash-pairs", default="",
                    help="<prefix>_pairs.tsv.zst of the protein axis. When given, the "
                         "sketch distance of every pair is joined onto its AAI, which "
                         "is what turns the validation into a correspondence between "
                         "the two scales. Pairs beyond the screening cutoff have no "
                         "distance to join and are written with an empty field rather "
                         "than dropped.")
    ap.add_argument("--keep", action="store_true", help="keep the diamond output")
    a = ap.parse_args()

    lab = load_labels(a.labels)
    have = sorted(acc for acc in lab
                  if os.path.exists(os.path.join(a.faa, acc + ".faa.gz")))
    print("proteomes available: %d" % len(have))

    accs, genera = choose(lab, have, a.genera, a.per_genus, a.random, a.seed)
    print("subset: %d genomes from %d multi-species genera plus %d unrelated"
          % (len(accs), len(genera), a.random))
    print("genera drawn: %s" % ", ".join(sorted(genera)))

    work = a.out + "_work"
    os.makedirs(work, exist_ok=True)
    fasta = os.path.join(work, "subset.faa")
    n_prot = write_db(accs, a.faa, fasta)
    print("proteins in the subset: %d" % n_prot)

    dmnd = os.path.join(work, "subset.dmnd")
    hits = os.path.join(work, "hits.tsv")
    if not os.path.exists(hits):
        subprocess.run([a.diamond, "makedb", "--in", fasta, "-d", dmnd,
                        "--threads", str(a.threads), "--quiet"], check=True)
        subprocess.run([a.diamond, "blastp", "-q", fasta, "-d", dmnd, "-o", hits,
                        "--threads", str(a.threads), "--quiet",
                        "--max-target-seqs", "5000", "--id", str(a.min_id),
                        "--query-cover", str(a.min_cov), "--subject-cover",
                        str(a.min_cov),
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore"],
                       check=True)
    print("diamond output: %s (%s)"
          % (hits, subprocess.run(["du", "-h", hits], capture_output=True,
                                  text=True).stdout.split()[0]))

    # best hit of every protein in every target genome
    best = {}                       # (query protein, target genome) -> (bits, id, sseq)
    with open(hits) as fh:
        for line in fh:
            q, sname, pid, bits = line.rstrip("\n").split("\t")
            qa = q.split("|")[0]
            sa = sname.split("|")[0]
            if qa == sa:
                continue
            k = (q, sa)
            b = float(bits)
            if k not in best or b > best[k][0]:
                best[k] = (b, float(pid), sname)

    # reciprocal best hits, and the AAI they define
    pair_ids = collections.defaultdict(list)
    for (q, sa), (b, pid, sname) in best.items():
        qa = q.split("|")[0]
        back = best.get((sname, qa))
        if back and back[2] == q:
            key = tuple(sorted((qa, sa)))
            if qa < sa:                       # count each reciprocal pair once
                pair_ids[key].append(pid)

    rows = []
    for (x, y), ids in pair_ids.items():
        rows.append((x, y, len(ids), sum(ids) / len(ids)))
    rows.sort(key=lambda r: -r[3])

    def cls(x, y):
        rx, ry = lab[x], lab[y]
        if rx["species_taxid"] == ry["species_taxid"]:
            return "same_species"
        return "congeneric" if rx["genus_taxid"] == ry["genus_taxid"] else "inter_genus"

    dist = {}
    if a.mash_pairs:
        want = {tuple(sorted((x, y))) for x, y, _, _ in rows}
        proc = subprocess.Popen(["zstd", "-dc", a.mash_pairs], stdout=subprocess.PIPE)
        next(proc.stdout)
        for line in proc.stdout:
            f = line.decode().split("\t")
            k = tuple(sorted((f[0], f[1])))
            if k in want:
                dist[k] = float(f[2])
        proc.wait()
        print("sketch distances joined: %d of %d pairs (the rest are beyond the "
              "screening cutoff)" % (len(dist), len(rows)))

    with open(a.out + ".tsv", "w") as fh:
        fh.write("acc_a\tacc_b\tclass\tgenus_a\tgenus_b\tn_bbh\taai\tprot_dist\n")
        for x, y, n, aai in rows:
            d = dist.get(tuple(sorted((x, y))))
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%.4f\t%s\n"
                     % (x, y, cls(x, y), lab[x]["genus"], lab[y]["genus"], n, aai,
                        "%.6f" % d if d is not None else ""))
    with open(a.out + "_subset.txt", "w") as fh:
        fh.write("\n".join(accs) + "\n")

    n_pairs = len(accs) * (len(accs) - 1) // 2
    print("\npairs with a reciprocal-best-hit AAI: %d of %d possible"
          % (len(rows), n_pairs))
    by = collections.Counter(cls(x, y) for x, y, _, _ in rows)
    print("  by class: %s" % dict(by))
    for c in ("same_species", "congeneric", "inter_genus"):
        v = sorted(r[3] for r in rows if cls(r[0], r[1]) == c)
        if v:
            print("  %-14s n=%4d  AAI  min %.1f  mediana %.1f  max %.1f"
                  % (c, len(v), v[0], statistics.median(v), v[-1]))
    print("\nwrote %s.tsv" % a.out)
    if not a.keep:
        print("(the diamond working directory is kept at %s; --keep is the default "
              "for reruns)" % work)


if __name__ == "__main__":
    main()
