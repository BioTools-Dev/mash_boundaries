#!/usr/bin/env python3
"""The second formal genus criterion: percentage of conserved proteins, measured.

Taxonomy delimits the genus by two genome-based criteria, not one. §3.9 places
the first, amino acid identity at ~65 %, on the axes of this study. This step
places the second, the percentage of conserved proteins at 50 % (Qin et al.,
2014), and asks the question the AAI validation could not: the AAI criterion is
an identity and therefore lives naturally on an identity axis, while POCP counts
shared genes and is independent of how similar they are. If both criteria land
in the same region of Mash distance, the genus cutoff proposed by this work is
supported by two indices that measure different things.

POCP is defined pairwise as the sum of the conserved proteins of the two genomes
over the sum of their protein counts, a protein counting as conserved when it has
a match in the other genome with an E-value below 1e-5, identity above 40 % and an
alignable region above 50 % of its own length. Those three thresholds are the
definition and are not tuned here.

Two departures from the original are declared rather than hidden. The alignment
is diamond and not BLASTP, which the reference implementation of POCP reports as
equivalent in ultra-sensitive mode to within ~0.16 percentage points (Hölzer,
2024); the sensitivity actually used is written into the output. And, as in
§3.9, all the proteomes go into one database and one search rather than a run
per pair, so a protein's matches in every other genome are read off the same
alignment.

The subset is the one the AAI validation drew, taken from its output rather than
redrawn, so that every pair carries both criteria and the two are comparable
genome by genome.
"""
import argparse
import collections
import csv
import gzip
import math
import os
import random
import statistics
import subprocess
import sys


def load_labels(path):
    lab = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab[r["accession"]] = r
    return lab


def write_db(accs, faa_dir, out):
    """One FASTA with every protein tagged by the genome it came from.

    Returns the protein count of each genome, which is the denominator of POCP
    and has to be the number of proteins actually searched, not the number the
    assembly report claims.
    """
    total = {}
    with open(out, "w") as fh:
        for i, acc in enumerate(accs):
            n = 0
            with gzip.open(os.path.join(faa_dir, acc + ".faa.gz"), "rt") as g:
                for line in g:
                    if line.startswith(">"):
                        n += 1
                        fh.write(">%d|%d\n" % (i, n))
                    else:
                        fh.write(line)
            total[acc] = n
    return total


def run_diamond(dmnd, fasta, work, threads, sensitivity, min_id, min_cov, evalue,
                max_target):
    """One all-versus-all search under the POCP thresholds.

    Only the two identifiers are asked for: the three thresholds of the
    definition are applied by diamond itself, so a line in the output *is* a
    conserved-protein event and nothing has to be re-filtered afterwards.
    """
    db = os.path.join(work, "pocp.dmnd")
    hits = os.path.join(work, "hits.tsv")
    if os.path.exists(hits):
        print("reusing %s" % hits)
        return hits
    subprocess.run([dmnd, "makedb", "--in", fasta, "-d", db,
                    "--threads", str(threads), "--quiet"], check=True)
    cmd = [dmnd, "blastp", "-q", fasta, "-d", db, "-o", hits + ".part",
           "--threads", str(threads), "--quiet",
           "--max-target-seqs", str(max_target),
           "--evalue", str(evalue),
           "--id", str(min_id), "--query-cover", str(min_cov),
           "--outfmt", "6", "qseqid", "sseqid"]
    if sensitivity:
        cmd.append("--" + sensitivity)
    print("running: %s" % " ".join(cmd[:2] + cmd[6:]), file=sys.stderr)
    subprocess.run(cmd, check=True)
    os.rename(hits + ".part", hits)     # a partial file must never look finished
    return hits


def conserved(hits, n_genomes, cap):
    """For every ordered pair of genomes, how many proteins of the first have a
    match in the second.

    A protein counts once per target genome however many matches it has there,
    so the per-query target set is collected as a bitmask and only then folded
    into the counts. That also makes the pass independent of whether diamond
    groups its output by query.
    """
    seen = collections.defaultdict(int)          # query protein -> bitmask
    per_query = collections.Counter()            # to detect a truncated search
    for line in open(hits):
        q, s = line.rstrip("\n").split("\t")
        per_query[q] += 1
        qi = q.split("|", 1)[0]
        si = int(s.split("|", 1)[0])
        if qi == s.split("|", 1)[0]:
            continue                             # a genome against itself
        seen[q] |= 1 << si
    # If any query fills the --max-target-seqs quota its hit list was cut, and a
    # genome missing from it would silently lower POCP. The cap is generous for
    # this subset size but the check costs nothing and the failure would be
    # invisible otherwise.
    if per_query and max(per_query.values()) >= cap:
        n = sum(1 for v in per_query.values() if v >= cap)
        print("WARNING: %d query proteins reached the --max-target-seqs cap of %d; "
              "raise it or POCP is underestimated" % (n, cap), file=sys.stderr)
    c = [[0] * n_genomes for _ in range(n_genomes)]
    for q, mask in seen.items():
        i = int(q.split("|", 1)[0])
        while mask:
            low = mask & -mask
            c[i][low.bit_length() - 1] += 1
            mask ^= low
    return c


def mash_distances(accs, manifest, sketch_dir, tag, mash, k, s, threads):
    """Every pair of the subset, with no screening cutoff.

    The all-versus-all of §3.6 only reports what falls inside its window, and
    most inter-genus pairs of this subset fall outside it. Re-sketching the 60
    genomes with the same parameters gives the same sketches and therefore the
    same distances, and lets every pair carry a value instead of a censored one.
    """
    path = {}
    with open(manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            path[r["accession"]] = r["path"]
    lst = os.path.join(sketch_dir, "%s_subset.txt" % tag)
    with open(lst, "w") as fh:
        fh.write("\n".join(path[a] for a in accs) + "\n")
    msh = os.path.join(sketch_dir, "%s_subset_s%d" % (tag, s))
    if not os.path.exists(msh + ".msh"):
        subprocess.run([mash, "sketch", "-p", str(threads), "-k", str(k),
                        "-s", str(s), "-o", msh, "-l", lst],
                       check=True, capture_output=True)
    out = subprocess.run([mash, "dist", "-p", str(threads), msh + ".msh",
                          msh + ".msh"], check=True, capture_output=True, text=True)
    base = {os.path.basename(path[a]): a for a in accs}
    d, censored = {}, set()
    for line in out.stdout.splitlines():          # mash reports each pair twice
        f = line.split("\t")
        a, b = base[os.path.basename(f[0])], base[os.path.basename(f[1])]
        if a == b:
            continue
        # A pair with no shared hash gets distance 1 from mash. That is the
        # estimator saying it cannot see the pair, not a measurement of it, and
        # writing it as a number would put a fabricated value in the table.
        if f[4].split("/")[0] == "0":
            censored.add(tuple(sorted((a, b))))
            continue
        d[tuple(sorted((a, b)))] = float(f[2])
    return d, len(censored)


def join(path, keys, col):
    """A column of an existing result table, keyed by the accession pair."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = tuple(sorted((r["acc_a"], r["acc_b"])))
            if k in keys and r.get(col):
                out[k] = float(r[col])
    return out


def stratified(pairs_file, lab, faa_dir, bands, per_band, cap, seed):
    """Genomes chosen so that congeneric pairs land inside the genus window.

    The AAI subset was drawn to span the range of AAI, and inside the proposed
    genus window it is nearly empty: a genus contributes one genome per species,
    and two species of the same genus are rarely that close. Testing whether the
    window ever admits a pair the criterion would reject needs pairs *in* the
    window, so they are sampled directly from the all-versus-all, stratified by
    distance band.

    One pair per genus per band, because the ten largest genera hold most of the
    congeneric signal and letting them fill the sample would make the answer a
    statement about *Streptomyces*.
    """
    rng = random.Random(seed)
    cand = collections.defaultdict(list)
    proc = subprocess.Popen(["zstd", "-dc", pairs_file], stdout=subprocess.PIPE)
    next(proc.stdout)
    for line in proc.stdout:
        f = line.decode().split("\t")
        if f[4] != "congeneric" or f[5].strip() != "1":
            continue                      # gold-standard genomes only, as in §2.1
        d = float(f[2])
        for i, (lo, hi) in enumerate(bands):
            if lo <= d < hi:
                cand[i].append((f[0], f[1], d))
                break
    proc.kill()

    picked, seen_genus = [], set()
    for i, (lo, hi) in enumerate(bands):
        rng.shuffle(cand[i])
        taken = 0
        for x, y, d in cand[i]:
            if taken >= per_band or len(set(picked)) >= cap:
                break
            g = lab[x]["genus"]
            if (i, g) in seen_genus:
                continue
            if not all(os.path.exists(os.path.join(faa_dir, a + ".faa.gz"))
                       for a in (x, y)):
                continue
            seen_genus.add((i, g))
            picked.extend([x, y])
            taken += 1
        print("  band %.2f - %.2f: %d candidate pairs, %d drawn"
              % (lo, hi, len(cand[i]), taken))
    return sorted(set(picked))


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def crossing(rows, level, lo, hi):
    """The Mash distance at which POCP crosses `level`, read off the data.

    Two figures are reported and they answer different questions: the largest
    distance at which every pair is still above the criterion, and the smallest
    at which none is. Between them the criterion and the distance disagree, and
    the width of that interval is the honest form of the answer — a single
    crossing point would hide it.
    """
    above = [d for d, p in rows if p >= level]
    below = [d for d, p in rows if p < level]
    if not above or not below:
        return None
    return (max([d for d in above if all(p >= level for dd, p in rows if dd <= d)]
                or [lo]),
            min([d for d in below if all(p < level for dd, p in rows if dd >= d)]
                or [hi]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--manifest", required=True, help="data/genome_manifest.tsv")
    ap.add_argument("--faa", required=True, help="directory of <accession>.faa.gz")
    ap.add_argument("--subset", default="",
                    help="results/aai_subset.txt — the genomes the AAI validation "
                         "drew, so that both criteria describe the same pairs")
    ap.add_argument("--stratified-from", default="",
                    help="<prefix>_pairs.tsv.zst; instead of reusing the AAI "
                         "subset, draw congeneric pairs stratified by distance "
                         "band so that the genus window is populated")
    ap.add_argument("--bands", default="0.02:0.04,0.04:0.06,0.06:0.08,0.08:0.10,"
                                       "0.10:0.12,0.12:0.15,0.15:0.20",
                    help="distance bands for --stratified-from. They control "
                         "which pairs are drawn, not where the window is read: "
                         "the window analysis filters on distance afterwards, so "
                         "moving a cutoff does not require redrawing the subset")
    ap.add_argument("--per-band", type=int, default=6)
    ap.add_argument("--max-genomes", type=int, default=80)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--aai", default="", help="results/aai.tsv, to carry the AAI "
                                              "and the protein distance across")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--diamond", default="diamond")
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--sketch-dir", default=".")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--sensitivity", default="ultra-sensitive",
                    choices=["", "sensitive", "very-sensitive", "ultra-sensitive"],
                    help="diamond sensitivity mode; recorded in the summary "
                         "because POCP is a count of detected homologues and so "
                         "depends on it")
    ap.add_argument("--min-id", type=float, default=40.0,
                    help="identity floor, per the definition of POCP")
    ap.add_argument("--min-cov", type=float, default=50.0,
                    help="alignable fraction of the query, per the definition")
    ap.add_argument("--evalue", default="1e-5", help="per the definition")
    ap.add_argument("--max-target-seqs", type=int, default=5000,
                    help="hits kept per query protein. It has to exceed the number "
                         "of homologues one protein can have across the whole "
                         "subset, or a genome drops out of a query's list and its "
                         "POCP falls silently; the run warns if any query reaches it")
    ap.add_argument("--level", type=float, default=50.0,
                    help="the published genus criterion, in percent")
    ap.add_argument("--mash-k", type=int, default=21)
    ap.add_argument("--mash-s", type=int, default=10000)
    a = ap.parse_args()

    ceiling = -math.log(2.0 / (a.mash_s + 1)) / a.mash_k

    lab = load_labels(a.labels)
    if a.stratified_from:
        bands = [tuple(float(v) for v in b.split(":"))
                 for b in a.bands.split(",")]
        print("drawing congeneric pairs by distance band:")
        accs = stratified(a.stratified_from, lab, a.faa, bands, a.per_band,
                          a.max_genomes, a.seed)
        print("subset: %d genomes drawn from %s" % (len(accs), a.stratified_from))
    elif a.subset:
        accs = [l.strip() for l in open(a.subset) if l.strip()]
        print("subset: %d genomes, taken from %s" % (len(accs), a.subset))
    else:
        sys.exit("give either --subset or --stratified-from")
    idx = {acc: i for i, acc in enumerate(accs)}

    work = a.out + "_work"
    os.makedirs(work, exist_ok=True)
    fasta = os.path.join(work, "subset.faa")
    if os.path.exists(fasta + ".counts"):
        total = {}
        with open(fasta + ".counts") as fh:
            for line in fh:
                acc, n = line.split()
                total[acc] = int(n)
    else:
        total = write_db(accs, a.faa, fasta)
        with open(fasta + ".counts", "w") as fh:
            for acc in accs:
                fh.write("%s\t%d\n" % (acc, total[acc]))
    print("proteins in the subset: %d" % sum(total.values()))

    hits = run_diamond(a.diamond, fasta, work, a.threads, a.sensitivity,
                       a.min_id, a.min_cov, a.evalue, a.max_target_seqs)
    c = conserved(hits, len(accs), a.max_target_seqs)

    dist, censored = mash_distances(accs, a.manifest, a.sketch_dir,
                                    os.path.basename(a.out), a.mash,
                                    a.mash_k, a.mash_s, a.threads)
    print("Mash distances for the subset: %d measured, %d with no shared hash "
          "at s = %d" % (len(dist), censored, a.mash_s))

    def cls(x, y):
        rx, ry = lab[x], lab[y]
        if rx["species_taxid"] == ry["species_taxid"]:
            return "same_species"
        return "congeneric" if rx["genus_taxid"] == ry["genus_taxid"] else "inter_genus"

    keys = {tuple(sorted((x, y))) for i, x in enumerate(accs) for y in accs[i + 1:]}
    aai = join(a.aai, keys, "aai")
    pdist = join(a.aai, keys, "prot_dist")

    rows = []
    for i, x in enumerate(accs):
        for y in accs[i + 1:]:
            k = tuple(sorted((x, y)))
            cx, cy = c[idx[x]][idx[y]], c[idx[y]][idx[x]]
            pocp = 100.0 * (cx + cy) / (total[x] + total[y])
            rows.append((x, y, cls(x, y), cx, total[x], cy, total[y], pocp,
                         dist.get(k), aai.get(k), pdist.get(k)))

    with open(a.out + ".tsv", "w") as fh:
        fh.write("acc_a\tacc_b\tclass\tgenus_a\tgenus_b\tconserved_a\tproteins_a\t"
                 "conserved_b\tproteins_b\tpocp\tmash_dist\taai\tprot_dist\n")
        for r in sorted(rows, key=lambda r: -r[7]):
            x, y, cl, cx, tx, cy, ty, pocp, d, ai, pd = r
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%s\t%s\t%s\n"
                     % (x, y, cl, lab[x]["genus"], lab[y]["genus"], cx, tx, cy, ty,
                        pocp, "" if d is None else "%.6f" % d,
                        "" if ai is None else "%.4f" % ai,
                        "" if pd is None else "%.6f" % pd))

    with open(a.out + "_summary.txt", "w") as fh:
        def out(s=""):
            print(s)
            fh.write(s + "\n")

        out("POCP over %s"
            % ("congeneric pairs drawn by distance band from %s"
               % a.stratified_from if a.stratified_from
               else "the subset of %s" % a.subset))
        out("  genomes %d, pairs %d, proteins %d"
            % (len(accs), len(rows), sum(total.values())))
        out("  diamond %s, id >= %g %%, query cover >= %g %%, E < %s"
            % (a.sensitivity or "default", a.min_id, a.min_cov, a.evalue))
        out("  Mash k = %d, s = %d: %d pairs measured, %d with no shared hash "
            "and therefore left without a distance"
            % (a.mash_k, a.mash_s, len(dist), censored))
        out()
        out("%-14s %6s %8s %8s %8s" % ("class", "pairs", "min", "median", "max"))
        for cl in ("same_species", "congeneric", "inter_genus"):
            v = sorted(r[7] for r in rows if r[2] == cl)
            if v:
                out("%-14s %6d %8.2f %8.2f %8.2f"
                    % (cl, len(v), v[0], statistics.median(v), v[-1]))
        cong = [r[7] for r in rows if r[2] in ("congeneric", "same_species")]
        inter = [r[7] for r in rows if r[2] == "inter_genus"]
        if cong and inter:
            out()
            out("separation at the %g %% criterion:" % a.level)
            out("  same-genus pairs below it : %d of %d"
                % (sum(1 for v in cong if v < a.level), len(cong)))
            out("  inter-genus pairs above it: %d of %d"
                % (sum(1 for v in inter if v >= a.level), len(inter)))
            out("  lowest same-genus POCP %.2f, highest inter-genus POCP %.2f"
                % (min(cong), max(inter)))

        # The two formal genus criteria placed on the axis the cutoffs are
        # stated in. They measure different things — POCP counts shared genes,
        # AAI measures how similar they are — so where they land relative to
        # each other is the check, and it is reported for both or for neither.
        for level, key, name, tol in ((a.level, 7, "%g %% POCP" % a.level, 5.0),
                                      (65.0, 9, "65 % AAI", 2.0)):
            have = sorted([(r[8], r[key]) for r in rows
                           if r[8] is not None and r[key] is not None])
            if not have:
                continue
            out()
            out("where the %s criterion falls in Mash distance (%d pairs):"
                % (name, len(have)))
            cr = crossing(have, level, 0.0, 1.0)
            if cr:
                out("  every pair is above it up to d = %.4f" % cr[0])
                out("  no pair is above it from  d = %.4f" % cr[1])
            band = [d for d, p in have if abs(p - level) <= tol]
            if band:
                out("  pairs within %g points of it: %d, Mash distance "
                    "%.4f - %.4f (median %.4f)"
                    % (tol, len(band), min(band), max(band),
                       statistics.median(band)))

        # correlations, over the range where the DNA sketch still measures
        both = [(r[8], r[7], r[9]) for r in rows
                if r[8] is not None and r[8] < ceiling and r[9] is not None]
        if len(both) > 2:
            out()
            out("over the %d pairs below the measurable ceiling (%.4f at s = %d):"
                % (len(both), ceiling, a.mash_s))
            out("  POCP against Mash distance   r = %+.3f"
                % pearson([x[0] for x in both], [x[1] for x in both]))
            out("  POCP against AAI             r = %+.3f"
                % pearson([x[2] for x in both], [x[1] for x in both]))
        out()
        out("POCP by Mash distance band:")
        out("%-16s %6s %8s %8s %8s" % ("band", "pairs", "min", "median", "max"))
        have = sorted([(r[8], r[7]) for r in rows if r[8] is not None])
        edges = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.15, 0.20, 0.25, 0.30, 1.01]
        for lo, hi in zip(edges, edges[1:]):
            v = sorted(p for d, p in have if lo <= d < hi)
            if v:
                out("%-16s %6d %8.2f %8.2f %8.2f"
                    % ("%.2f - %.2f" % (lo, hi), len(v), v[0],
                       statistics.median(v), v[-1]))

    with open(a.out + "_subset.txt", "w") as fh:
        fh.write("\n".join(accs) + "\n")
    print("\nwrote %s.tsv, %s_summary.txt and %s_subset.txt" % (a.out, a.out, a.out))


if __name__ == "__main__":
    main()
