#!/usr/bin/env python3
"""The ANI axis: what alignment-based identity says about the same pairs.

Mash compared only against itself measures nothing. This script reads the skani
all-vs-all against the Mash all-vs-all and answers three separate questions:

  coverage      of the pairs Mash reports at a given distance, how many does
                skani resolve at all — and are there pairs with an ANI that Mash
                never reported? The second half is the check that neither metric
                decides what the other gets to see (D6).

  correspondence
                where the two axes place each other: the ANI of the pairs Mash
                puts at each distance, and the distance of the pairs skani puts
                at each identity. The conversion d = 1 - ANI that Mash's own
                paper proposes is a claim about closely related genomes; how far
                it holds into the genus band is measured here, not assumed.

  the boundary  the leave-one-out view of the genus cutoff rebuilt on the ANI
                axis, so the question stops being "does Mash approximate ANI"
                and becomes "where do two independent indices put the same
                boundary, and does either resolve the genus".

Alignment-based ANI stops being defined near 80 % identity, so a query whose
nearest usable relative is beyond that has no ANI at all. Those queries are
counted and reported as such: the silence of the axis is a property of the axis
and is part of the answer, not a gap to be filled.
"""
import argparse
import collections
import csv
import re
import statistics
import subprocess
import sys

CLASSES = ("same_strain", "same_species", "congeneric", "inter_genus")
SAME_GENUS = CLASSES[:3]
ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")
SCENARIOS = {
    "as_reported": lambda c: True,
    "no_self_strain": lambda c: c != "same_strain",
    "novel_species": lambda c: c not in ("same_strain", "same_species"),
}


def load_labels(path):
    """accession -> (species_taxid, genus_taxid, strain_taxid or '', genus, domain)."""
    lab = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            strain = r["taxid"] if r["taxid"] != r["species_taxid"] else ""
            lab[r["accession"]] = (r["species_taxid"], r["genus_taxid"], strain,
                                   r["genus"], r["domain"])
    return lab


def pair_class(a, b, lab):
    sa, ga, ta, _, _ = lab[a]
    sb, gb, tb, _, _ = lab[b]
    if sa == sb:
        return "same_strain" if (ta and ta == tb) else "same_species"
    return "congeneric" if ga == gb else "inter_genus"


def zopen(path):
    """Line iterator over a plain or zstd-compressed file, in bytes."""
    if path.endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)
        return p.stdout, p
    return open(path, "rb"), None


def read_skani(path, lab):
    """Pair table and per-genome neighbour lists from `skani triangle -E` output.

    The key is the ordered accession pair as bytes, in both directions, so the
    Mash stream can be probed with a single dictionary lookup per line without
    building anything per line.
    """
    fh, proc = zopen(path)
    pairs, nbr = {}, collections.defaultdict(list)
    unknown = 0
    next(fh)                                    # header
    for line in fh:
        f = line.split(b"\t")
        ma = ACC_RE.search(f[0].decode())
        mb = ACC_RE.search(f[1].decode())
        if not (ma and mb):
            sys.exit("ABORT: cannot parse accessions from %r" % line[:120])
        a, b = ma.group(1), mb.group(1)
        if a not in lab or b not in lab:
            unknown += 1
            continue
        ani = float(f[2])
        afmin = min(float(f[3]), float(f[4]))
        cls = pair_class(a, b, lab)
        ka, kb = a.encode(), b.encode()
        pairs[ka + b"|" + kb] = (ani, afmin)
        pairs[kb + b"|" + ka] = (ani, afmin)
        nbr[a].append((ani, b, cls, afmin))
        nbr[b].append((ani, a, cls, afmin))
    if proc:
        proc.wait()
    else:
        fh.close()
    if unknown:
        print("skani pairs skipped for want of a label: %d" % unknown, file=sys.stderr)
    return pairs, nbr


def scan_mash(path, pairs, bw, progress):
    """Join the Mash edge stream against the ANI table in one pass.

    Returns per (distance bin, class): pairs seen, pairs with an ANI, pairs with
    an ANI at the stricter alignment fraction, and the ANI values themselves.
    """
    fh, proc = zopen(path)
    seen = collections.Counter()          # (bin, class) -> mash pairs
    hit = collections.Counter()           # (bin, class) -> with ANI
    hit15 = collections.Counter()         # (bin, class) -> with ANI, af >= 15
    anis = collections.defaultdict(list)  # (bin, class) -> ANI values
    dist_of = {}                          # ani bin -> mash distances
    matched = set()
    n = 0
    next(fh)                              # header
    for line in fh:
        f = line.split(b"\t")
        d = float(f[2])
        cls = f[4].decode()
        b = int(d / bw)
        seen[(b, cls)] += 1
        seen[(b, "all")] += 1          # the conversion between axes is over every pair
        key = f[0] + b"|" + f[1]
        v = pairs.get(key)
        if v is not None:
            ani, afmin = v
            hit[(b, cls)] += 1
            hit[(b, "all")] += 1
            if afmin >= 15.0:
                hit15[(b, cls)] += 1
                hit15[(b, "all")] += 1
            anis[(b, cls)].append(ani)
            anis[(b, "all")].append(ani)
            dist_of.setdefault(int(ani * 2) / 2.0, []).append(d)
            matched.add(key)
        n += 1
        if progress and n % progress == 0:
            print("  %d M edges" % (n // 1_000_000), file=sys.stderr, flush=True)
    if proc:
        proc.wait()
    else:
        fh.close()
    return seen, hit, hit15, anis, dist_of, matched, n


def quart(v):
    v = sorted(v)
    if not v:
        return (float("nan"),) * 3
    n = len(v)
    return (v[n // 4], statistics.median(v), v[(3 * n) // 4])


def conversion_fit(anis, seen, hit, bw, upto):
    """Fit the corrected distance-to-identity conversion on the band medians.

    The customary conversion is ANI = 1 - d. Read over these pairs it is
    optimistic, and the excess grows with distance rather than sitting at a
    constant offset, so the correction is a slope and not a shift: ANI = 1 - k d.
    The fit is a least squares through the origin, which is not a convenience —
    identity is 100 % at distance zero by construction, and an intercept would
    only absorb the curvature of the far bands into a quantity with no meaning.

    Only bands where the ANI axis resolves most of what Mash reports enter the
    fit; beyond that the median is drawn from a surviving minority and describes
    the pairs skani still aligns rather than the pairs that exist.

    Returns (k, bands used, and a function giving the largest deviation in
    points of any coefficient from the band medians, so that the rounded form
    and the coefficient reported elsewhere can both be checked here).
    """
    pts = []
    for (b, cls), vals in anis.items():
        if cls != "all" or b == 0:
            continue
        if hit[(b, cls)] / seen[(b, cls)] < 0.5:
            continue
        d = (b + 0.5) * bw
        if d <= upto:
            pts.append((d, statistics.median(vals) / 100.0))
    k = sum(d * (1 - a) for d, a in pts) / sum(d * d for d, _ in pts)

    def deviation(coef):
        return 100.0 * max(abs((1 - coef * d) - a) for d, a in pts)

    return k, len(pts), deviation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skani", required=True, help="skani triangle -E output (.tsv or .zst)")
    ap.add_argument("--mash-pairs", required=True, help="<prefix>_pairs.tsv.zst from the accumulator")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--relabel", default="",
                    help="TSV with accession/species_key/genus_key; when given the "
                         "leave-one-out view is recomputed under those labels on "
                         "the same ANI values, as in the GTDB comparison")
    ap.add_argument("--mash-calls", default="",
                    help="loo_*_calls.tsv, to put the two axes side by side on the "
                         "same queries")
    ap.add_argument("--out", required=True, help="output prefix, e.g. results/ani")
    ap.add_argument("--bin", type=float, default=0.005, help="Mash distance bin")
    ap.add_argument("--af-floors", default="5,15,50",
                    help="floors on the smaller of the two alignment fractions "
                         "for the leave-one-out view; the first is the one the "
                         "call table is written for")
    ap.add_argument("--conversion-upto", type=float, default=0.13,
                    help="upper distance for the corrected conversion fit; the "
                         "default is the genus cutoff, the far end of the range "
                         "the conversion is meant to be used over")
    ap.add_argument("--conversion-k", type=float, default=1.12,
                    help="the rounded coefficient the paper states; every run "
                         "reports how far it sits from its own band medians, so "
                         "a single stated form can be checked at every sketch size")
    ap.add_argument("--progress", type=int, default=20_000_000)
    a = ap.parse_args()

    lab = load_labels(a.labels)
    print("labels: %d genomes" % len(lab), file=sys.stderr)

    pairs, nbr = read_skani(a.skani, lab)
    n_ani = len(pairs) // 2
    print("skani: %d pairs, %d genomes with at least one" % (n_ani, len(nbr)),
          file=sys.stderr)

    seen, hit, hit15, anis, dist_of, matched, n_mash = scan_mash(
        a.mash_pairs, pairs, a.bin, a.progress)
    orphan = n_ani - len(matched)          # the Mash stream reports each pair once
    print("mash: %d edges; ANI pairs the Mash window does not contain: %d"
          % (n_mash, orphan), file=sys.stderr)

    # --- coverage and correspondence, by Mash distance ------------------------
    with open(a.out + "_coverage.tsv", "w") as fh:
        fh.write("dist_lo\tdist_hi\tclass\tmash_pairs\twith_ani\tcoverage\t"
                 "with_ani_af15\tcoverage_af15\tani_q1\tani_median\tani_q3\n")
        for (b, cls) in sorted(seen, key=lambda k: (k[0], k[1])):
            q1, md, q3 = quart(anis.get((b, cls), []))
            fh.write("%.4f\t%.4f\t%s\t%d\t%d\t%.6f\t%d\t%.6f\t%.2f\t%.2f\t%.2f\n"
                     % (b * a.bin, (b + 1) * a.bin, cls, seen[(b, cls)],
                        hit[(b, cls)], hit[(b, cls)] / seen[(b, cls)],
                        hit15[(b, cls)], hit15[(b, cls)] / seen[(b, cls)],
                        q1, md, q3))

    # the same correspondence read the other way round
    with open(a.out + "_dist_by_ani.tsv", "w") as fh:
        fh.write("ani_lo\tani_hi\tpairs\tdist_q1\tdist_median\tdist_q3\timplied_dist\n")
        for lo in sorted(dist_of):
            q1, md, q3 = quart(dist_of[lo])
            fh.write("%.1f\t%.1f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n"
                     % (lo, lo + 0.5, len(dist_of[lo]), q1, md, q3, 1.0 - lo / 100.0))

    # --- the correction to the customary conversion ----------------------------
    k, nb, deviation = conversion_fit(anis, seen, hit, a.bin, a.conversion_upto)
    ref = a.conversion_k                     # the form the paper states, checked here too
    with open(a.out + "_conversion.txt", "w") as fh:
        fh.write("corrected conversion, fitted on the band medians of every pair\n")
        fh.write("range      0 < d <= %.3f, %d bands of width %.3f\n"
                 % (a.conversion_upto, nb, a.bin))
        fh.write("fit        ANI = 1 - %.4f d   (rounded: 1 - %.2f d, off by at "
                 "most %.2f points)\n" % (k, round(k, 2), deviation(round(k, 2))))
        fh.write("inverse    d = %.4f (1 - ANI)\n" % (1.0 / k))
        fh.write("stated     ANI = 1 - %.2f d, off by at most %.2f points\n"
                 % (ref, deviation(ref)))
        for d in (0.02, 0.043, 0.05, 0.08, 0.10, 0.13):
            fh.write("  d = %.3f   ANI = %.2f %%   (1 - d would give %.2f %%)\n"
                     % (d, 100 * (1 - ref * d), 100 * (1 - d)))
        for ani in (95.0, 92.0, 90.0, 86.0):
            fh.write("  ANI = %.1f %%   d = %.4f\n"
                     % (ani, (1 - ani / 100.0) / ref))
    print("conversion: ANI = 1 - %.4f d over %d bands; the stated 1 - %.2f d is "
          "off by at most %.2f points" % (k, nb, ref, deviation(ref)),
          file=sys.stderr)

    # --- the boundary on the ANI axis -----------------------------------------
    relabel = {}
    if a.relabel:
        with open(a.relabel) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["species_key"] and r["genus_key"]:
                    relabel[r["accession"]] = (r["species_key"], r["genus_key"])

    def reclass(q, n, original):
        if original == "same_strain":
            return original                      # strain identity is not taxonomic
        qs, qg = relabel[q]
        ns, ng = relabel[n]
        if qs == ns:
            return "same_species"
        return "congeneric" if qg == ng else "inter_genus"

    # A high ANI over a small aligned fraction is not evidence of relatedness at
    # the rank in question: two genomes of different genera can share a plasmid
    # or a ribosomal region at 99 % identity over 5 % of their length. The
    # boundary is therefore rebuilt under several floors on the alignment
    # fraction, and the floor turns out to matter more than the identity cutoff.
    floors = [float(x) for x in a.af_floors.split(",")]
    primary = floors[0]
    curve = collections.defaultdict(lambda: [0, 0])   # (floor, scenario, bin) -> [same genus, n]
    silent = collections.Counter()
    best = collections.defaultdict(dict)              # (floor, scenario) -> acc -> (ani, partner, cls)
    with open(a.out + "_loo_calls.tsv", "w") as fh:
        fh.write("scenario\taccession\tgenus\tdomain\tneighbour\tani\taf_min\tclass\tsame_genus\n")
        for acc in lab:
            lst = nbr.get(acc, [])
            if relabel:
                lst = [(ani, n, reclass(acc, n, c), af) for ani, n, c, af in lst
                       if acc in relabel and n in relabel]
            lst.sort(key=lambda x: -x[0])
            genus = relabel[acc][1] if (relabel and acc in relabel) else lab[acc][3]
            for floor in floors:
                for name, keep in SCENARIOS.items():
                    top = next((t for t in lst if keep(t[2]) and t[3] >= floor), None)
                    if top is None:
                        silent[(floor, name)] += 1
                        continue
                    ani, part, cls, afmin = top
                    same = cls in SAME_GENUS
                    b = int(ani * 2)                   # 0.5-point bins
                    curve[(floor, name, b)][1] += 1
                    curve[(floor, name, b)][0] += int(same)
                    best[(floor, name)][acc] = (ani, part, cls)
                    if floor == primary:
                        fh.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%d\n"
                                 % (name, acc, genus, lab[acc][4], part, ani,
                                    afmin, cls, int(same)))

    with open(a.out + "_loo.tsv", "w") as fh:
        fh.write("af_floor\tscenario\tani_lo\tani_hi\tsame_genus\ttotal\tp_same_genus\t"
                 "cum_total\tcum_p\n")
        for floor in floors:
            for name in SCENARIOS:
                bins = sorted((b for (f, s, b) in curve if f == floor and s == name),
                              reverse=True)
                ch = ct = 0
                for b in bins:
                    sg, tot = curve[(floor, name, b)]
                    ch += sg
                    ct += tot
                    fh.write("%.0f\t%s\t%.1f\t%.1f\t%d\t%d\t%.6f\t%d\t%.6f\n"
                             % (floor, name, b / 2.0, b / 2.0 + 0.5, sg, tot,
                                sg / tot, ct, ch / ct))

    # --- console and file summary ---------------------------------------------
    out = []
    p = out.append
    p("skani pairs with an ANI: %d" % n_ani)
    p("Mash edges in the window: %d" % n_mash)
    p("ANI pairs absent from the Mash window: %d (%.4f %%)"
      % (orphan, 100.0 * orphan / n_ani))
    p("genomes with at least one ANI neighbour: %d of %d" % (len(nbr), len(lab)))
    p("")
    p("Coverage of the Mash bands by the ANI axis (all classes):")
    p("%-14s %14s %14s %9s %9s" % ("Mash band", "mash pairs", "with ANI", "cov %", "median ANI"))
    edges = [0.0, 0.05, 0.08, 0.13, 0.15, 0.20, 0.25, 0.30, 0.40]
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        s = sum(v for (b, c), v in seen.items() if c == "all" and lo <= b * a.bin < hi)
        h = sum(v for (b, c), v in hit.items() if c == "all" and lo <= b * a.bin < hi)
        vals = [x for (b, c), lst in anis.items()
                if c == "all" and lo <= b * a.bin < hi for x in lst]
        p("%-14s %14d %14d %8.2f%% %9s"
          % ("%.2f-%.2f" % (lo, hi), s, h, 100.0 * h / s if s else 0.0,
             "%.2f" % statistics.median(vals) if vals else "-"))
    p("")
    p("Leave-one-out on the ANI axis, scenario novel_species, by floor on the")
    p("alignment fraction:")
    p("%-10s %s" % ("", "  ".join("%22s" % ("af >= %g" % f) for f in floors)))
    p("%-10s %s" % ("ANI cutoff", "  ".join("%10s %11s" % ("calls", "precision")
                                            for _ in floors)))
    for t in (98.0, 96.0, 95.0, 94.0, 93.0, 92.0, 90.0, 88.0, 86.0, 85.0, 83.0, 80.0):
        cells = []
        for floor in floors:
            tot = sum(v[1] for (f, s, b), v in curve.items()
                      if f == floor and s == "novel_species" and b / 2.0 >= t)
            sg = sum(v[0] for (f, s, b), v in curve.items()
                     if f == floor and s == "novel_species" and b / 2.0 >= t)
            cells.append("%10d %10.2f%%" % (tot, 100.0 * sg / tot) if tot
                         else "%10s %11s" % ("-", "-"))
        p("%-10.1f %s" % (t, "  ".join(cells)))
    p("")
    p("queries with no usable ANI neighbour, by floor and scenario:")
    for floor in floors:
        for name in SCENARIOS:
            p("  af >= %-4g %-16s %6d silent, %6d with a call"
              % (floor, name, silent[(floor, name)], len(best[(floor, name)])))

    # --- the two axes on the same queries --------------------------------------
    if a.mash_calls:
        mash = {}
        with open(a.mash_calls) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["scenario"] == "novel_species":
                    mash[r["accession"]] = (float(r["dist"]), r["neighbour"],
                                            int(r["same_genus"]))
        ani_best = best[(primary, "novel_species")]
        both = set(mash) & set(ani_best)
        agree = sum(1 for q in both if mash[q][1] == ani_best[q][1])
        p("")
        p("Both axes, scenario novel_species, on the queries both can answer:")
        p("  queries answered by Mash: %d" % len(mash))
        p("  queries answered by ANI : %d" % len(ani_best))
        p("  answered by both        : %d" % len(both))
        p("  same nearest neighbour  : %d (%.2f %%)"
          % (agree, 100.0 * agree / len(both) if both else 0.0))
        p("  Mash-only queries       : %d, of which same genus %.2f %%"
          % (len(set(mash) - set(ani_best)),
             100.0 * sum(mash[q][2] for q in set(mash) - set(ani_best))
             / max(1, len(set(mash) - set(ani_best)))))
        p("  on the shared queries, precision of the call at the same coverage:")
        m_ok = sum(mash[q][2] for q in both)
        a_ok = sum(1 for q in both if ani_best[q][2] in SAME_GENUS)
        p("    Mash %.2f %%    ANI %.2f %%"
          % (100.0 * m_ok / len(both), 100.0 * a_ok / len(both)))

        # Disagreeing on which genome is nearest is not the same as disagreeing
        # on the genus, and the two are worth separating.
        diff = [q for q in both if mash[q][1] != ani_best[q][1]]
        verdict = sum(1 for q in diff
                      if mash[q][2] == int(ani_best[q][2] in SAME_GENUS))
        p("  where the two pick a different neighbour (%d queries), they still"
          % len(diff))
        p("    agree on the genus verdict in %d of them (%.2f %%)"
          % (verdict, 100.0 * verdict / len(diff) if diff else 0.0))

        # Where the ANI axis is silent, what was Mash's own call worth?
        only = [q for q in mash if q not in ani_best]
        if only:
            dd = sorted(mash[q][0] for q in only)
            p("  the %d queries only Mash answers sit at Mash distance"
              % len(only))
            p("    q1 %.3f  median %.3f  q3 %.3f — beyond any usable cutoff"
              % (dd[len(dd) // 4], dd[len(dd) // 2], dd[(3 * len(dd)) // 4]))

        # The comparison that matters: at the same number of calls, which index
        # is more precise? A cutoff is only interpretable next to its coverage.
        p("")
        p("Matched call volume, scenario novel_species:")
        p("%-12s %8s %10s   %s" % ("Mash cutoff", "calls", "precision",
                                   "  ".join("%-26s" % ("ANI, af >= %g" % f)
                                             for f in floors)))
        for cut in (0.05, 0.07, 0.08, 0.10, 0.12, 0.13, 0.15, 0.20):
            sel = [q for q in mash if mash[q][0] <= cut]
            if not sel:
                continue
            m_p = 100.0 * sum(mash[q][2] for q in sel) / len(sel)
            cells = []
            for floor in floors:
                bins = sorted((b for (f, s, b) in curve
                               if f == floor and s == "novel_species"), reverse=True)
                run_t = run_s = 0
                pick = None
                for b in bins:
                    sg, tot = curve[(floor, "novel_species", b)]
                    run_t += tot
                    run_s += sg
                    if pick is None or abs(run_t - len(sel)) < abs(pick[1] - len(sel)):
                        pick = (b / 2.0, run_t, 100.0 * run_s / run_t)
                cells.append("ANI>=%.1f %7d %7.2f%%" % pick)
            p("%-12.2f %8d %9.2f%%   %s" % (cut, len(sel), m_p, "  ".join(cells)))

    text = "\n".join(out)
    print(text)
    with open(a.out + "_summary.txt", "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
