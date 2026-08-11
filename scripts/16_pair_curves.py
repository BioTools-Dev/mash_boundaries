#!/usr/bin/env python3
"""The pairwise view: P(same species | d) and P(same genus | d) over all pairs.

This is the other half of D8, and the one that carries the base rate. A random
pair of type-material bacteria is inter-genus 239 times out of 240, so a curve
read over all pairs answers a different question from the leave-one-out view of
§3.4: not "where should a cutoff go for a query against a database" but "what
does a Mash distance, on its own, say about two genomes".

Both halves of the curve are computed:

  per bin        P(class | d) among the pairs whose distance falls in the bin.
                 This is the quantity the literature's "between 80 and 95 %
                 identity they almost always share a genus" is a claim about,
                 and where it crosses one half is the sharpest single statement
                 of where the boundary lies in the pairwise view.
  cumulative     the same over every pair at or below d, which is what a cutoff
                 applied to a pair list would actually deliver.

The genus-resolved histogram gives the macro version: the curve is computed
inside each genus and the genera are then averaged with equal weight, because
the ten largest genera contribute 70 % of all congeneric pairs and *Streptomyces*
alone 37.5 % (§3.1). Archaea are averaged separately, per D3.

The counts come from the accumulator's histograms, so nothing here re-reads the
pair stream. Pairs above the screening cutoff are absent by construction; the
census totals are recomputed here to report what fraction of each class the
window contains, rather than leaving the truncation implicit.
"""
import argparse
import collections
import csv
import random
import statistics
import sys

CLASSES = ("same_strain", "same_species", "congeneric", "inter_genus")
SAME_SPECIES = ("same_strain", "same_species")
SAME_GENUS = ("same_strain", "same_species", "congeneric")
BANDS = (0.0, 0.05, 0.08, 0.13, 0.15, 0.20, 0.25, 0.40)


def census(rows):
    """Pair counts by class from group sizes, as in step 04."""
    n = len(rows)
    total = n * (n - 1) // 2
    by_strain = collections.Counter(r["taxid"] for r in rows if r["strain_dup"] == "1")
    by_species = collections.Counter(r["species_taxid"] for r in rows)
    by_genus = collections.Counter(r["genus_taxid"] for r in rows)
    ss = sum(v * (v - 1) // 2 for v in by_strain.values())
    sp = sum(v * (v - 1) // 2 for v in by_species.values())
    gn = sum(v * (v - 1) // 2 for v in by_genus.values())
    return {"same_strain": ss, "same_species": sp - ss, "congeneric": gn - sp,
            "inter_genus": total - gn, "total": total}


def curve_rows(counts, bins):
    """Per-bin and cumulative probabilities from {(bin, class): count}."""
    cum = collections.Counter()
    out = []
    for b in bins:
        per = {c: counts.get((b, c), 0) for c in CLASSES}
        tot = sum(per.values())
        if not tot:
            continue
        for c in CLASSES:
            cum[c] += per[c]
        ctot = sum(cum.values())
        out.append((b, per, tot,
                    sum(per[c] for c in SAME_SPECIES) / tot,
                    sum(per[c] for c in SAME_GENUS) / tot,
                    ctot,
                    sum(cum[c] for c in SAME_SPECIES) / ctot,
                    sum(cum[c] for c in SAME_GENUS) / ctot))
    return out


def crossing(rows, idx, bw):
    """First bin, going outward, at which a per-bin probability drops below 0.5."""
    prev = None
    for b, per, tot, ps, pg, ctot, cps, cpg in rows:
        p = (ps, pg)[idx]
        if prev is not None and prev >= 0.5 > p:
            return b * bw, p
        prev = p
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist", required=True, help="<prefix>_hist.tsv")
    ap.add_argument("--hist-genus", default="", help="<prefix>_hist_genus.tsv")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--bin", type=float, default=0.005)
    ap.add_argument("--min-pairs", type=int, default=20,
                    help="least pairs a genus needs in a bin to enter the macro average")
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()

    lab = list(csv.DictReader(open(a.labels), delimiter="\t"))
    dom_of_genus = {r["genus_taxid"]: r["domain"] for r in lab}
    name_of_genus = {r["genus_taxid"]: r["genus"] for r in lab}

    # --- the pairwise curve, pooled and on the gold stratum -------------------
    pooled = collections.Counter()
    gold = collections.Counter()
    with open(a.hist) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            b = int(round(float(r["dist_lo"]) / a.bin - 0.5))
            n = int(r["count"])
            pooled[(b, r["class"])] += n
            if r["gold_pair"] == "1":
                gold[(b, r["class"])] += n

    bins = sorted({b for (b, _) in pooled})
    rows_pooled = curve_rows(pooled, bins)
    rows_gold = curve_rows(gold, sorted({b for (b, _) in gold}))

    with open(a.out + "_pairs.tsv", "w") as fh:
        fh.write("stratum\tdist_lo\tdist_hi\t%s\ttotal\tp_same_species\tp_same_genus\t"
                 "cum_total\tcum_p_same_species\tcum_p_same_genus\n" % "\t".join(CLASSES))
        for name, rows in (("all", rows_pooled), ("gold", rows_gold)):
            for b, per, tot, ps, pg, ctot, cps, cpg in rows:
                fh.write("%s\t%.4f\t%.4f\t%s\t%d\t%.6f\t%.6f\t%d\t%.6f\t%.6f\n"
                         % (name, b * a.bin, (b + 1) * a.bin,
                            "\t".join(str(per[c]) for c in CLASSES),
                            tot, ps, pg, ctot, cps, cpg))

    # --- macro over genera, by domain -----------------------------------------
    macro_rows, macro_bands = [], []
    if a.hist_genus:
        per_genus = collections.defaultdict(collections.Counter)   # genus -> (bin, class)
        with open(a.hist_genus) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                b = int(round(float(r["dist_lo"]) / a.bin - 0.5))
                per_genus[r["genus_taxid"]][(b, r["class"])] += int(r["count"])

        rng = random.Random(42)

        def macro(genera, sel_bins):
            """Equal-weight mean over genera of P(same genus | d in the range)."""
            vals, hit, tot = [], 0, 0
            for g in genera:
                c = per_genus[g]
                n = sum(c.get((b, k), 0) for b in sel_bins for k in CLASSES)
                if n < a.min_pairs:
                    continue
                s = sum(c.get((b, k), 0) for b in sel_bins for k in SAME_GENUS)
                vals.append(s / n)
                hit += s
                tot += n
            if len(vals) < 2:
                return None
            boot = sorted(statistics.mean(rng.choices(vals, k=len(vals)))
                          for _ in range(a.boot))
            q = sorted(vals)
            return (len(vals), tot, hit / tot, statistics.mean(vals),
                    boot[int(0.025 * a.boot)], boot[int(0.975 * a.boot)],
                    q[len(q) // 4], statistics.median(q), q[(3 * len(q)) // 4])

        for domain in ("Bacteria", "Archaea"):
            genera = [g for g in per_genus if dom_of_genus.get(g) == domain]
            for b in bins:
                m = macro(genera, (b,))
                if m:
                    macro_rows.append((domain, b * a.bin) + m)
            # A single 0.005 bin leaves few genera above the minimum, so the
            # bands the rest of the study uses are averaged as well.
            for i in range(len(BANDS) - 1):
                lo, hi = BANDS[i], BANDS[i + 1]
                m = macro(genera, [b for b in bins if lo <= b * a.bin < hi])
                if m:
                    macro_bands.append((domain, lo, hi) + m)

        with open(a.out + "_pairs_macro.tsv", "w") as fh:
            fh.write("domain\tdist_lo\tgenera\tpairs\tmicro_p\tmacro_p\tci_lo\tci_hi\t"
                     "q1\tmedian\tq3\n")
            for r in macro_rows:
                fh.write("%s\t%.4f\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % r)

        with open(a.out + "_pairs_macro_bands.tsv", "w") as fh:
            fh.write("domain\tdist_lo\tdist_hi\tgenera\tpairs\tmicro_p\tmacro_p\t"
                     "ci_lo\tci_hi\tq1\tmedian\tq3\n")
            for r in macro_bands:
                fh.write("%s\t%.4f\t%.4f\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % r)

    # --- summary ---------------------------------------------------------------
    out = []
    p = out.append
    bac = [r for r in lab if r["domain"] == "Bacteria"]
    cen = census(bac)
    p("Bacterial pair census (recomputed): %s"
      % "  ".join("%s %d" % (c, cen[c]) for c in CLASSES))
    p("")
    p("Pairwise view, all reported pairs, by band:")
    p("%-14s %14s %12s %12s %12s   %10s %10s"
      % ("band", "pairs", "same sp.", "congeneric", "inter-genus",
         "P(sp.|d)", "P(gen.|d)"))
    for i in range(len(BANDS) - 1):
        lo, hi = BANDS[i], BANDS[i + 1]
        sel = [r for r in rows_pooled if lo <= r[0] * a.bin < hi]
        if not sel:
            continue
        per = {c: sum(r[1][c] for r in sel) for c in CLASSES}
        tot = sum(per.values())
        p("%-14s %14d %12d %12d %12d   %9.2f%% %9.2f%%"
          % ("%.2f-%.2f" % (lo, hi), tot,
             per["same_strain"] + per["same_species"], per["congeneric"],
             per["inter_genus"],
             100.0 * sum(per[c] for c in SAME_SPECIES) / tot,
             100.0 * sum(per[c] for c in SAME_GENUS) / tot))
    p("")
    for idx, what in ((0, "same species"), (1, "same genus")):
        d, pv = crossing(rows_pooled, idx, a.bin)
        p("P(%s | d) falls below one half in the bin starting at %s"
          % (what, "%.3f (to %.3f there)" % (d, pv) if d is not None else "no bin"))
    p("")
    p("What a cutoff on a pair list would deliver (cumulative, all reported pairs):")
    p("%-10s %14s %14s %12s" % ("d <=", "pairs", "same genus", "precision"))
    for cut in (0.05, 0.08, 0.12, 0.13, 0.15, 0.20, 0.25):
        sel = [r for r in rows_pooled if r[0] * a.bin < cut]
        if not sel:
            continue
        last = sel[-1]
        p("%-10.2f %14d %14d %11.2f%%"
          % (cut, last[5], round(last[7] * last[5]), 100.0 * last[7]))
    p("")
    all_cen = census(lab)
    p("Window: the reported pairs are %.2f %% of all congeneric prokaryotic pairs"
      % (100.0 * sum(r[1]["congeneric"] for r in rows_pooled) / all_cen["congeneric"]))
    p("")
    if macro_bands:
        p("Macro over genera (equal weight per genus), genera with >= %d pairs in the band:"
          % a.min_pairs)
        p("%-10s %-12s %8s %10s %10s %22s   %s"
          % ("domain", "band", "genera", "micro", "macro", "IC 95 %", "Q1 / med / Q3"))
        for domain, lo, hi, ng, tot, mi, ma, cl, ch, q1, md, q3 in macro_bands:
            p("%-10s %-12s %8d %9.2f%% %9.2f%% [%9.2f%%, %9.2f%%]   %.1f / %.1f / %.1f %%"
              % (domain, "%.2f-%.2f" % (lo, hi), ng, 100 * mi, 100 * ma,
                 100 * cl, 100 * ch, 100 * q1, 100 * md, 100 * q3))

    text = "\n".join(out)
    print(text)
    with open(a.out + "_pairs_summary.txt", "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
