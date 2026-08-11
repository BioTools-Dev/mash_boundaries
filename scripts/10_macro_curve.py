#!/usr/bin/env python3
"""Genus-weighted (macro) version of the leave-one-out precision curve.

The micro curve gives every query genome one vote, so a genus with 775 sequenced
species contributes 775 times the weight of a genus with one. Since *Streptomyces*
alone accounts for 37.5 % of all congeneric pairs, a curve read off the raw
queries is substantially a curve of the few largest genera.

The macro curve computes the precision separately inside each genus and then
averages the genera with equal weight. If the two agree, the imbalance did not
matter and that is worth stating; if they diverge, the divergence is the finding,
and it would explain why published thresholds do not transfer between studies
that sampled the tree differently.

Genera contribute only above a minimum number of calls, since a genus with two
queries yields a precision of 0, 50 or 100 % and nothing in between. The
threshold is a declared parameter, not a tuned one, and the number of genera
retained at each cutoff is reported alongside the estimate.

The spread across genera is reported as quartiles rather than a standard error:
the per-genus precisions are bounded and strongly skewed, so a symmetric interval
would misdescribe them. A percentile bootstrap over genera gives the interval of
the macro mean itself.
"""
import argparse
import collections
import csv
import random
import statistics
import sys

CUTOFFS = (0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.13, 0.15, 0.18, 0.20, 0.25)
# The cutoffs above are on the DNA scale; an axis with a different alphabet
# needs its own, so they can be replaced from the command line.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="the _calls.tsv written by step 09")
    ap.add_argument("--out", required=True)
    ap.add_argument("--scenario", default="novel_species")
    ap.add_argument("--min-calls", type=int, default=5,
                    help="least calls a genus needs to enter the macro average")
    ap.add_argument("--boot", type=int, default=2000, help="bootstrap resamples over genera")
    ap.add_argument("--domain", default="Bacteria")
    ap.add_argument("--cutoffs", default="",
                    help="comma-separated cutoffs, for an axis whose scale is "
                         "not the one of the DNA sketch")
    a = ap.parse_args()

    # per genus: list of (distance, same_genus) for the queries of that genus
    per_genus = collections.defaultdict(list)
    with open(a.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["scenario"] != a.scenario or r["domain"] != a.domain:
                continue
            per_genus[r["genus"]].append((float(r["dist"]), r["same_genus"] == "1"))

    if not per_genus:
        sys.exit("ABORT: no rows for scenario %s / domain %s" % (a.scenario, a.domain))

    rng = random.Random(42)
    cutoffs = ([float(x) for x in a.cutoffs.split(",")]
               if a.cutoffs else CUTOFFS)
    rows = []
    for D in cutoffs:
        micro_hit = micro_tot = 0
        pg = []
        for g, lst in per_genus.items():
            hit = sum(1 for d, s in lst if d <= D and s)
            tot = sum(1 for d, _ in lst if d <= D)
            micro_hit += hit
            micro_tot += tot
            if tot >= a.min_calls:
                pg.append(hit / tot)
        if not pg or not micro_tot:
            continue

        macro = statistics.fmean(pg)
        q1, med, q3 = (statistics.quantiles(pg, n=4) if len(pg) > 3
                       else (min(pg), statistics.median(pg), max(pg)))
        boots = sorted(statistics.fmean(rng.choices(pg, k=len(pg))) for _ in range(a.boot))
        lo = boots[int(0.025 * len(boots))]
        hi = boots[int(0.975 * len(boots)) - 1]

        rows.append(dict(cutoff=D, identity=1 - D, n_genera=len(pg), n_calls=micro_tot,
                         micro=micro_hit / micro_tot, macro=macro,
                         boot_lo=lo, boot_hi=hi, q1=q1, median=med, q3=q3))

    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print("scenario %s, %s, genera with >= %d calls\n"
          % (a.scenario, a.domain, a.min_calls))
    print("%-7s %8s %8s %9s %9s %19s %23s"
          % ("cutoff", "genera", "calls", "micro", "macro", "macro 95 % CI", "per-genus Q1/med/Q3"))
    for r in rows:
        print("%-7.2f %8d %8d %8.2f%% %8.2f%%   [%6.2f%%, %6.2f%%]   %6.1f%% %6.1f%% %6.1f%%"
              % (r["cutoff"], r["n_genera"], r["n_calls"], 100 * r["micro"], 100 * r["macro"],
                 100 * r["boot_lo"], 100 * r["boot_hi"],
                 100 * r["q1"], 100 * r["median"], 100 * r["q3"]))


if __name__ == "__main__":
    main()
