#!/usr/bin/env python3
"""Name the cases where the nearest neighbour of a genome is another genus.

Under the novel_species scenario every genome of the query's own species is
removed from the database, so a wrong genus call means the closest remaining
genome carries a different genus name. Below distance 0.05 that is a strong
statement: two genomes at over 95 % identity, i.e. inside the operational species
boundary, bearing different genus names. One of the two names is very probably
wrong, and the pair is a taxonomic conflict rather than a failure of the metric.

The distinction matters for how the curves of §3.4 are read. A residual error
rate of 4 % is not 4 % of Mash being wrong: an unknown share of it is the
reference taxonomy disagreeing with genomic evidence. Listing the cases is what
turns that from an excuse into a measurement, and the same list is the input to
the GTDB comparison, where most of these pairs should already be reconciled.

Three products:

  <out>_pairs.tsv    every conflicting query, with both names and the distance
  <out>_genera.tsv   the conflicts aggregated by unordered genus pair
  <out>_worst.tsv    the genera whose queries are called wrongly most often
"""
import argparse
import collections
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="the _calls.tsv written by step 09")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--scenario", default="novel_species")
    ap.add_argument("--gtdb", default="",
                    help="accession/species_key/genus_key table from step 12. When "
                         "given, every conflict is checked against that second "
                         "taxonomy and reported as reconciled when it places the two "
                         "genomes in one genus. This is the falsifiable half of the "
                         "claim that the residual error is nomenclatural: if these "
                         "pairs are disagreements about names rather than about "
                         "genomes, a taxonomy built on genomic coherence should have "
                         "merged most of them already.")
    ap.add_argument("--max-dist", type=float, default=0.05,
                    help="distance below which a different-genus hit is treated as "
                         "a taxonomic conflict rather than a distant relative")
    a = ap.parse_args()

    lab = {}
    with open(a.labels) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab[r["accession"]] = (r["genus"], r["species"], r["domain"],
                                   r["taxcheck"], r["organism"])

    gtdb = {}
    if a.gtdb:
        with open(a.gtdb) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["genus_key"]:
                    gtdb[r["accession"]] = r["genus_key"]

    rows = []
    per_genus = collections.Counter()
    per_genus_bad = collections.Counter()
    with open(a.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["scenario"] != a.scenario:
                continue
            per_genus[r["genus"]] += 1
            if r["same_genus"] == "1":
                continue
            per_genus_bad[r["genus"]] += 1
            d = float(r["dist"])
            if d > a.max_dist:
                continue
            qg, qs, qd, qt, qo = lab[r["accession"]]
            ng, ns, _, nt, no = lab.get(r["neighbour"], ("?", "?", "?", "?", "?"))
            row = dict(dist=d, domain=qd,
                       query=r["accession"], query_genus=qg, query_species=qs,
                       neighbour=r["neighbour"], neighbour_genus=ng,
                       neighbour_species=ns,
                       query_taxcheck=qt, neighbour_taxcheck=nt)
            if a.gtdb:
                gq, gn = gtdb.get(r["accession"], ""), gtdb.get(r["neighbour"], "")
                row["gtdb_query_genus"] = gq or "-"
                row["gtdb_neighbour_genus"] = gn or "-"
                row["gtdb_status"] = ("absent" if not (gq and gn)
                                      else "reconciled" if gq == gn else "kept apart")
            rows.append(row)

    rows.sort(key=lambda r: r["dist"])
    with open(a.out + "_pairs.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    pairs = collections.Counter()
    closest = {}
    for r in rows:
        key = tuple(sorted((r["query_genus"], r["neighbour_genus"])))
        pairs[key] += 1
        if key not in closest or r["dist"] < closest[key]:
            closest[key] = r["dist"]
    # what the second taxonomy did with each conflicting genus pair
    verdict = collections.defaultdict(collections.Counter)
    for r in rows:
        if a.gtdb:
            key = tuple(sorted((r["query_genus"], r["neighbour_genus"])))
            verdict[key][r["gtdb_status"]] += 1

    def status_of(key):
        """One word for the pair: what GTDB did with the majority of its cases."""
        if not a.gtdb:
            return ""
        c = verdict[key]
        return c.most_common(1)[0][0] if c else "absent"

    with open(a.out + "_genera.tsv", "w") as fh:
        extra = "\tgtdb_verdict\tgtdb_genus" if a.gtdb else ""
        fh.write("genus_a\tgenus_b\tconflicts\tmin_dist%s\n" % extra)
        for (ga, gb), n in pairs.most_common():
            tail = ""
            if a.gtdb:
                st = status_of((ga, gb))
                # the merged name is only meaningful when the pair was merged
                names = collections.Counter(
                    r["gtdb_query_genus"] for r in rows
                    if tuple(sorted((r["query_genus"], r["neighbour_genus"]))) == (ga, gb)
                    and r["gtdb_query_genus"] != "-") if st == "reconciled" else None
                tail = "\t%s\t%s" % (st, names.most_common(1)[0][0] if names else "-")
            fh.write("%s\t%s\t%d\t%.6f%s\n" % (ga, gb, n, closest[(ga, gb)], tail))

    with open(a.out + "_worst.tsv", "w") as fh:
        fh.write("genus\tqueries\twrong_genus_calls\trate\n")
        worst = sorted(((g, per_genus[g], per_genus_bad[g], per_genus_bad[g] / per_genus[g])
                        for g in per_genus if per_genus[g] >= 5 and per_genus_bad[g]),
                       key=lambda x: (-x[2], -x[3]))
        for g, n, bad, rate in worst:
            fh.write("%s\t%d\t%d\t%.4f\n" % (g, n, bad, rate))

    print("conflicts below distance %.2f: %d queries\n" % (a.max_dist, len(rows)))
    if a.gtdb:
        tot = collections.Counter(r["gtdb_status"] for r in rows)
        pair_tot = collections.Counter(status_of(k) for k in pairs)
        print("against the second taxonomy: %d of %d conflicting queries reconciled, "
              "%d kept apart, %d absent" % (tot["reconciled"], len(rows),
                                            tot["kept apart"], tot["absent"]))
        print("                            %d of %d genus pairs reconciled\n"
              % (pair_tot["reconciled"], len(pairs)))
    print("%-26s %-26s %6s %9s %s" % ("genus A", "genus B", "cases", "min dist",
                                      "GTDB" if a.gtdb else ""))
    for (ga, gb), n in pairs.most_common(20):
        print("%-26s %-26s %6d %9.4f %s" % (ga, gb, n, closest[(ga, gb)],
                                            status_of((ga, gb))))
    print("\nclosest individual conflicts:\n")
    print("%-9s %-34s %-34s" % ("dist", "query", "nearest neighbour"))
    for r in rows[:15]:
        print("%-9.5f %-34s %-34s" % (r["dist"], r["query_species"][:33],
                                      r["neighbour_species"][:33]))
    print("\nwrote %s_{pairs,genera,worst}.tsv" % a.out)


if __name__ == "__main__":
    main()
