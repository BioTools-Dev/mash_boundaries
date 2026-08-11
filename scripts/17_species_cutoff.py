#!/usr/bin/env python3
"""Where the species cutoff goes, in the leave-one-out view.

Sections 5.7 to 5.13 all answer the genus question. The species cutoff has so
far been taken from the literature — 95 % ANI, converted to Mash distance as
1 - ANI — and never measured in the condition a tool operates in: a query
against a database, where what gets evaluated is the minimum over the database
and not a random pair.

The scenario is `no_self_strain`: replicate assemblies of the query's own strain
are removed, because a type-material collection holds up to twelve assemblies of
the same strain and a hit against one of those measures nothing. Every other
genome stays, so a query whose species is represented can find it.

Two quantities are reported at each cutoff, and they pull in opposite directions:

  precision  of the queries called at or below the cutoff, the fraction whose
             neighbour really is the same species. This is what a false species
             claim costs.
  recall     of the queries that *have* a conspecific genome in the database at
             all, the fraction the cutoff actually calls. A cutoff that is
             precise because it calls almost nothing is not a good cutoff, and
             the pair of numbers is the only honest way to show that.

The errors are broken out by what they really are: a call that is congeneric is
the right genus and the wrong species, which is a mild error; one that is
inter-genus is a different kind of failure.
"""
import argparse
import collections
import csv

# The console table stays readable; the TSV carries a dense grid so the curve can
# be drawn without re-running the analysis at every point.
REPORT = (0.02, 0.03, 0.04, 0.0426, 0.045, 0.05, 0.06, 0.07, 0.08, 0.10)
CUTOFFS = sorted(set(REPORT) | {round(0.002 * i, 4) for i in range(1, 76)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="the _calls.tsv written by step 09")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--relabel", default="",
                    help="accession/species_key/genus_key table; when given, the "
                         "denominator of recall is recomputed under those labels "
                         "too, since a second taxonomy merges and splits species")
    ap.add_argument("--scenario", default="no_self_strain")
    ap.add_argument("--domain", default="Bacteria")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    lab = {r["accession"]: r for r in csv.DictReader(open(a.labels), delimiter="\t")}

    # how many genomes of the query's own species remain once its own strain is
    # taken out — the denominator of recall
    if a.relabel:
        key = {}
        with open(a.relabel) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["species_key"]:
                    key[r["accession"]] = r["species_key"]
        by_species = collections.Counter(key.values())
        by_strain = collections.Counter()          # strain identity stays taxid-based
        for acc, r in lab.items():
            if acc in key and r["taxid"] != r["species_taxid"]:
                by_strain[r["taxid"]] += 1

        def available(acc):
            if acc not in key:
                return None
            n = by_species[key[acc]] - 1
            r = lab[acc]
            if r["taxid"] != r["species_taxid"]:
                n -= by_strain[r["taxid"]] - 1
            return n
    else:
        by_species = collections.Counter(r["species_taxid"] for r in lab.values())
        by_strain = collections.Counter(r["taxid"] for r in lab.values()
                                        if r["taxid"] != r["species_taxid"])

        def available(acc):
            r = lab[acc]
            n = by_species[r["species_taxid"]] - 1
            if r["taxid"] != r["species_taxid"]:
                n -= by_strain[r["taxid"]] - 1
            return n

    rows = []
    with open(a.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["scenario"] != a.scenario or r["domain"] != a.domain:
                continue
            rows.append((float(r["dist"]), r["class"], r["accession"]))

    queried = {acc for _, _, acc in rows}
    with_conspecific = sum(1 for acc in queried
                           if (available(acc) or 0) > 0)

    out = []
    p = out.append
    p("scenario %s, domain %s" % (a.scenario, a.domain))
    p("queries with a usable neighbour: %d" % len(rows))
    p("of those, queries whose species is represented by another strain: %d"
      % with_conspecific)
    p("")
    p("%-8s %9s %11s %10s %12s %12s" % ("d <=", "llamadas", "precisión", "cobertura",
                                        "err. género", "err. otro gén."))
    table = []
    for D in CUTOFFS:
        sel = [c for d, c, _ in rows if d <= D]
        if not sel:
            continue
        n = len(sel)
        ok = sum(1 for c in sel if c in ("same_species", "same_strain"))
        cong = sum(1 for c in sel if c == "congeneric")
        inter = sum(1 for c in sel if c == "inter_genus")
        table.append((D, n, ok / n, ok / with_conspecific, cong, inter))
        if D in REPORT:
            p("%-8.4f %9d %10.2f%% %9.2f%% %12d %12d"
              % (D, n, 100.0 * ok / n, 100.0 * ok / with_conspecific, cong, inter))

    with open(a.out + ".tsv", "w") as fh:
        fh.write("cutoff\tcalls\tprecision\trecall\terr_congeneric\terr_inter_genus\n")
        for D, n, pr, rc, cong, inter in table:
            fh.write("%.4f\t%d\t%.6f\t%.6f\t%d\t%d\n" % (D, n, pr, rc, cong, inter))

    text = "\n".join(out)
    print(text)
    with open(a.out + "_summary.txt", "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
