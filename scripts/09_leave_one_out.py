#!/usr/bin/env python3
"""Leave-one-out view: P(same genus | the nearest neighbour is at distance d).

This is the quantity a tool actually operates on. It does not evaluate a random
pair of genomes — it evaluates the *minimum* over a reference database, and the
base rate of that minimum bears no relation to the base rate over all pairs. The
pairwise view answers "how does Mash distance relate to taxonomy"; this one
answers "where should a cutoff go".

Three scenarios are computed, each removing more of the query's own taxonomy from
the database, because they answer different questions:

  as_reported   the nearest neighbour, whatever it is. Includes replicate
                assemblies of the same strain, so it overstates how close a real
                query lands.
  no_self_strain
                replicate assemblies of the query's own strain are removed from
                the database. This is the ordinary case: the genome is new, its
                species may or may not be represented.
  novel_species
                every genome of the query's own species is removed, so the
                nearest possible hit is a different species. This is the
                out-of-domain case a genus cutoff exists to serve, and the one
                the mirror experiment probes at n=14.

The neighbour table keeps only the k nearest per genome, so a query whose k
neighbours are all excluded by a scenario has no usable neighbour left. Those are
counted and reported rather than silently dropped: they are a truncation
artefact, not data.
"""
import argparse
import collections
import csv
import os
import sys

SAME_GENUS = ("same_strain", "same_species", "congeneric")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nn", required=True, help="<prefix>_nn.tsv from the accumulator")
    ap.add_argument("--out", required=True, help="output TSV of the binned curves")
    ap.add_argument("--calls", default="",
                    help="if set, write one row per query and scenario here; the "
                         "genus-weighted macro curve and the conflict listing are "
                         "both built from this table rather than recomputing the "
                         "scenario logic")
    ap.add_argument("--bin", type=float, default=0.005)
    ap.add_argument("--relabel", default="",
                    help="TSV with accession/species_key/genus_key columns. When "
                         "given, the class of every pair is recomputed from those "
                         "keys instead of the one the accumulator wrote, which is "
                         "how the same distances are read against a second "
                         "taxonomy. A genome absent from the table cannot be "
                         "classified and is removed from the database, exactly as "
                         "a scenario exclusion would remove it.")
    ap.add_argument("--domain", default="",
                    help="restrict the queries to one domain. Neighbours are left "
                         "alone: a real database holds both domains, and a query "
                         "whose best hit is in the other one is a result, not an "
                         "exclusion.")
    ap.add_argument("--labels", default="",
                    help="genome_labels.tsv; needed only with --stratum")
    ap.add_argument("--stratum", default="", choices=["", "gold", "quality", "taxcheck"],
                    help="restrict both the queries and the genomes they may match "
                         "to a quality stratum, which is the leave-one-out "
                         "equivalent of restricting the reference database. The "
                         "three gates are kept apart on purpose: `gold` is the full "
                         "gate, `quality` is the CheckM criterion alone, and "
                         "`taxcheck` is NCBI's own ANI-based check of the declared "
                         "name. Only the last one is entangled with what is being "
                         "measured, so reporting them separately is what tells an "
                         "improvement in assembly quality apart from an improvement "
                         "obtained by discarding the genomes whose names were "
                         "already suspect. Queries whose k kept neighbours all fall "
                         "outside the stratum lose them: with a truncated neighbour "
                         "table the restriction cannot invent a farther hit that was "
                         "never stored, and the loss is reported rather than hidden.")
    a = ap.parse_args()

    GATES = {"gold": ("pass_gold", "1"),
             "quality": ("low_quality", "0"),
             "taxcheck": ("taxcheck_bad", "0")}
    keep_set = set()
    if a.stratum:
        if not a.labels:
            sys.exit("ABORT: --stratum needs --labels")
        col, want = GATES[a.stratum]
        with open(a.labels) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r[col] == want:
                    keep_set.add(r["accession"])

    relabel = {}
    if a.relabel:
        with open(a.relabel) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["species_key"] and r["genus_key"]:
                    relabel[r["accession"]] = (r["species_key"], r["genus_key"])

    def reclass(q, n, original):
        """Class of a pair under the replacement labels."""
        if original == "same_strain":
            return original                      # strain identity is not taxonomic
        qs, qg = relabel[q]
        ns, ng = relabel[n]
        if qs == ns:
            return "same_species"
        return "congeneric" if qg == ng else "inter_genus"

    # neighbours of each genome, already ordered by rank
    nbr = collections.defaultdict(list)
    meta = {}
    dropped = 0
    off_stratum = 0
    with open(a.nn) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, n, cls = r["accession"], r["neighbour"], r["class"]
            genus = r["genus"]
            if a.domain and r["domain"] != a.domain:
                continue
            if keep_set and (q not in keep_set or n not in keep_set):
                if q in keep_set:
                    off_stratum += 1          # the query stays, this neighbour goes
                continue
            if relabel:
                if q not in relabel or n not in relabel:
                    dropped += 1
                    continue
                cls = reclass(q, n, cls)
                genus = relabel[q][1]
            nbr[q].append((float(r["dist"]), cls, n))
            meta[q] = (genus, r["domain"])
    if relabel:
        print("relabelled with %s: %d neighbours dropped for want of a label\n"
              % (os.path.basename(a.relabel), dropped))
    if keep_set:
        print("stratum %s: %d genomes pass; %d neighbours of a passing query were "
              "discarded for falling outside it\n"
              % (a.stratum, len(keep_set), off_stratum))
    if a.domain:
        print("queries restricted to %s: %d with at least one neighbour\n"
              % (a.domain, len(nbr)))

    scenarios = {
        "as_reported": lambda c: True,
        "no_self_strain": lambda c: c != "same_strain",
        "novel_species": lambda c: c not in ("same_strain", "same_species"),
    }

    # counts[(scenario, bin)] = [same-genus hits, total]
    counts = collections.defaultdict(lambda: [0, 0])
    exhausted = collections.Counter()
    used = collections.Counter()

    calls = None
    if a.calls:
        calls = open(a.calls, "w")
        calls.write("scenario\taccession\tgenus\tdomain\tneighbour\tdist\tclass\tsame_genus\n")

    for acc, lst in nbr.items():
        for name, keep in scenarios.items():
            hit = next(((d, c, n) for d, c, n in lst if keep(c)), None)
            if hit is None:
                exhausted[name] += 1
                continue
            d, c, part = hit
            b = int(d / a.bin)
            counts[(name, b)][1] += 1
            same = c in SAME_GENUS
            if same:
                counts[(name, b)][0] += 1
            used[name] += 1
            if calls is not None:
                g, dm = meta[acc]
                calls.write("%s\t%s\t%s\t%s\t%s\t%.6f\t%s\t%d\n"
                            % (name, acc, g, dm, part, d, c, int(same)))

    if calls is not None:
        calls.close()

    with open(a.out, "w") as fh:
        fh.write("scenario\tdist_lo\tdist_hi\tsame_genus\ttotal\tp_same_genus\tcum_total\tcum_p\n")
        for name in scenarios:
            bins = sorted(b for (s, b) in counts if s == name)
            ch = ct = 0
            for b in bins:
                sg, tot = counts[(name, b)]
                ch += sg
                ct += tot
                fh.write("%s\t%.4f\t%.4f\t%d\t%d\t%.6f\t%d\t%.6f\n"
                         % (name, b * a.bin, (b + 1) * a.bin, sg, tot,
                            sg / tot, ct, ch / ct))

    # --- console summary at the cutoffs the tool uses -------------------------
    print("genomes with at least one neighbour: %d\n" % len(nbr))
    edges = [0.0, 0.05, 0.08, 0.13, 0.20, 1.0]   # 0.13 is the genus cutoff
    print("%-16s %8s %8s   %s" % ("scenario", "usable", "no nbr",
                                  "  ".join("%-18s" % ("[%.2f,%.2f)" % (edges[i], edges[i + 1]))
                                            for i in range(len(edges) - 1))))
    for name in scenarios:
        cells = []
        for i in range(len(edges) - 1):
            lo, hi = edges[i], edges[i + 1]
            sg = sum(counts[(n, b)][0] for (n, b) in counts
                     if n == name and lo <= b * a.bin < hi)
            tt = sum(counts[(n, b)][1] for (n, b) in counts
                     if n == name and lo <= b * a.bin < hi)
            cells.append("%7d %7.2f%%" % (tt, 100.0 * sg / tt) if tt else "%7s %7s" % ("-", "-"))
        print("%-16s %8d %8d   %s" % (name, used[name], exhausted[name],
                                      "  ".join("%-18s" % c for c in cells)))
    print("\nEach cell is: queries whose nearest usable neighbour falls in that")
    print("distance band, and the percentage of those whose neighbour is the same genus.")


if __name__ == "__main__":
    main()
