#!/usr/bin/env python3
"""Attach the GTDB lineage to every genome of the collection.

GTDB is the second taxonomy the study reports against (D2). It matters here for
one specific reason: §3.7 showed that the residual error of the leave-one-out
curves is dominated by pairs of genomes at over 95 % identity carrying different
genus names — recent, and in several cases contested, generic splits. GTDB
reclassifies most of those, so the same curves computed on GTDB labels are a
falsifiable prediction, not a second opinion.

GTDB keys a genome by its RefSeq accession when one exists (`RS_GCF_...`) and by
its GenBank accession otherwise (`GB_GCA_...`). The collection is indexed by
GenBank accession, so the RefSeq accession carried in the assembly report is what
makes most of the lookups succeed; both are tried and which one matched is
recorded.

Coverage is never complete: GTDB excludes genomes that fail its quality criteria,
so some type material has no GTDB lineage at all. Those genomes are kept in the
table with empty fields and counted, so that the comparison in step 13 can be
restricted to genomes present in both taxonomies rather than quietly changing
denominator.
"""
import argparse
import collections
import csv
import gzip
import os
import sys

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
PREFIX = dict(zip("dpcofgs", RANKS))


def parse_lineage(s):
    out = {}
    for part in s.strip().split(";"):
        if len(part) > 3 and part[1:3] == "__":
            out[PREFIX.get(part[0], part[0])] = part[3:]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--gtdb-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    gtdb = {}
    for name in ("bac120_taxonomy.tsv.gz", "ar53_taxonomy.tsv.gz"):
        path = os.path.join(a.gtdb_dir, name)
        if not os.path.exists(path):
            sys.exit("ABORT: %s missing" % path)
        with gzip.open(path, "rt") as fh:
            for line in fh:
                acc, lineage = line.rstrip("\n").split("\t")
                gtdb[acc.split("_", 1)[1]] = lineage      # drop the RS_/GB_ prefix

    version = "unknown"
    vpath = os.path.join(a.gtdb_dir, "VERSION.txt")
    if os.path.exists(vpath):
        version = open(vpath).readline().strip()
    print("GTDB %s: %d genomes in the taxonomy tables" % (version, len(gtdb)))

    rows = list(csv.DictReader(open(a.labels), delimiter="\t"))
    matched = collections.Counter()
    disagree = collections.Counter()

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "matched_by"] + ["gtdb_" + r for r in RANKS]
                   + ["ncbi_genus", "ncbi_species", "genus_agrees", "species_agrees"])
        for r in rows:
            lineage, how = None, "none"
            if r["accession"] in gtdb:
                lineage, how = gtdb[r["accession"]], "gca"
            elif r["paired_gcf"] and r["paired_gcf"] in gtdb:
                lineage, how = gtdb[r["paired_gcf"]], "gcf"
            matched[how] += 1

            lin = parse_lineage(lineage) if lineage else {}
            ga = sa = ""
            if lin:
                ga = int(lin.get("genus", "") == r["genus"])
                sa = int(lin.get("species", "") == r["species"])
                disagree["genus"] += 1 - ga
                disagree["species"] += 1 - sa
                disagree["n"] += 1
            w.writerow([r["accession"], how] + [lin.get(k, "") for k in RANKS]
                       + [r["genus"], r["species"], ga, sa])

    print("\nlookup: %d by GenBank accession, %d by RefSeq accession, %d not in GTDB"
          % (matched["gca"], matched["gcf"], matched["none"]))
    print("coverage: %.2f %% of %d assemblies"
          % (100.0 * (len(rows) - matched["none"]) / len(rows), len(rows)))
    if disagree["n"]:
        print("\namong the %d matched:" % disagree["n"])
        print("  genus name differs from NCBI:   %6d  (%5.2f %%)"
              % (disagree["genus"], 100.0 * disagree["genus"] / disagree["n"]))
        print("  species name differs from NCBI: %6d  (%5.2f %%)"
              % (disagree["species"], 100.0 * disagree["species"] / disagree["n"]))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
