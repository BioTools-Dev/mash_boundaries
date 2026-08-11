#!/usr/bin/env python3
"""Count how many genome pairs fall in each taxonomic class, per stratum.

The class of a pair is decided on rank taxids, never on names: same species,
same genus but different species, or different genus. Pairs of assemblies that
share a below-species taxid are the same strain and are reported separately —
they are not independent within-species observations.

These counts fix the base rate the whole study has to be read against, so they
are produced before any distance is computed.
"""
import collections
import csv
import itertools
import os
import sys

data = os.environ.get("DATA")
if not data:
    sys.exit("ABORT: source config.sh first")

rows = list(csv.DictReader(open(os.path.join(data, "genome_labels.tsv")), delimiter="\t"))

STRATA = (
    ("all prokaryotes", lambda r: True),
    ("bacteria", lambda r: r["domain"] == "Bacteria"),
    ("bacteria, gold", lambda r: r["domain"] == "Bacteria" and r["pass_gold"] == "1"),
    ("archaea", lambda r: r["domain"] == "Archaea"),
    ("archaea, gold", lambda r: r["domain"] == "Archaea" and r["pass_gold"] == "1"),
)


def census(sub):
    """Pair counts by class, computed from group sizes rather than enumeration."""
    n = len(sub)
    total = n * (n - 1) // 2

    by_strain = collections.Counter(r["taxid"] for r in sub if r["strain_dup"] == "1")
    same_strain = sum(v * (v - 1) // 2 for v in by_strain.values())

    by_species = collections.Counter(r["species_taxid"] for r in sub)
    same_species = sum(v * (v - 1) // 2 for v in by_species.values())

    by_genus = collections.Counter(r["genus_taxid"] for r in sub)
    same_genus = sum(v * (v - 1) // 2 for v in by_genus.values())

    return dict(n=n, species=len(by_species), genera=len(by_genus), total=total,
                same_strain=same_strain,
                same_species=same_species - same_strain,
                congeneric=same_genus - same_species,
                inter_genus=total - same_genus)


hdr = ("%-16s %7s %7s %6s %14s %11s %11s %12s %10s"
       % ("stratum", "genomes", "species", "genera", "pairs",
          "same strain", "same sp.", "congeneric", "inter-gen."))
print(hdr)
print("-" * len(hdr))
for name, keep in STRATA:
    c = census([r for r in rows if keep(r)])
    print("%-16s %7d %7d %6d %14d %11d %11d %12d %10d"
          % (name, c["n"], c["species"], c["genera"], c["total"],
             c["same_strain"], c["same_species"], c["congeneric"], c["inter_genus"]))
    if c["congeneric"]:
        print("%-16s %s" % ("", "  base rate inter-genus : congeneric = %.0f : 1"
                            % (c["inter_genus"] / c["congeneric"])))

# genus size distribution drives how much of the congeneric signal a few large
# genera contribute; a boundary fitted on Streptomyces alone is not a boundary.
bac = [r for r in rows if r["domain"] == "Bacteria"]
gsp = collections.defaultdict(set)
for r in bac:
    gsp[r["genus"]].add(r["species_taxid"])
sizes = sorted(((len(v), g) for g, v in gsp.items()), reverse=True)
pairs = [(k * (k - 1) // 2, g) for k, g in sizes]
tot = sum(p for p, _ in pairs)
print("\nbacterial genera: %d (%d with >= 2 species)"
      % (len(sizes), sum(1 for k, _ in sizes if k > 1)))
print("congeneric species-pairs contributed by the 10 largest genera: %.1f %%"
      % (100.0 * sum(p for p, _ in pairs[:10]) / tot))
for p, g in pairs[:10]:
    print("  %-18s %3d species  %8d pairs  (%4.1f %%)"
          % (g, dict((b, a) for a, b in sizes)[g], p, 100.0 * p / tot))
