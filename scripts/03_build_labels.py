#!/usr/bin/env python3
"""Join assembly metadata with the resolved lineage into the gold-standard table.

Writes $DATA/genome_labels.tsv, one row per assembly, carrying the rank taxids a
pair of genomes is classified on (same species / same genus, different species /
different genus) plus every flag needed to define a stratum later.

Nothing is dropped here. Filtering is a decision of the analysis, not of the
label build, so each disqualifying condition becomes a column and the counts are
reported. The flags are:

  dup_taxid       another assembly shares this taxid
  strain_dup      the shared taxid is below species rank, so the assemblies are
                  the same strain sequenced more than once — such a pair is not
                  an independent within-species observation and must not be
                  counted as one. When the shared taxid IS the species taxid the
                  assemblies are different strains of that species, which is a
                  legitimate within-species pair and is kept.
  no_genus        the lineage has no genus rank (unclassified / candidate taxa)
  name_mismatch   the species from the lineage disagrees with the species NCBI
                  parsed from the submitted name
  taxcheck_bad    NCBI's own ANI check did not return OK
  low_quality     CheckM completeness < 90 % or contamination > 5 %

`pass_gold` is true when none of no_genus, name_mismatch, taxcheck_bad or
low_quality is set. The duplicate flags are deliberately left out of it: those
assemblies are valid genomes and only their mutual pairs need special handling.
"""
import csv
import collections
import os
import sys

data = os.environ.get("DATA")
if not data:
    sys.exit("ABORT: source config.sh first")

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")

# --- lineage keyed by strain taxid -------------------------------------------
lineage = {}
with open(os.path.join(data, "lineage.tsv")) as fh:
    rd = csv.reader(fh, delimiter="\t")
    next(rd)
    for taxid, names, taxids in rd:
        nm = names.split("|")
        tx = taxids.split("|")
        lineage[taxid] = dict(zip(RANKS, nm), **{r + "_taxid": t for r, t in zip(RANKS, tx)})

# --- archaeal accessions, as listed at download time -------------------------
archaea = set()
with open(os.environ["ARCHAEA_LIST"]) as fh:
    for line in fh:
        acc = line.strip().rstrip("/")
        if acc.startswith("GCA_"):
            archaea.add(acc)

rows = list(csv.DictReader(open(os.path.join(data, "genomes_raw.tsv")), delimiter="\t"))
dup = collections.Counter(r["taxid"] for r in rows)


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


out_cols = (["accession", "taxid", "organism", "domain", "domain_listed"]
            + [c for r in RANKS[1:] for c in (r, r + "_taxid")]
            + ["species_name_ncbi", "paired_gcf", "level", "length", "n_contigs",
               "completeness", "contamination", "taxcheck",
               "dup_taxid", "strain_dup", "no_genus", "name_mismatch",
               "taxcheck_bad", "low_quality", "pass_gold"])

FLAG_ORDER = ("dup_taxid", "strain_dup", "no_genus", "name_mismatch",
              "taxcheck_bad", "low_quality", "pass_gold")

flags = collections.Counter()
out_path = os.path.join(data, "genome_labels.tsv")
with open(out_path, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(out_cols)
    for r in rows:
        lin = lineage.get(r["taxid"], {})
        comp, cont = fnum(r["completeness"]), fnum(r["contamination"])

        f = {
            "dup_taxid": dup[r["taxid"]] > 1,
            "strain_dup": dup[r["taxid"]] > 1 and r["taxid"] != lin.get("species_taxid"),
            "no_genus": not lin.get("genus_taxid"),
            "name_mismatch": bool(lin.get("species")) and bool(r["species_name"])
                             and lin["species"] != r["species_name"],
            "taxcheck_bad": r["taxcheck"] != "OK",
            "low_quality": (comp is None or comp < 90.0
                            or cont is None or cont > 5.0),
        }
        f["pass_gold"] = not (f["no_genus"] or f["name_mismatch"]
                              or f["taxcheck_bad"] or f["low_quality"])
        for k, v in f.items():
            if v:
                flags[k] += 1

        w.writerow([r["accession"], r["taxid"], r["organism"],
                    lin.get("domain", ""),
                    "Archaea" if r["accession"] in archaea else "Bacteria"]
                   + [lin.get(c, "") for rk in RANKS[1:] for c in (rk, rk + "_taxid")]
                   + [r["species_name"], r["paired_gcf"], r["level"], r["length"], r["n_contigs"],
                      r["completeness"], r["contamination"], r["taxcheck"]]
                   + [int(f[k]) for k in FLAG_ORDER])

print("wrote %s (%d assemblies)\n" % (out_path, len(rows)))

dom = collections.Counter((r["accession"] in archaea, lineage.get(r["taxid"], {}).get("domain"))
                          for r in rows)
print("domain (listed as archaeal, taxonkit domain):")
for k, v in sorted(dom.items(), key=lambda x: -x[1]):
    print("  %-28s %6d" % ("listed=%s / %s" % ("Archaea" if k[0] else "Bacteria", k[1] or "-"), v))

print("\nflags set (assemblies, not mutually exclusive):")
for k in FLAG_ORDER:
    print("  %-14s %6d  (%5.2f %%)" % (k, flags[k], 100.0 * flags[k] / len(rows)))
