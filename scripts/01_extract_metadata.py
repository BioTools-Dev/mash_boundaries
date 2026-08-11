#!/usr/bin/env python3
"""Flatten the NCBI datasets assembly report into one row per assembly.

Reads $META_JSONL and writes $DATA/genomes_raw.tsv. Only the fields the study
uses are kept:

  accession       GCA accession, the key everything else joins on
  taxid           taxid as submitted — a STRAIN taxid for most type material
  organism        organism name as submitted, strain suffix included
  species_name    NCBI's own parse of the species binomial
  species_taxid   species taxid resolved by the CheckM marker-set assignment
  taxcheck        NCBI's ANI-based check that the submitted name is consistent
  completeness    CheckM completeness (%)
  contamination   CheckM contamination (%)
  level           assembly level (Complete Genome / Chromosome / Scaffold / Contig)
  length          total sequence length (bp)
  n_contigs       number of contigs
  paired_gcf      the RefSeq accession of the same assembly, when it has one.
                  GTDB indexes a genome by its RefSeq accession whenever one
                  exists and by its GenBank accession otherwise, so without this
                  column most of the collection cannot be looked up in GTDB.

species_taxid and species_name come from two independent places in the report
and are kept side by side on purpose: disagreement between them flags a record
whose identity is not settled, and those are excluded from the gold standard
downstream rather than silently reconciled.
"""
import csv
import json
import os
import sys

src = os.environ.get("META_JSONL")
data = os.environ.get("DATA")
if not src or not data:
    sys.exit("ABORT: source config.sh first (META_JSONL and DATA must be set)")

out_path = os.path.join(data, "genomes_raw.tsv")
cols = ("accession taxid organism species_name species_taxid taxcheck "
        "completeness contamination level length n_contigs paired_gcf").split()

n = 0
with open(out_path, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for line in open(src):
        r = json.loads(line)
        ani = r.get("averageNucleotideIdentity") or {}
        cm = r.get("checkmInfo") or {}
        org = r.get("organism") or {}
        st = r.get("assemblyStats") or {}
        ai = r.get("assemblyInfo") or {}
        w.writerow([r.get("accession"), org.get("taxId"), org.get("organismName"),
                    ani.get("submittedSpecies"), cm.get("checkmSpeciesTaxId"),
                    ani.get("taxonomyCheckStatus"), cm.get("completeness"),
                    cm.get("contamination"), ai.get("assemblyLevel"),
                    st.get("totalSequenceLength"), st.get("numberOfContigs"),
                    (ai.get("pairedAssembly") or {}).get("accession")])
        n += 1

print("wrote %s (%d assemblies)" % (out_path, n))
