# Supplementary tables

The supplementary tables of the manuscript, under the names the paper
uses. This directory is assembled from `results/` and `data/` by
`scripts/33_supplementary.py`, which reads the manuscript itself, so the
two views cannot disagree. Each entry lists the files it was built from,
which is where the pipeline writes them.

## Table S1

The complete catalogue of taxonomic conflicts below distance 0.05: all 71

**Here:** `Table_S1/ (2 files) and Table_S1.xlsx`

**Built from:** `results/conflicts_s10000_genera.tsv`, `results/conflicts_s10000_pairs.tsv`

## Table S2

The pairwise view in bins of 0.005 over the whole measured range, by class

**Here:** `Table_S2.tsv`

**Built from:** `results/s10000_pairs.tsv`

## Table S3

The per-genus macro estimates with their bootstrap intervals and quartiles,

**Here:** `Table_S3/ (7 files) and Table_S3.xlsx`

**Built from:** `results/macro_prot_k7.tsv`, `results/macro_prot_k7_gtdb.tsv`, `results/macro_s100000.tsv`, `results/macro_s10000_gtdb.tsv`, `results/macro_s10000_novel.tsv`, `results/macro_s1000_novel.tsv`, `results/s10000_pairs_macro.tsv`

## Table S4

The leave-one-out curves at every cutoff for the three sketch sizes, the two

**Here:** `Table_S4/ (17 files) and Table_S4.xlsx`

**Built from:** `results/loo_prot_k7.tsv`, `results/loo_prot_k7_archaea.tsv`, `results/loo_prot_k7_bacteria.tsv`, `results/loo_prot_k7_bacteria_gold.tsv`, `results/loo_prot_k7_gtdb.tsv`, `results/loo_prot_k9.tsv`, `results/loo_s10000.tsv`, `results/loo_s100000.tsv`, `results/loo_s100000_gtdb.tsv`, `results/loo_s10000_archaea.tsv`, `results/loo_s10000_archaea_gold.tsv`, `results/loo_s10000_bacteria.tsv`, `results/loo_s10000_bacteria_gold.tsv`, `results/loo_s10000_bacteria_quality.tsv`, `results/loo_s10000_bacteria_taxcheck.tsv`, `results/loo_s10000_gtdb.tsv`, `results/loo_s1000_d0.28.tsv`

## Table S5

The alignment-based ANI axis: band coverage, the measured distance–identity

**Here:** `Table_S5/ (18 files) and Table_S5.xlsx`

**Built from:** `results/ani_conversion.txt`, `results/ani_coverage.tsv`, `results/ani_dist_by_ani.tsv`, `results/ani_gtdb_coverage.tsv`, `results/ani_gtdb_dist_by_ani.tsv`, `results/ani_gtdb_loo.tsv`, `results/ani_gtdb_summary.txt`, `results/ani_loo.tsv`, `results/ani_s1000_conversion.txt`, `results/ani_s1000_coverage.tsv`, `results/ani_s1000_dist_by_ani.tsv`, `results/ani_s1000_loo.tsv`, `results/ani_s1000_summary.txt`, `results/ani_summary.txt`, `results/ani_vs_prot_coverage.tsv`, `results/ani_vs_prot_dist_by_ani.tsv`, `results/ani_vs_prot_loo.tsv`, `results/ani_vs_prot_summary.txt`

## Table S6

The AAI validation: all 1,770 subset pairs with their reciprocal-best-hit

**Here:** `Table_S6.tsv`

**Built from:** `results/aai.tsv`

## Table S7

The simulated-bin grid, of which the manuscript reports the clean series, with

**Here:** `Table_S7/ (2 files) and Table_S7.xlsx`

**Built from:** `results/bins_s10000.tsv`, `results/bins_s10000_calls.tsv`

## Table S8

The gold-standard label table: one row per assembly with its resolved

**Here:** `Table_S8.tsv`

**Built from:** `data/genome_labels.tsv`

## Table S9

The percentage of conserved proteins for every pair of both subsets, with the

**Here:** `Table_S9/ (4 files) and Table_S9.xlsx`

**Built from:** `results/pocp.tsv`, `results/pocp_subset.txt`, `results/pocp_window.tsv`, `results/pocp_window_subset.txt`

## Table S10

The *Fusobacterium* test case: for each of the 533 genomes, its nearest type

**Here:** `Table_S10.tsv`

**Built from:** `results/fusobacterium_calls.tsv`

## Table S11

The three sketch sizes at matched call volume. Scenario `novel_species`, all

**Here:** `Table_S11.tsv`

**Built from:** `results/three_sketch_sizes.tsv`

## Table S12

Genus pairs in conflict below distance 0.05, and what the second taxonomy does
with them. Scenario `novel_species`, s = 10,000. The twelve leading pairs, ordered by number

**Here:** `Table_S12.tsv`

**Built from:** `results/conflicts_s10000_genera.tsv`

