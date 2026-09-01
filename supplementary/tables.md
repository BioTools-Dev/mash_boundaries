# Supplementary tables

The catalogue of the supplementary tables: what each one is, and the files of `results/` and
`data/` it is assembled from. `scripts/33_supplementary.py` reads this file to build the
directory, so the names the paper uses and the files deposited here cannot disagree.

Tables S1 to S10 are the exhaustive form of what the paper reports in summary; each is what
the corresponding figure or table was computed from, so any curve can be rebuilt without
rerunning the comparison. Tables S11 and S12 are printed in full in the paper and deposited
here in their complete form. A table backed by a single file is one file; one backed by
several is both a directory of those files and a spreadsheet with one sheet per file.

- **Table S1.** The complete catalogue of taxonomic conflicts below distance 0.05: all 71
  genus pairs with their minimum distance and GTDB verdict, and all 199 individual queries
  with both accessions, both names and both GTDB genera
  in `supplementary/Table_S1/`, assembled by the pipeline from `results/conflicts_s10000_{genera,pairs}.tsv`.
- **Table S2.** The pairwise view in bins of 0.005 over the whole measured range, by class
  and by stratum, from which Table 2 and Figure 1a are aggregated
  in `supplementary/Table_S2.tsv`, assembled by the pipeline from `results/s10000_pairs.tsv`.
- **Table S3.** The per-genus macro estimates with their bootstrap intervals and quartiles,
  for the pairwise and the leave-one-out views, in `supplementary/Table_S3/`, assembled by the pipeline from `results/s10000_pairs_macro.tsv`,
  `results/macro_*.tsv`.
- **Table S4.** The leave-one-out curves at every cutoff for the three sketch sizes, the two
  taxonomies, both domains and each quality stratum, in the three scenarios
  in `supplementary/Table_S4/`, assembled by the pipeline from `results/loo_*[!s].tsv`. The
  per-query call tables the same step writes, `results/loo_*_calls.tsv`, are the input to
  steps 10, 11 and 17 rather than a curve, and are deposited with the other result tables.
- **Table S5.** The alignment-based ANI axis: band coverage, the measured distance–identity
  equivalence, and the leave-one-out curve at each alignment-fraction floor
  in `supplementary/Table_S5/`, assembled by the pipeline from `results/ani_*[!s].tsv`,
  `results/ani_*summary.txt`, `results/ani_conversion.txt` and
  `results/ani_s1000_conversion.txt`, the last two carrying the fitted correction to the
  customary conversion with the coefficient obtained at each sketch size and its largest
  departure from a band median.
- **Table S6.** The AAI validation: all 1,770 subset pairs with their reciprocal-best-hit
  count, their AAI and, where it falls inside the screening window, their protein sketch
  distance, in `supplementary/Table_S6.tsv`, assembled by the pipeline from `results/aai.tsv`.
- **Table S7.** The simulated-bin grid, of which the manuscript reports the clean series, with
  the per-bin calls behind every cell, in `supplementary/Table_S7/`, assembled by the pipeline from `results/bins_s10000{,_calls}.tsv`.
- **Table S8.** The gold-standard label table: one row per assembly with its resolved
  lineage, its rank taxids, its assembly statistics, the six disqualifying conditions kept as
  separate columns — `dup_taxid`, `strain_dup`, `no_genus`, `name_mismatch`, `taxcheck_bad`
  and `low_quality` — and the aggregate `pass_gold`, which is set when none of the last four
  is. The duplicate flags are deliberately left out of the aggregate (D9 of the repository):
  sharing a taxid disqualifies a pair, not a genome. The table is in `supplementary/Table_S8.tsv`, assembled by the pipeline from `data/genome_labels.tsv`.
- **Table S9.** The percentage of conserved proteins for every pair of both subsets, with the
  conserved and total protein counts of each genome, the Mash distance, and — for the first
  subset — the AAI and protein sketch distance of the same pair, so that the two genus criteria
  can be read side by side, in `supplementary/Table_S9/`, assembled by the pipeline from `results/pocp.tsv`,
  `results/pocp_window.tsv` and `results/pocp{,_window}_subset.txt`, the last of which gives
  the exact genome list of each subset.
- **Table S10.** The *Fusobacterium* test case: for each of the 533 genomes, its nearest type
  strain and the distance to it, the species this work calls at d ≤ 0.043, the name it carried
  in NCBI and the name assigned by the revision it is tested against
  in `supplementary/Table_S10.tsv`, assembled by the pipeline from `results/fusobacterium_calls.tsv`. The revised assignments are Supplementary Data 1 of
  Bi et al. (2026); the accessions of the 533 genomes are the first column of the same
  table.
- **Table S11.** The three sketch sizes at matched call volume, scenario `novel_species`,
  over all prokaryotes: for every target call volume and every sketch size, the cutoff whose
  call volume comes closest to the target, the volume actually reached and the precision of
  the genus call there. The paper prints the nine target volumes; the file carries every
  cutoff of the search. It is in `supplementary/Table_S11.tsv`, assembled by the pipeline
  from `results/three_sketch_sizes.tsv`.
- **Table S12.** The genus pairs in conflict below distance 0.05 and what the second taxonomy
  does with them, scenario `novel_species`, s = 10,000: for each pair, the number of
  conflicting queries, the closest pair the two genera contain, the GTDB verdict and the
  genus GTDB assigns. The paper prints the twelve leading pairs; the file carries all 71, in
  `supplementary/Table_S12.tsv`, assembled by the pipeline from
  `results/conflicts_s10000_genera.tsv`. The same rows, together with the 199 individual
  queries behind them, are Table S1.
