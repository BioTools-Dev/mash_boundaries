# mash_boundaries

Where the species and the genus boundaries fall in Mash distance space, measured over every
pair of the prokaryotic type material: **30,209 genomes, 4.56 × 10⁸ pairs**.

The measured result is a pair of boundaries — **species at d ≤ 0.043, genus at
0.043 < d ≤ 0.13, and abstention beyond** — together with three corrections for anyone using
Mash thresholds today: identity derived as 1 − d overstates it by 12 d points and should be
read as **ANI = 1 − 1.12 d**; a distance threshold is uninterpretable without the sketch size
that produced it; and ten thousand hashes exhaust the estimator, so the error left above 0.15
is real overlap between genera rather than something a larger sketch can buy away.

The manuscript is in `manuscript/mash_boundaries_paper.md`. **Section numbers of the form
§n.n cited from the code refer to that manuscript.** This file documents how to reproduce
every number and every figure in it.

## Layout

```
mash_boundaries/
├── README.md              this file
├── LICENSE                MIT for the code, CC-BY-4.0 for data, results and figures
├── config.sh.example      path template; config.sh is local and is not versioned
├── manuscript/            the paper, its figures' captions and its working documents
├── data/                  gold-standard labels and manifests (versioned)
├── scripts/               all the code, numbered in dependency order
├── results/               every output table the paper is computed from (versioned)
├── figures/               the seven figures, in SVG + PDF + PNG at 300 dpi (versioned)
├── supplementary/         the twelve supplementary tables, under the names the paper uses
├── fusobacterium/         the test case's genomes and intermediates — not versioned
├── proteome/              RefSeq proteomes — regenerable, not versioned
├── sketch/                Mash sketches — regenerable, not versioned
├── dist/                  pairwise distances — regenerable, not versioned
└── logs/                  run logs, not versioned
```

`results/` and `figures/` are versioned in full, so **every claim in the paper can be checked
without rerunning anything**. The regenerable intermediates are not: the raw pairwise
distances alone are hundreds of gigabytes, and each is rebuilt by the numbered script that
produced it, from the manifests that are versioned.

## Requirements

| Tool | Version used | What for |
|---|---|---|
| Mash | 2.3 | sketching and all-versus-all distances |
| skani | 0.1.4 | alignment-based ANI, the independent axis |
| DIAMOND | 2.2.4 | AAI by reciprocal best hits, and POCP |
| TaxonKit | 0.20.0 | taxid → lineage, with a pinned NCBI dump |
| NCBI `datasets` | — | fetching proteomes |
| zstd | 1.5.7 | compressing the pairwise streams |
| Python | 3.14, with matplotlib for the figures | analysis and figures |
| pandoc | — | only for the `.docx` of the manuscript |

No Python package beyond matplotlib is needed: the analysis uses the standard library.

The whole study runs on a single machine. The costly steps, timed on 16 threads on the runs
actually performed rather than extrapolated from the benchmark of step 06, are the
all-versus-all — 6 min 34 s at s = 1,000, 1 h 01 m at s = 10,000 and 11 h 13 m at
s = 100,000 — the protein all-versus-all at 12 min 50 s, the skani triangle at 52 min 04 s,
and the bin sweep at about 1 hour. **Every conclusion of
the paper is reachable at s = 10,000.** Peak disk use is dominated by the skani index (24 GB)
and the s = 100,000 pair stream.

## Configuration

No script hardcodes a path: they all come from `config.sh`, which is local and is not
versioned. What is shared is the template.

```bash
cp config.sh.example config.sh    # then edit the local paths at the top
source ./config.sh
```

The paths that must be set are the genome directory, the NCBI assembly report, the pinned
taxonomy dump and the working directory. `THREADS` defaults to 16.

## Input data

The 30,209 assemblies are the prokaryotic **type material** of NCBI, which is what makes a
gold standard possible: a type strain is the name-bearing reference of its species. They are
not redistributed here — the accessions are in `data/genome_manifest.tsv` and the assemblies
are fetched from NCBI under their terms.

The genus and species of a genome are taken from its **taxid** through TaxonKit, never from
its name, which breaks on `Candidatus` and on names in brackets. The taxonomy dump is pinned
to a date, and GTDB r232 is carried alongside as a second taxonomy so that every curve can be
read under both.

Two of the versioned tables carry a `path` column pointing at the machine they were built on:
`data/genome_manifest.tsv` and `data/proteome_manifest.tsv`. The durable content of both is
the accession; the paths are rewritten for the local machine by step 05 and step 19
respectively, and nothing downstream depends on the values committed here.

## Reproducing, step by step

The scripts are numbered in dependency order. Steps 01–04 are deterministic and cheap.

```bash
python3 scripts/01_extract_metadata.py   # NCBI jsonl        -> data/genomes_raw.tsv
bash    scripts/02_taxonkit_lineage.sh   # taxids            -> data/lineage.tsv
python3 scripts/03_build_labels.py       # join + flags      -> data/genome_labels.tsv
python3 scripts/04_pair_census.py        # the base rate     -> results/pair_census.txt
bash    scripts/05_build_sketches.sh     # manifest + the three sketches
```

The all-versus-all takes the sketch sizes as arguments and the screening cutoff as a variable.
The cutoff differs by sketch size because the small sketch cannot express a distance above
0.296:

```bash
THREADS=16 bash scripts/08_allvsall.sh 1000            # common window, d <= 0.28
CUTOFF=0.40 bash scripts/08_allvsall.sh 10000 100000   # past the small sketch's ceiling
```

Then the two views of the same measurement — pairwise, and leave-one-out nearest neighbour —
with the genus-weighted curve, the second taxonomy and the conflict catalogue. Every step takes
its inputs explicitly, so the sketch size and the scenario are always visible in the command:

```bash
# the leave-one-out curve and its per-query call table, for each sketch size
for s in 1000:0.28 10000:0.40 100000:0.40; do
    S=${s%%:*}; D=${s##*:}
    python3 scripts/09_leave_one_out.py \
        --nn dist/s${S}_d${D}_nn.tsv \
        --calls results/loo_s${S}_calls.tsv \
        --out results/loo_s${S}.tsv
done

# the same distances relabelled with GTDB, and the per-stratum curves
python3 scripts/12_gtdb_labels.py --labels data/genome_labels.tsv \
    --gtdb-dir data/gtdb --out data/gtdb_labels.tsv
python3 scripts/09_leave_one_out.py --nn dist/s10000_d0.40_nn.tsv \
    --relabel data/gtdb_keys.tsv --calls results/loo_s10000_gtdb_calls.tsv \
    --out results/loo_s10000_gtdb.tsv
python3 scripts/09_leave_one_out.py --nn dist/s10000_d0.40_nn.tsv \
    --domain Bacteria --out results/loo_s10000_bacteria.tsv        # and Archaea
python3 scripts/09_leave_one_out.py --nn dist/s10000_d0.40_nn.tsv \
    --labels data/genome_labels.tsv --domain Bacteria --stratum gold \
    --out results/loo_s10000_bacteria_gold.tsv                     # and quality, taxcheck

# genus-weighted curve, conflict catalogue, pairwise view, species cutoff
python3 scripts/10_macro_curve.py --calls results/loo_s10000_calls.tsv \
    --domain Bacteria --out results/macro_s10000_novel.tsv
python3 scripts/11_taxonomic_conflicts.py --calls results/loo_s10000_calls.tsv \
    --labels data/genome_labels.tsv --gtdb data/gtdb_keys.tsv \
    --out results/conflicts_s10000
python3 scripts/16_pair_curves.py --hist dist/s10000_d0.40_hist.tsv \
    --hist-genus dist/s10000_d0.40_hist_genus.tsv \
    --labels data/genome_labels.tsv --out results/s10000
python3 scripts/17_species_cutoff.py --calls results/loo_s1000_calls.tsv \
    --labels data/genome_labels.tsv --out results/species_s1000
```

The alignment-based ANI axis is independent of the Mash all-versus-all and can be run in any
order; step 15 is what crosses them:

```bash
bash scripts/14_skani_allvsall.sh                      # 2 min sketch + 52 min triangle
python3 scripts/15_ani_axis.py \
    --skani dist/skani_ani.tsv.zst \
    --mash-pairs dist/s10000_d0.40_pairs.tsv.zst \
    --labels data/genome_labels.tsv \
    --mash-calls results/loo_s10000_calls.tsv \
    --out results/ani                                  # add --relabel data/gtdb_keys.tsv for GTDB
```

Step 15 also fits the corrected distance-to-identity conversion and writes it to
`results/ani_conversion.txt`. It walks the 100.8 M edges in a little over a minute.

The protein axis needs the proteomes first:

```bash
bash scripts/19_fetch_proteomes.sh
bash scripts/20_build_protein_sketches.sh

python3 scripts/21_aai_validation.py \
    --labels data/genome_labels.tsv --faa proteome/faa \
    --mash-pairs dist/prot_k7_d0.50_pairs.tsv.zst --out results/aai

python3 scripts/22_three_axes.py \
    --dna results/loo_s10000.tsv --protein results/loo_prot_k7.tsv \
    --ani results/ani_loo.tsv --ani-floor 15 --out results/three_axes.tsv
```

The percentage of conserved proteins, the second formal genus criterion, is step 26. It is run
**twice with the same script**, because the two runs answer different questions: the first
reuses the AAI subset so that every pair carries both criteria, the second draws its own,
stratified by distance band, to populate the genus window that the first leaves nearly empty.

```bash
python3 scripts/26_pocp_validation.py \
    --labels data/genome_labels.tsv --manifest data/genome_manifest.tsv \
    --faa proteome/faa --sketch-dir sketch \
    --diamond $DIAMOND --mash $MASH \
    --subset results/aai_subset.txt --aai results/aai.tsv \
    --out results/pocp                                 # the 60 genomes of the AAI validation

python3 scripts/26_pocp_validation.py \
    --labels data/genome_labels.tsv --manifest data/genome_manifest.tsv \
    --faa proteome/faa --sketch-dir sketch \
    --diamond $DIAMOND --mash $MASH \
    --stratified-from dist/s10000_d0.40_pairs.tsv.zst --aai results/aai.tsv \
    --out results/pocp_window                          # 80 genomes inside the window
```

The subset actually used is written to `<prefix>_subset.txt`, and the distances come from a
sketch of that subset rather than from the all-versus-all, because most of its inter-genus
pairs fall outside the screening window and would have no distance to read. Pairs sharing no
hash at all are reported as having **no distance** rather than as distance 1.

The *Fusobacterium* test case starts from the article's Europe PMC identifier — its
supplementary data carries the revised species of each genome, which is the gold standard the
cutoffs are tested against — and needs no file placed by hand:

```bash
python3 scripts/27_fuso_supplementary.py --work $FUSO/data/supp \
    --labels-out data/fusobacterium_bi2026_labels.tsv \
    --out $FUSO/data/bi2026_anib.tsv

$DATASETS download genome accession --inputfile data/fusobacterium_accessions.txt \
    --include genome --filename $FUSO/data/genomes.zip     # 533 of 540; 7 are suppressed
bash scripts/28_fuso_axes.sh                               # Mash and skani over the same genomes

python3 scripts/29_fuso_bridge.py --anib $FUSO/data/bi2026_anib.tsv \
    --mash s10000=$FUSO/data/mash_s10000.raw s100000=$FUSO/data/mash_s100000.raw \
    --skani $FUSO/data/skani_ani.tsv --out $FUSO/results/bridge
python3 scripts/30_fuso_depth.py --pairs $FUSO/results/bridge_pairs.tsv \
    --column anib --identity --out results/fusobacterium_depth_anib.tsv
python3 scripts/31_fuso_calls.py --labels data/genome_labels.tsv \
    --bi data/fusobacterium_bi2026_labels.tsv --dist $FUSO/data/query_vs_type.raw \
    --out results/fusobacterium_calls.tsv
```

Seven of the 540 genomes the article used have since been suppressed from RefSeq and cannot be
retrieved from either RefSeq or GenBank; the analysis runs on the remaining 533 and none of the
seven defines an edge of the gap. The revised assignments themselves are versioned here, so the
comparison can be checked without re-downloading anything.

Finally, what an incomplete metagenome bin sees:

```bash
python3 scripts/23_make_bins.py \
    --labels data/genome_labels.tsv --manifest data/genome_manifest.tsv \
    --out bins/bins_s10000 --work bins/work --mash $MASH --queries 400

$MASH dist -p $THREADS -d 0.40 sketch/type_s10000.msh bins/bins_s10000.msh \
  | python3 scripts/24_bin_curves.py \
        --bins bins/bins_s10000_manifest.tsv --labels data/genome_labels.tsv \
        --out results/bins_s10000
```

The sweep streams the distances rather than storing them: at 400 queries over the completeness
and contamination grid it is the ~1 hour step. The genus and species cutoffs are parameters
(`--genus-cut`, `--species-cut`) and default to the ones this work proposes.

### What each script writes

| Script | Output |
|---|---|
| `01_extract_metadata.py` | `data/genomes_raw.tsv` — one row per assembly |
| `02_taxonkit_lineage.sh` | `data/lineage.tsv` — standard-rank lineage per taxid |
| `03_build_labels.py` | `data/genome_labels.tsv` — **the gold standard**, with quality flags |
| `04_pair_census.py` | `results/pair_census.txt` — counts by class and stratum |
| `05_build_sketches.sh` | `data/genome_manifest.tsv`, `sketch/type_s{1000,10000,100000}.msh` |
| `06_triangle_benchmark.sh` | `logs/06_benchmark.log` — measured and extrapolated cost |
| `07_accumulate.py` | consumes the edge stream; invoked by step 08 |
| `08_allvsall.sh` | `dist/s<S>_d<cutoff>_{hist,hist_genus,nn,summary,pairs}` |
| `09_leave_one_out.py` | `results/loo_s<S>_*.tsv` — leave-one-out curves and call tables |
| `10_macro_curve.py` | `results/macro_s<S>_*.tsv` — genus-weighted curve with bootstrap |
| `11_taxonomic_conflicts.py` | `results/conflicts_s<S>_{pairs,genera,worst}.tsv` |
| `12_gtdb_labels.py` | `data/gtdb_labels.tsv` — GTDB lineage per genome, and coverage |
| `13_skani_benchmark.sh` | `logs/13_skani_benchmark.log` — cost, and the effect of `--min-af` |
| `14_skani_allvsall.sh` | `dist/skani_ani.tsv.zst` — ANI for every resolvable pair |
| `15_ani_axis.py` | `results/ani_{coverage,dist_by_ani,conversion,loo,loo_calls,summary}` |
| `16_pair_curves.py` | `results/s<S>_pairs{,_macro,_macro_bands,_summary}` |
| `17_species_cutoff.py` | `results/species_s<S>[_gtdb]{.tsv,_summary.txt}` |
| `18_figures.py` | `figures/fig<N>_*.{svg,pdf,png}` |
| `19_fetch_proteomes.sh` | `proteome/faa/<accession>.faa.gz`, `data/proteome_manifest.tsv` |
| `20_build_protein_sketches.sh` | `sketch/prot_k<K>_s<S>.msh` |
| `21_aai_validation.py` | `results/aai.tsv` — AAI by reciprocal best hits |
| `22_three_axes.py` | `results/{three_axes,prot_k7_vs_k9,three_sketch_sizes}.tsv` |
| `23_make_bins.py` | `bins/` — simulated bins with controlled completeness |
| `24_bin_curves.py` | `results/bins_s<S>{,_calls}.tsv` |
| `25_manuscript_docx.sh` | the `.docx` of any Markdown document of `manuscript/` |
| `26_pocp_validation.py` | `results/pocp{,_window}{.tsv,_summary.txt,_subset.txt}` |
| `27_fuso_supplementary.py` | `data/fusobacterium_bi2026_labels.tsv` and the ANI matrix in long form |
| `28_fuso_axes.sh` | Mash at two sketch sizes and skani over the same 533 genomes |
| `29_fuso_bridge.py` | the three scales joined pair by pair, and the conversion refitted |
| `30_fuso_depth.py` | `results/fusobacterium_depth_{anib,mash}.tsv` — the gap against sampling depth |
| `31_fuso_calls.py` | `results/fusobacterium_calls.tsv` — Table S10 |
| `33_supplementary.py` | `supplementary/` — the ten tables under the names the paper uses |
| `xlsx.py` | a spreadsheet reader, so the supplementary needs no extra dependency |

## The figures

The seven figures are built from the tables in `results/`, never by recomputing: a figure
cannot disagree with the number it illustrates.

```bash
$PYVIZ scripts/18_figures.py      # PYVIZ is the interpreter with matplotlib, see config.sh
```

Each is written as **SVG, PDF and PNG at 300 dpi**, with the vector as the master and the text
live rather than outlined, so it stays editable. The run is **byte-reproducible**: the SVG hash
salt is fixed and creation dates are suppressed in all three formats, so regenerating an
unchanged figure produces an identical file and a real change is visible in a diff.

| Figure | What it shows |
|---|---|
| `fig1_vista_de_pares` | P(same species \| d) and P(same genus \| d), unweighted and weighted by genus |
| `fig2_vecino_mas_cercano` | leave-one-out precision and call volume, three sketch sizes and two taxonomies |
| `fig3_dominio_y_calidad` | the same curve by domain and by quality gate |
| `fig4_eje_de_ani` | Mash distance against alignment-based ANI, both conversions, and how far the ANI axis reaches |
| `fig5_tres_ejes` | the three axes at matched coverage, AAI against the sketch that stands in for it, and POCP against Mash distance |
| `fig6_corolario_de_los_cortes` | precision and coverage of the species and genus calls, with the two proposed windows |
| `fig7_bins_de_metagenoma` | precision and call volume as completeness falls |

## The manuscript

```bash
bash scripts/25_manuscript_docx.sh          # the paper; TABLE_PT sets the table font size
```

The Markdown is the master and the `.docx` is derived, so a correction is never made twice.
The conversion embeds the 300 dpi PNGs and fits the wide tables to the page; a journal should
receive the vector masters in `figures/*.svg` and `figures/*.pdf` alongside.

## Supplementary tables

The twelve supplementary tables of the manuscript are in `supplementary/`, under the names
the paper uses:

```
supplementary/
├── README.md                      what each table is, and the files it was built from
├── Table_S1/  + Table_S1.xlsx     taxonomic conflicts below 0.05
├── Table_S2.tsv                   the pairwise view in bins of 0.005
├── Table_S3/  + Table_S3.xlsx     per-genus macro estimates, with bootstrap intervals
├── Table_S4/  + Table_S4.xlsx     the leave-one-out curves, every sketch size and scenario
├── Table_S5/  + Table_S5.xlsx     the alignment-based ANI axis and the fitted conversion
├── Table_S6.tsv                   the AAI validation, pair by pair
├── Table_S7/  + Table_S7.xlsx     the simulated-bin grid and its per-bin calls
├── Table_S8.tsv                   the gold-standard label table
├── Table_S9/  + Table_S9.xlsx     the percentage of conserved proteins, both subsets
├── Table_S10.tsv                  the Fusobacterium test case, genome by genome
├── Table_S11.tsv                  the three sketch sizes at matched call volume
└── Table_S12.tsv                  the genus pairs in conflict, all 71
```

A table backed by one file is one file. One backed by several is both a directory of those
files and a spreadsheet with one sheet per file, which is the form a journal expects: the
files do not share a schema, so merging them into a single flat table would invent one, and
separate sheets keep each as it was written.

The directory is assembled by `scripts/33_supplementary.py`, which reads the manuscript's own
supplementary section, so the two cannot disagree. The files keep the names the pipeline gave
them, so each traces back to the step that wrote it, and the spreadsheets are written without
any dependency beyond the standard library.

## Design decisions referenced from the code

The code cites these by number.

| | Decision, and why |
|---|---|
| **D1** | The main question is answered on real pairs, not on simulated reads. Simulation would measure the simulator. |
| **D2** | Type material as the gold standard: a type strain is the name-bearing reference of its species, so its label is as good as a label gets. |
| **D3** | Archaea are separated, not discarded. They are 1,014 genomes, 3.4 % of the collection, and reporting their curve as its own stratum costs nothing: if the boundaries differ between domains, that is a finding. |
| **D4** | Genus and species come from the taxid, never from the name. |
| **D5** | Three sketch sizes, so that the answer's dependence on the sketch is measured rather than assumed. |
| **D6** | Alignment-based ANI as an independent axis: skani screens with its own k-mer filter, so neither metric decides what the other sees. |
| **D7** | A third axis on amino acid identity, which is the index taxonomy actually uses to delimit the genus. |
| **D8** | Two views of the same measurement, and they differ: pairwise P(class \| d), and leave-one-out nearest neighbour. Both are reported for every curve. |
| **D9** | Filtering is an analysis decision, not a label-construction one. The label table discards nothing: every disqualifying condition is a column with its counts, so strata are defined downstream and sensitivities are measured rather than assumed. |
| **D10** | Sketches are rebuilt from a canonical manifest, so the identifiers cannot drift from the labels. |
| **D11** | The all-versus-all is screened at the source and what is not reported is recovered by subtraction, because retaining 4.56 × 10⁸ raw pairs is not viable. |
| **D12** | The AAI axis is built on RefSeq proteomes and at k = 7, chosen on a measurement rather than inherited from the conventional protein k of 9. |
| **D13** | Bins are simulated from the type material itself and measured against their own intact genome, so the correct answer is known exactly. |
| **D14** | The corrected conversion is a single coefficient fitted through the origin: identity is 100 % at distance zero by construction, so the line has no intercept to spend. |
| **D16** | The genus cutoff is the far edge of the precision plateau, located on the 0.005 grid the curve was computed on rather than on the coarser grid it was first reported at. Within the plateau coverage is free, so the cutoff goes at its wide end. |
| **D15** | The second genus criterion is measured on the same pairs as the first, so that POCP and AAI can be compared genome by genome, plus a second subset drawn by distance band to populate the genus window. |

## Data availability

The genomes are the NCBI prokaryotic type-material collection, retrieved with the `datasets`
client under the accessions in `data/genome_manifest.tsv`. The taxonomy is a pinned NCBI dump
and GTDB r232. The gold-standard label table, the pair census and every result table
underlying the figures and tables are in this repository. The large, regenerable intermediates
— sketches, raw pairwise distances and the DIAMOND working directories — are not deposited;
each is rebuilt by the numbered script that produced it.

## Licence

MIT for everything in `scripts/`; CC-BY-4.0 for `data/`, `results/` and `figures/`. See
`LICENSE`. The genome assemblies themselves are not redistributed.

## Citation

The manuscript is not yet published. Until it is, cite this repository; the accession and DOI
will be added here once the archived version exists.
