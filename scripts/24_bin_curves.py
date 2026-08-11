#!/usr/bin/env python3
"""What a metagenome bin sees, and whether containment gives it back.

Reads a `mash dist` stream of simulated bins against the reference collection
and answers three questions per cell of the completeness x contamination grid:

  1. does a fixed distance cutoff still work,
  2. does the nearest neighbour still carry the right genus even when the cutoff
     has stopped working,
  3. does containment recover what the Jaccard index loses.

The three are different questions and they are predicted to have different
answers. A bin holding the fraction c of its genome's k-mers has Jaccard index
exactly c against that genome, so its distance to *itself* is fixed by c alone —
0.019 at half completeness, 0.052 at a fifth — and a fixed cutoff is spent on
incompleteness before it is spent on taxonomy. But that penalty applies to every
candidate at once, so the ranking may survive intact while the cutoff fails.

Containment is the third answer. C(bin in reference) = |B and R| / |B| does not
depend on how much of the genome the bin holds, only on how much of what it does
hold is explained by the reference, and it is therefore the quantity a bin should
have been queried with all along. It is not reported by `mash dist`, but it
follows from what is:

    C = J (|B| + |R|) / ((1 + J) |B|)

with J the Jaccard estimate behind the reported distance. Sizes come from the
assembly length of the reference and from the bases the generator recorded for
the bin, so this is a size-corrected containment rather than a direct estimate;
`--screen-check` validates it against `mash screen`, which measures containment
without the correction. Note that the |B| factor cancels out of any comparison
within one bin, so the *ranking* by containment does not depend on the bin size
estimate at all — only the fixed-threshold columns do.

The scenario is `novel_species` throughout, as in the leave-one-out view: every
genome of the query's own species is removed from the database, including the
genome the bin was cut from. Anything else would let a bin find its own source
and would measure nothing.

A query whose genus has no representative left in the database after the
exclusions cannot be called correctly, and counting it as an error would measure
the composition of the sample instead of the effect of the bin. Every cell is
therefore reported over its *answerable* bins, with the total alongside so the
reader sees where the restriction bites. It bites in exactly one place: a
congeneric contaminant whose species is also removed takes a second species out
of the query's own genus, which empties a two-species genus outright. Without
the restriction that cell reads 62 % where the measurement is 92.7 %.

Contamination is reported under two conditions, because the first measurement
made on this grid showed they are not the same experiment. A contaminant does
not merely dilute the query: if the organism it came from is itself in the
reference, it stands there as an exact match competing with the query's own
congenerics, which are by definition not exact. So the grid is read once with
the donor's species present in the database — the pessimistic case, and the real
one whenever the contaminating organism happens to be a described species — and
once with the donor's species removed as well, which isolates dilution from
decoying. Reporting only the first would blame incompleteness for a failure that
belongs to the reference; reporting only the second would understate what a bin
actually meets.
"""
import argparse
import collections
import csv
import math
import re
import sys

ACC = re.compile(r"(GC[AF]_\d+\.\d+)")
KEEP = 40           # neighbours retained per bin, enough to survive the exclusion


def jaccard_of(d, k):
    """The Jaccard index behind a Mash distance."""
    if d <= 0:
        return 1.0
    x = math.exp(-k * d)
    return x / (2.0 - x)


def containment_cut(d, k):
    """The containment threshold matching a distance cutoff between equals.

    Two genomes of the same size at Mash distance d share a fraction
    2J/(1+J) of their k-mers, and that fraction is what a bin cut from one of
    them would show as containment in the other, whatever its completeness.
    """
    j = jaccard_of(d, k)
    return 2.0 * j / (1.0 + j)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", default="-", help="mash dist output, or - for stdin")
    ap.add_argument("--bins", required=True, help="<prefix>_manifest.tsv from step 23")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--kmer", type=int, default=21)
    ap.add_argument("--genus-cut", type=float, default=0.13)
    ap.add_argument("--species-cut", type=float, default=0.043)
    a = ap.parse_args()

    lab = {}
    species_of_genus = collections.defaultdict(set)
    with open(a.labels) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lab[r["accession"]] = r
            if r["low_quality"] == "0":
                species_of_genus[r["genus"]].add(r["species_taxid"])

    bins = {}
    with open(a.bins) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["size"] = int(r["native_bases"]) + int(r["foreign_bases"])
            bins[r["bin"] + ".fna"] = r

    # --- one pass over the stream, keeping the nearest few per bin ------------
    near = collections.defaultdict(list)
    fh = sys.stdin if a.dist == "-" else open(a.dist)
    n_lines = 0
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        n_lines += 1
        b = bins.get(f[1])
        if b is None:
            continue
        m = ACC.search(f[0])
        if not m:
            continue
        d = float(f[2])
        v = near[f[1]]
        v.append((d, m.group(1)))
        if len(v) > 4 * KEEP:                 # trimmed rarely, kept cheap
            v.sort()
            del v[KEEP:]
    if a.dist != "-":
        fh.close()
    print("aristas leídas: %d; bins con al menos un vecino: %d"
          % (n_lines, len(near)), file=sys.stderr)

    gcut_c = containment_cut(a.genus_cut, a.kmer)
    scut_c = containment_cut(a.species_cut, a.kmer)
    print("umbrales de contención equivalentes: género %.4f, especie %.4f"
          % (gcut_c, scut_c), file=sys.stderr)

    # --- per bin: the two rankings, under novel_species -----------------------
    cells = collections.defaultdict(lambda: collections.Counter())
    calls = []
    for name, v in near.items():
        b = bins[name]
        sp = b["species_taxid"]
        donor_sp = ""
        if b["donor"] and b["donor"] in lab:
            donor_sp = lab[b["donor"]]["species_taxid"]
        cand = []
        for d, acc in sorted(set(v)):
            r = lab.get(acc)
            if r is None or r["species_taxid"] == sp:
                continue                      # novel_species: own species removed
            j = jaccard_of(d, a.kmer)
            size_r = float(r["length"])
            cont = j * (b["size"] + size_r) / ((1.0 + j) * b["size"])
            cand.append((d, min(cont, 1.0), acc, r))

        for donor_state in ("present", "absent"):
            pool = cand
            if donor_state == "absent" and donor_sp:
                pool = [t for t in cand if t[3]["species_taxid"] != donor_sp]
            key = (b["completeness"], b["contamination"], b["donor_class"],
                   donor_state)
            c = cells[key]
            c["bins"] += 1

            # Could this bin have been called correctly at all? Only if its own
            # genus still holds a species after the exclusions.
            gone = {sp}
            if donor_state == "absent" and donor_sp and \
                    lab[b["donor"]]["genus"] == b["genus"]:
                gone.add(donor_sp)
            if not (species_of_genus[b["genus"]] - gone):
                c["unanswerable"] += 1
                continue
            c["answerable"] += 1
            if not pool:
                c["no_neighbour"] += 1
                continue

            by_d = min(pool, key=lambda t: t[0])
            by_c = max(pool, key=lambda t: t[1])
            ok_d = by_d[3]["genus"] == b["genus"]
            ok_c = by_c[3]["genus"] == b["genus"]
            for tag, ok in (("nn_d", ok_d), ("nn_c", ok_c)):
                c[tag + "_calls"] += 1
                c[tag + "_right"] += 1 if ok else 0
            if by_d[0] <= a.genus_cut:
                c["cut_d_calls"] += 1
                c["cut_d_right"] += 1 if ok_d else 0
            if by_c[1] >= gcut_c:
                c["cut_c_calls"] += 1
                c["cut_c_right"] += 1 if ok_c else 0
            calls.append((name, b["source"], b["completeness"], b["contamination"],
                          b["donor_class"], donor_state, "%.6f" % by_d[0], by_d[2],
                          ok_d, "%.6f" % by_c[1], by_c[2], ok_c))

    def pct(c, tag):
        n = c[tag + "_calls"]
        return 100.0 * c[tag + "_right"] / n if n else float("nan")

    order = sorted(cells, key=lambda t: (t[3] != "present", -float(t[0]),
                                         float(t[1]), t[2]))
    with open(a.out + ".tsv", "w") as fh:
        fh.write("completeness\tcontamination\tdonor\tdonor_species\tbins\t"
                 "answerable\tno_neighbour\tnn_dist_calls\tnn_dist_p\t"
                 "nn_cont_calls\tnn_cont_p\tcut_dist_calls\tcut_dist_p\t"
                 "cut_cont_calls\tcut_cont_p\n")
        for key in order:
            c = cells[key]
            fh.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t"
                     "%d\t%.4f\n"
                     % (key[0], key[1], key[2], key[3], c["bins"],
                        c["answerable"], c["no_neighbour"],
                        c["nn_d_calls"], pct(c, "nn_d"),
                        c["nn_c_calls"], pct(c, "nn_c"),
                        c["cut_d_calls"], pct(c, "cut_d"),
                        c["cut_c_calls"], pct(c, "cut_c")))

    with open(a.out + "_calls.tsv", "w") as fh:
        fh.write("bin\tsource\tcompleteness\tcontamination\tdonor\tdonor_species\t"
                 "nn_dist\tnn_dist_acc\tnn_dist_ok\tnn_cont\tnn_cont_acc\t"
                 "nn_cont_ok\n")
        for r in calls:
            fh.write("\t".join(str(x) for x in r) + "\n")

    out = ["Lo que ve un bin, escenario novel_species",
           "umbral de género: distancia <= %.3f  |  contención >= %.4f"
           % (a.genus_cut, gcut_c), "",
           "las celdas se calculan sobre los bins medibles: aquellos cuyo género "
           "conserva\nalguna especie en la base tras las exclusiones", "",
           "%-7s %-6s %-6s %-8s %6s %6s   %-18s %-18s   %-18s %-18s"
           % ("compl", "cont", "donante", "sp.don.", "bins", "medib.",
              "vecino por dist.", "vecino por cont.",
              "corte de dist.", "corte de cont.")]
    for key in order:
        c = cells[key]
        out.append("%-7s %-6s %-6s %-8s %6d %6d   %6d %10.2f%%  %6d %10.2f%%   "
                   "%6d %10.2f%%  %6d %10.2f%%"
                   % (key[0], key[1], key[2], key[3], c["bins"], c["answerable"],
                      c["nn_d_calls"], pct(c, "nn_d"),
                      c["nn_c_calls"], pct(c, "nn_c"),
                      c["cut_d_calls"], pct(c, "cut_d"),
                      c["cut_c_calls"], pct(c, "cut_c")))
    text = "\n".join(out)
    print(text)
    with open(a.out + "_summary.txt", "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
