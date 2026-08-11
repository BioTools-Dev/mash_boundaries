#!/usr/bin/env python3
"""The three axes compared where they can be compared: at equal coverage.

Mash over DNA, alignment-based ANI and Mash over proteomes carry three different
scales, so their cutoffs are not comparable and a table indexed by cutoff would
be meaningless. What is comparable is the operating point: how precise each axis
is when it answers for the same number of genomes. That is the frame this table
is built on, and it is also the frame in which the interesting difference shows
up — the axes agree on precision while they all still answer, and differ in how
far they keep answering at all.

Each axis contributes its cumulative curve from the leave-one-out view under the
same scenario, and for every target coverage the table reports the cutoff that
reaches it and the precision there. An axis that cannot reach a coverage is left
empty rather than extrapolated: running out is a property of the index, not a
missing value.

The frame is not particular to the three axes: any set of indices measured over
the same genomes is comparable this way, so `--names` relabels the columns and
the same table serves the comparison of two values of k and of the three sketch
sizes.
"""
import argparse
import csv


def curve(path, scenario, extra=None):
    """(cutoff, calls, precision) triples, ascending in coverage."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["scenario"] != scenario:
                continue
            if extra and any(r[k] != v for k, v in extra.items()):
                continue
            cut = float(r["dist_hi"]) if "dist_hi" in r else float(r["ani_lo"])
            out.append((cut, int(r["cum_total"]), 100.0 * float(r["cum_p"])))
    out.sort(key=lambda t: t[1])
    return out


def at(c, target):
    """The point of a curve closest to a target coverage, or None if out of reach."""
    if not c or target > c[-1][1] * 1.02:
        return None
    return min(c, key=lambda t: abs(t[1] - target))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dna", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--ani", required=True)
    ap.add_argument("--ani-floor", default="15",
                    help="alignment-fraction floor of the third curve; empty when "
                         "the third curve is not an ANI table and carries no such "
                         "column")
    ap.add_argument("--names", default="Mash sobre DNA,Mash sobre proteomas,"
                                       "ANI por alineamiento",
                    help="column headings, in the order --dna,--protein,--ani")
    ap.add_argument("--scenario", default="novel_species")
    ap.add_argument("--out", required=True)
    ap.add_argument("--targets", default="4000,8000,12000,15000,18500,21000,23500,26500,29000")
    a = ap.parse_args()

    names = [s.strip() for s in a.names.split(",")]
    floor = {"af_floor": a.ani_floor} if a.ani_floor else None
    axes = list(zip(names, [curve(a.dna, a.scenario),
                            curve(a.protein, a.scenario),
                            curve(a.ani, a.scenario, floor)]))
    targets = [int(x) for x in a.targets.split(",")]

    with open(a.out + ".tsv", "w") as fh:
        fh.write("target_calls\taxis\tcutoff\tcalls\tprecision\n")
        for t in targets:
            for name, c in axes:
                p = at(c, t)
                if p:
                    fh.write("%d\t%s\t%.4f\t%d\t%.4f\n" % (t, name, p[0], p[1], p[2]))
                else:
                    fh.write("%d\t%s\t\t\t\n" % (t, name))

    out = []
    out.append("%s, a igual cobertura, escenario %s"
               % (" frente a ".join(n for n, _ in axes), a.scenario))
    out.append("")
    out.append("%-10s %s" % ("llamadas", "  ".join("%-28s" % n for n, _ in axes)))
    for t in targets:
        cells = []
        for name, c in axes:
            p = at(c, t)
            cells.append("corte %6.3f  %6d  %6.2f%%" % (p[0], p[1], p[2]) if p
                         else "%-28s" % "  — (el eje se acaba)")
        out.append("%-10d %s" % (t, "  ".join(cells)))
    out.append("")
    for name, c in axes:
        if c:
            out.append("%-22s alcance máximo: %d llamadas, %.2f %% de precisión"
                       % (name, c[-1][1], c[-1][2]))
    text = "\n".join(out)
    print(text)
    with open(a.out + "_summary.txt", "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
