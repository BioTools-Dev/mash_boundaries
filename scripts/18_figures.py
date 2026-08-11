#!/usr/bin/env python3
"""The study's figures, from the result tables written by the earlier steps.

Nothing here recomputes anything: every panel reads a TSV from results/, so a
figure cannot disagree with the number it illustrates. Each figure is written as
SVG, PDF and PNG at 300 dpi, with the vector as the master.

Design decisions that apply to all four, and the reason for each:

  three colours, in fixed order   blue, orange, aqua. Assigned to entities, not
                                  to rank, so the same series keeps its hue from
                                  one figure to the next. The set is checked for
                                  colour-vision separation rather than eyeballed.
  no second y axis                a panel carries one scale. Precision and call
                                  count are different quantities and get their
                                  own panels; aligning two scales on one frame
                                  invents a relation the data does not have.
  legend plus selective labels    identity is never carried by colour alone, and
                                  no value is printed on every point.
  recessive chrome                hairline grid on one axis, no top or right
                                  spine, text in ink rather than in the series
                                  colour.
"""
import argparse
import collections
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8981"
GRID = "#e8e7e3"
AXIS = "#c9c8c3"
BAND = "#eeedea"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8.5,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.8,
    "lines.solid_capstyle": "round",
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "legend.labelcolor": INK2,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    # Fixed salt: without it matplotlib draws the SVG element ids at random and
    # every regeneration rewrites the whole file. With it the figures are
    # byte-identical run to run.
    "svg.hashsalt": "mash_boundaries",
})

PCT = FuncFormatter(lambda v, _: "%g %%" % v)


def read(path, **where):
    """Rows of a TSV, optionally filtered by exact column values."""
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if all(r[k] == v for k, v in where.items()):
                yield r


def panel(ax, title, subtitle=None, xlabel=None, ylabel=None, pct=True, head=0.0):
    """Frame a panel. `head` reserves room above it for a bracket row."""
    ax.set_title(title, loc="left", color=INK, fontsize=9.5,
                 pad=(16 if subtitle else 8) + head * 260)
    if subtitle:
        ax.text(0, 1.02 + head, subtitle, transform=ax.transAxes, color=INK2,
                fontsize=8, va="bottom")
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, length=3, width=0.8, labelsize=8)
    if xlabel:
        ax.set_xlabel(xlabel, color=INK2, fontsize=8.5)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK2, fontsize=8.5)
    if pct:
        ax.yaxis.set_major_formatter(PCT)


def window(ax, lo, hi, label=None, y=0.955):
    """Shade a proposed distance window behind the marks."""
    ax.axvspan(lo, hi, color=BAND, lw=0, zorder=0)
    if label:
        ax.text((lo + hi) / 2, y, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", color=MUTED, fontsize=7.5)


def bracket(ax, lo, hi, label, y=1.045):
    """A labelled span above the frame, for a proposed window."""
    ax.annotate("", (lo, y), (hi, y), xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="|-|,widthA=0.35,widthB=0.35",
                                color=MUTED, lw=0.9, shrinkA=0, shrinkB=0))
    ax.text((lo + hi) / 2, y + 0.02, label, transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", color=INK2, fontsize=7.5)


def tag(ax, x, y, text, color, dx=6, dy=0):
    """A direct label: a colour dot carries identity, the text stays ink."""
    ax.plot([x], [y], marker="o", ms=4.5, color=color, mec=SURFACE, mew=1.2,
            zorder=5, clip_on=False)
    ax.annotate(text, (x, y), textcoords="offset points", xytext=(dx, dy),
                color=INK2, fontsize=7.5, va="center", zorder=5)


def save(fig, outdir, name):
    """The vector formats are the masters; the PNG is for reading drafts.

    The creation date is suppressed in every format so that regenerating an
    unchanged figure produces a byte-identical file: otherwise the whole figure
    set shows up as modified on every run and a real change becomes invisible
    among them.
    """
    stamp = {"svg": {"Date": None},
             "pdf": {"CreationDate": None},
             "png": {"Software": None}}
    for ext in ("svg", "pdf", "png"):
        fig.savefig(os.path.join(outdir, "%s.%s" % (name, ext)),
                    dpi=300, bbox_inches="tight", metadata=stamp[ext])
    plt.close(fig)
    print("  %s.{svg,pdf,png}" % name)


# --------------------------------------------------------------------------- 1
def fig_pairs(res, out):
    """The pairwise view: both curves, and what weighting by genus does to them."""
    rows = [r for r in read(os.path.join(res, "s10000_pairs.tsv"), stratum="all")
            if float(r["dist_hi"]) <= 0.30]
    x = [(float(r["dist_lo"]) + float(r["dist_hi"])) / 2 for r in rows]
    sp = [100 * float(r["p_same_species"]) for r in rows]
    gn = [100 * float(r["p_same_genus"]) for r in rows]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.4))

    panel(a1, "a  Any pair, at distance d",
          "P(class | d), all reported pairs",
          "Mash distance", "probability")
    a1.plot(x, gn, color=S1, label="same genus")
    a1.plot(x, sp, color=S2, label="same species")
    a1.axhline(50, color=AXIS, lw=0.8, zorder=1)
    a1.set_ylim(-3, 103)
    a1.set_xlim(0, 0.30)
    a1.text(0.297, 51, "one half", ha="right", va="bottom", color=MUTED, fontsize=7.5)
    tag(a1, 0.0325, 40.2, "0.030", S2, dx=7, dy=6)
    tag(a1, 0.2075, 43.1, "0.205", S1, dx=7, dy=6)
    a1.legend(loc="center left", bbox_to_anchor=(0.16, 0.62), fontsize=8)

    # Bins resting on a handful of genera swing wildly and say nothing; the
    # threshold is declared rather than smoothed away.
    mac = [r for r in read(os.path.join(res, "s10000_pairs_macro.tsv"),
                           domain="Bacteria")
           if float(r["dist_lo"]) <= 0.295 and int(r["genera"]) >= 5]
    mx = [float(r["dist_lo"]) + 0.0025 for r in mac]
    panel(a2, "b  The weight of the large genera",
          "P(same genus | d) in bacteria, bins with ≥ 5 genera",
          "Mash distance", None)
    a2.fill_between(mx, [100 * float(r["ci_lo"]) for r in mac],
                    [100 * float(r["ci_hi"]) for r in mac],
                    color=S2, alpha=0.16, lw=0)
    a2.plot(mx, [100 * float(r["micro_p"]) for r in mac], color=S1,
            label="unweighted (micro)")
    a2.plot(mx, [100 * float(r["macro_p"]) for r in mac], color=S2,
            label="weighted by genus (macro)")
    a2.set_ylim(-3, 103)
    a2.set_xlim(0, 0.30)
    a2.legend(loc="lower left", fontsize=8)
    a2.annotate("95 % CI of the macro", (0.168, 45), textcoords="offset points",
                xytext=(-96, 26), color=INK2, fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))

    fig.tight_layout()
    save(fig, out, "fig1_vista_de_pares")


# --------------------------------------------------------------------------- 2
def fig_loo(res, out):
    """The leave-one-out view, across sketch sizes and across taxonomies."""
    # The two large sketches are drawn one over the other on purpose: they
    # coincide, and that they coincide is the result of §3.6.
    series = [("s=1,000, NCBI", "loo_s1000_d0.28.tsv", S1),
              ("s=10,000, NCBI", "loo_s10000.tsv", S2),
              ("s=100,000, NCBI", "loo_s100000.tsv", S4),
              ("s=10,000, GTDB", "loo_s10000_gtdb.tsv", S3)]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.4))
    panel(a1, "a  Precision of the genus call",
          "nearest neighbour, scenario novel_species",
          "Mash distance cutoff", "cumulative precision")
    panel(a2, "b  How many genomes receive a call",
          "the same cutoffs, the same query", "Mash distance cutoff",
          "genomes with a call", pct=False)
    for ax in (a1, a2):
        window(ax, 0.043, 0.13, "proposed window")

    for label, fname, color in series:
        # A cumulative precision resting on a few hundred queries is noise, so
        # each curve starts where it has a thousand calls behind it.
        rows = [r for r in read(os.path.join(res, fname), scenario="novel_species")
                if float(r["dist_hi"]) <= 0.25 and int(r["cum_total"]) >= 1000]
        x = [float(r["dist_hi"]) for r in rows]
        a1.plot(x, [100 * float(r["cum_p"]) for r in rows], color=color, label=label)
        a2.plot(x, [int(r["cum_total"]) for r in rows], color=color, label=label)

    a1.set_ylim(88, 100.6)
    a1.set_xlim(0, 0.25)
    # Two columns, in the corridor between the NCBI plateau and the GTDB curve:
    # with four series there is no free vertical block left in this panel.
    a1.legend(loc="center", bbox_to_anchor=(0.54, 0.79), ncol=2, fontsize=8,
              columnspacing=1.2, handlelength=1.6)
    a2.set_xlim(0, 0.25)
    a2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d k" % (v / 1000)))
    a2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.88), fontsize=8)
    a1.axvline(0.13, color=MUTED, lw=0.8, zorder=1)
    a1.text(0.133, 88.3, "0.13", color=MUTED, fontsize=7.5, va="bottom")

    fig.tight_layout()
    save(fig, out, "fig2_vecino_mas_cercano")


# --------------------------------------------------------------------------- 3
def fig_ani(res, out):
    """The ANI axis: where it puts the same pairs, and where it stops existing."""
    rows = [r for r in read(os.path.join(res, "ani_coverage.tsv"))
            if r["class"] == "all" and float(r["dist_hi"]) <= 0.25]
    x = [(float(r["dist_lo"]) + float(r["dist_hi"])) / 2 for r in rows]
    med = [float(r["ani_median"]) for r in rows]
    q1 = [float(r["ani_q1"]) for r in rows]
    q3 = [float(r["ani_q3"]) for r in rows]
    cov = [100 * float(r["coverage"]) for r in rows]
    # solid while the axis resolves most of the band, faint while it resolves a
    # minority, and absent once the median rests on a handful of pairs
    k = max([i for i, c in enumerate(cov) if c >= 50] or [0]) + 1
    kf = max([i for i, c in enumerate(cov) if c >= 1] or [0]) + 1

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.4))
    panel(a1, "a  The two axes on the same pairs",
          "alignment-based ANI against Mash distance",
          "Mash distance", "ANI (%)")
    a1.plot([0, 0.25], [100, 75], color=MUTED, lw=0.9, ls=(0, (5, 4)), zorder=2,
            label="the customary conversion, 1 − d")
    # the corrected conversion, drawn only over the range it was fitted on
    a1.plot([0, 0.13], [100, 100 - 100 * 1.12 * 0.13], color=S2, lw=1.1,
            ls=(0, (1.6, 1.8)), zorder=4,
            label="the corrected conversion, 1 − 1.12 d")
    a1.fill_between(x[:k], q1[:k], q3[:k], color=S1, alpha=0.18, lw=0)
    a1.plot(x[:k], med[:k], color=S1, zorder=3)
    a1.plot(x[k - 1:kf], med[k - 1:kf], color=S1, lw=1.1, alpha=0.4, zorder=3)
    a1.set_xlim(0, 0.25)
    a1.set_ylim(76, 101)
    a1.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    a1.legend(loc="lower left", bbox_to_anchor=(0.0, 0.02), fontsize=7.5)
    tag(a1, 0.0446, 95.0, "95 % ANI at d = 0.045", S1, dx=8, dy=9)
    a1.annotate("coverage < 50 %,\nmedian is biased", (0.196, 81.6),
                textcoords="offset points", xytext=(16, 44), color=MUTED,
                fontsize=7.5, ha="right",
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))

    panel(a2, "b  How far the ANI axis reaches",
          "pairs with an ANI, of those Mash reports in each band",
          "Mash distance", "coverage")
    a2.fill_between(x, cov, color=S1, alpha=0.16, lw=0)
    a2.plot(x, cov, color=S1)
    a2.set_xlim(0, 0.25)
    a2.set_ylim(-3, 103)
    a2.axhline(50, color=AXIS, lw=0.8, zorder=1)
    tag(a2, 0.1975, 7.6, "beyond 0.20 the axis goes silent", S1, dx=-152, dy=14)

    fig.tight_layout()
    save(fig, out, "fig4_eje_de_ani")


# --------------------------------------------------------------------------- 4
def fig_cutoffs(res, out):
    """The corollary: what each candidate cutoff costs and buys."""
    sp = list(read(os.path.join(res, "species_s1000.tsv")))
    sp = [r for r in sp if float(r["cutoff"]) <= 0.15]
    sx = [float(r["cutoff"]) for r in sp]
    gn = [r for r in read(os.path.join(res, "loo_s1000_d0.28.tsv"),
                          scenario="novel_species") if float(r["dist_hi"]) <= 0.15]

    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    panel(ax, "The two cutoffs, and what each one costs",
          "s=1,000, bacteria, NCBI labels",
          "Mash distance cutoff", "precision or recall", head=0.085)
    ax.plot(sx, [100 * float(r["precision"]) for r in sp], color=S2,
            label="precision of the species call")
    ax.plot(sx, [100 * float(r["recall"]) for r in sp], color=S3,
            label="recall of the species call")
    ax.plot([float(r["dist_hi"]) for r in gn],
            [100 * float(r["cum_p"]) for r in gn], color=S1,
            label="precision of the genus call")
    ax.set_xlim(0, 0.15)
    ax.set_ylim(60, 102)
    bracket(ax, 0.001, 0.043, "species", y=1.03)
    bracket(ax, 0.046, 0.13, "genus", y=1.03)
    for d, lbl in ((0.043, "0.043"), (0.13, "0.13")):
        ax.axvline(d, color=MUTED, lw=0.8, zorder=1)
        ax.text(d + 0.0015, 60.8, lbl, color=MUTED, fontsize=7.5, va="bottom")
    ax.axvline(0.05, color=AXIS, lw=0.8, zorder=1)
    ax.annotate("0.05, the customary cutoff", (0.05, 88), textcoords="offset points",
                xytext=(8, 0), color=MUTED, fontsize=7.5, va="center")
    ax.legend(loc="lower left", bbox_to_anchor=(0.005, 0.02), fontsize=8)

    fig.tight_layout()
    save(fig, out, "fig6_corolario_de_los_cortes")


# --------------------------------------------------------------------------- 5
def fig_strata(res, out):
    """The same curve read on one domain at a time and on one quality gate at a time."""
    def curve(fname):
        rows = [r for r in read(os.path.join(res, fname), scenario="novel_species")
                if float(r["dist_hi"]) <= 0.25 and int(r["cum_total"]) >= 100]
        return ([float(r["dist_hi"]) for r in rows],
                [100 * float(r["cum_p"]) for r in rows])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.4))

    panel(a1, "a  By domain",
          "cumulative precision of the genus call, s=10,000",
          "Mash distance cutoff", "cumulative precision")
    for fname, color, label in (("loo_s10000_bacteria.tsv", S1, "Bacteria"),
                                ("loo_s10000_archaea.tsv", S2, "Archaea")):
        x, y = curve(fname)
        a1.plot(x, y, color=color, label=label)
    a1.set_xlim(0, 0.25)
    a1.set_ylim(88, 101.9)
    window(a1, 0.043, 0.13, "proposed window")
    a1.legend(loc="upper right", bbox_to_anchor=(1.0, 0.98), fontsize=8)

    panel(a2, "b  By quality stratum",
          "bacteria only; each gate on its own",
          "Mash distance cutoff", None)
    for fname, color, label in (
            ("loo_s10000_bacteria.tsv", S1, "unfiltered"),
            ("loo_s10000_bacteria_quality.tsv", S2, "CheckM: completeness and contamination"),
            ("loo_s10000_bacteria_taxcheck.tsv", S3, "NCBI ANI-based check"),
            ("loo_s10000_bacteria_gold.tsv", S4, "both (gold)")):
        x, y = curve(fname)
        a2.plot(x, y, color=color, label=label)
    a2.set_xlim(0, 0.25)
    a2.set_ylim(88, 101.9)
    window(a2, 0.043, 0.13, "proposed window")
    a2.legend(loc="upper right", bbox_to_anchor=(1.0, 0.98), fontsize=7.5)

    fig.tight_layout()
    save(fig, out, "fig3_dominio_y_calidad")


# --------------------------------------------------------------------------- 6
def fig_axes(res, out):
    """The three axes side by side, and the validation of the third one."""
    def operating(fname, key="cum_total"):
        rows = [r for r in read(os.path.join(res, fname), scenario="novel_species")
                if int(r["cum_total"]) >= 1000]
        return ([int(r[key]) for r in rows], [100 * float(r["cum_p"]) for r in rows])

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12.0, 3.6))

    # (a) precision against coverage: the only frame on which axes with different
    # scales can be compared at all
    panel(a1, "a  The three axes, at matched coverage",
          "precision against number of calls, novel_species",
          "genomes with a call", "cumulative precision")
    for fname, color, label in (("loo_s10000.tsv", S1, "Mash over DNA, s=10,000"),
                                ("loo_prot_k7.tsv", S2, "Mash over proteomes, k=7")):
        x, y = operating(fname)
        a1.plot(x, y, color=color, label=label)
    ani = [r for r in read(os.path.join(res, "ani_loo.tsv"),
                           scenario="novel_species", af_floor="15")
           if int(r["cum_total"]) >= 1000]
    a1.plot([int(r["cum_total"]) for r in ani],
            [100 * float(r["cum_p"]) for r in ani], color=S3,
            label="alignment-based ANI (af ≥ 15 %)")
    a1.set_ylim(85, 98)
    a1.set_xlim(0, 30000)
    a1.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d k" % (v / 1000)))
    a1.legend(loc="lower left", fontsize=8)
    end = max(ani, key=lambda r: int(r["cum_total"]))
    tag(a1, int(end["cum_total"]), 100 * float(end["cum_p"]),
        "the ANI axis runs out here", S3, dx=-24, dy=-56)

    # (b) the validation: a real AAI against the sketch that stands in for it
    rows = [r for r in read(os.path.join(res, "aai.tsv")) if r["prot_dist"]]
    panel(a2, "b  The proteome axis is AAI, measured",
          "AAI by reciprocal best hits, %d pairs" % len(rows),
          "Mash distance over proteomes", "AAI (%)", pct=False)
    for cl, color, label in (("congeneric", S1, "congeneric"),
                             ("inter_genus", S2, "inter-genus")):
        sel = [r for r in rows if r["class"] == cl]
        a2.plot([float(r["prot_dist"]) for r in sel],
                [float(r["aai"]) for r in sel], "o", ms=4, color=color,
                mec=SURFACE, mew=0.6, alpha=0.85, label=label, ls="none")
    sel = [r for r in rows if r["class"] == "same_species"]
    if sel:
        a2.plot([float(r["prot_dist"]) for r in sel],
                [float(r["aai"]) for r in sel], "o", ms=4, color=S3,
                mec=SURFACE, mew=0.6, label="conspecific", ls="none")
    a2.axhline(65, color=MUTED, lw=0.9)
    a2.text(0.012, 65.9, "genus criterion: 65 % AAI", color=MUTED, fontsize=7.5)
    a2.set_xlim(0, 0.5)
    a2.set_ylim(38, 100)
    a2.legend(loc="upper right", fontsize=8)

    # (c) the other formal criterion, on the axis the cutoffs are stated in. POCP
    # counts shared genes rather than measuring how similar they are, so it is
    # the one check of the genus cutoff that does not reduce to an identity.
    # Both subsets: the one drawn for the AAI validation reaches the far range
    # where the criterion is crossed, and the one drawn by distance band fills
    # the genus window, which the first leaves almost empty. Pairs with no
    # shared hash have no distance and are absent rather than placed at 1.
    XMAX = 0.32
    rows = [r for f in ("pocp.tsv", "pocp_window.tsv")
            for r in read(os.path.join(res, f))
            if r["mash_dist"] and float(r["mash_dist"]) <= XMAX]
    panel(a3, "c  The other genus criterion, POCP",
          "conserved proteins against Mash distance, %d pairs" % len(rows),
          "Mash distance over DNA", "POCP (%)", pct=False)
    for cl, color, label in (("same_species", S3, "conspecific"),
                             ("congeneric", S1, "congeneric"),
                             ("inter_genus", S2, "inter-genus")):
        sel = [r for r in rows if r["class"] == cl]
        if sel:
            a3.plot([float(r["mash_dist"]) for r in sel],
                    [float(r["pocp"]) for r in sel], "o", ms=4, color=color,
                    mec=SURFACE, mew=0.6, alpha=0.85, label=label, ls="none")
    a3.axhline(50, color=MUTED, lw=0.9)
    a3.axvline(0.13, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    a3.set_xlim(0, XMAX)
    a3.set_ylim(0, 100)
    a3.text(0.006, 51.6, "genus criterion: 50 % POCP", color=MUTED, fontsize=7.5)
    a3.text(0.136, 3, "the genus cutoff, 0.13", color=MUTED, fontsize=7.5)
    a3.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    save(fig, out, "fig5_tres_ejes")


# --------------------------------------------------------------------------- 7
def fig_bins(res, out):
    """What a metagenome bin sees, and what containment gives back."""
    path = os.path.join(res, "bins_s10000.tsv")
    if not os.path.exists(path):
        return
    rows = list(read(path))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.6))

    # (a) the clean bin: the cutoff fails while the ranking survives
    clean = sorted((r for r in rows
                    if r["contamination"] == "0.00" and r["donor_species"] == "present"),
                   key=lambda r: float(r["completeness"]))
    x = [100 * float(r["completeness"]) for r in clean]
    panel(a1, "a  A clean bin, as completeness falls",
          "genus call, scenario novel_species",
          "bin completeness", "precision")
    for key, color, label in (("nn_dist_p", S1, "nearest neighbour by distance"),
                              ("nn_cont_p", S3, "nearest neighbour by containment"),
                              ("cut_dist_p", S2, "fixed distance cutoff, 0.13"),
                              ("cut_cont_p", S4, "fixed containment threshold")):
        a1.plot(x, [float(r[key]) for r in clean], color=color, label=label,
                marker="o", ms=3.5, mec=SURFACE, mew=0.8)
    a1.set_xlim(0, 105)
    a1.legend(loc="lower right", fontsize=7.5)

    # (b) how many bins the fixed cutoff still answers for: the failure of a
    # cutoff shows up as silence, not as error, so it needs its own panel
    panel(a2, "b  How many bins each fixed threshold answers for",
          "the same bins, the same threshold", "bin completeness",
          "bins with a call", pct=False)
    for key, color, label in (("cut_dist_calls", S2, "distance cutoff, 0.13"),
                              ("cut_cont_calls", S4, "containment threshold"),
                              ("bins", MUTED, "bins in the cell")):
        a2.plot(x, [int(r[key]) for r in clean], color=color, label=label,
                marker="o", ms=3.5, mec=SURFACE, mew=0.8)
    a2.set_xlim(0, 105)
    a2.set_ylim(bottom=0)
    a2.legend(loc="lower right", fontsize=7.5)

    fig.tight_layout()
    save(fig, out, "fig7_bins_de_metagenoma")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.environ.get("RESULTS", "results"))
    ap.add_argument("--out", default=os.environ.get("FIGURES", "figures"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print("figuras en %s:" % a.out)
    # in the order the manuscript cites them
    fig_pairs(a.results, a.out)
    fig_loo(a.results, a.out)
    fig_strata(a.results, a.out)
    fig_ani(a.results, a.out)
    fig_axes(a.results, a.out)
    fig_cutoffs(a.results, a.out)
    fig_bins(a.results, a.out)


if __name__ == "__main__":
    main()
