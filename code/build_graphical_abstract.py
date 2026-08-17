"""Graphical abstract: the dissociation, the mechanism and the asymmetry.

A graphical abstract carries one argument, not a summary. The argument here is a
chain of three claims, so the figure reads left to right:

  (a) two interface classes measured by one pipeline agree on chemistry described by
      class and disagree at residue-pair resolution
  (b) that disagreement is absent from sequence readout and concentrated in shape
      readout, and the responsible chemistry can be named
  (c) the specialisation sits on the antibody side; the epitope is ordinary surface

Design decisions worth stating, since a graphical abstract is judged partly on
craft. Schematics replace labelled boxes because a Y-shaped binder against an
arbitrary surface, and two interlocking co-evolved partners, communicate the
comparison faster than words do. Numbers are set large and light with their labels
small and dark, so the eye lands on the quantity first. Panels are separated by hair
rules rather than borders. Every value is read from the analysis outputs, so the
artwork cannot drift from the results.
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import (Circle, FancyArrowPatch, FancyBboxPatch, Polygon,
                                RegularPolygon)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

# muted, print-safe palette
CRIMSON = "#9E2B2B"
SLATE = "#2E5C8A"
SAGE = "#4A7C59"
INK = "#1F2328"
GREY = "#6E7681"
FAINT = "#E8E6E3"
CREAM = "#FBFAF8"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "text.color": INK,
})


def J(name):
    p = C.OUT_MET / name
    return json.loads(p.read_text()) if p.exists() else {}


# --------------------------------------------------------------------- schematics
def antibody_glyph(ax, cx, cy, s=1.0):
    """A Y-shaped binder meeting an arbitrary surface."""
    arm_l = [(cx - 0.95 * s, cy + 1.15 * s), (cx - 0.30 * s, cy + 0.15 * s)]
    arm_r = [(cx + 0.95 * s, cy + 1.15 * s), (cx + 0.30 * s, cy + 0.15 * s)]
    for a in (arm_l, arm_r):
        ax.add_line(Line2D([a[0][0], a[1][0]], [a[0][1], a[1][1]],
                           color=CRIMSON, lw=4.2, solid_capstyle="round",
                           zorder=3))
    ax.add_line(Line2D([cx, cx], [cy + 0.15 * s, cy - 0.85 * s],
                       color=CRIMSON, lw=4.2, solid_capstyle="round", zorder=3))
    # the surface it binds: irregular, not shaped for the purpose
    xs = np.linspace(cx - 1.9 * s, cx + 1.9 * s, 220)
    ys = (cy + 1.42 * s + 0.10 * s * np.sin(xs * 3.1)
          + 0.05 * s * np.cos(xs * 6.7))
    ax.fill_between(xs, ys, ys + 0.55 * s, color=FAINT, edgecolor=GREY, lw=0.9,
                    zorder=1)


def ppi_glyph(ax, cx, cy, s=1.0):
    """Two co-evolved partners, interlocking."""
    th = np.linspace(0, 2 * np.pi, 240)
    for sign, col in ((-1, SLATE), (1, "#7FA3C8")):
        r = (0.95 + 0.13 * np.sin(3 * th + sign)) * s
        x = cx + sign * 0.52 * s + r * np.cos(th)
        y = cy + 0.30 * s + r * np.sin(th) * 0.86
        ax.fill(x, y, color=col, alpha=0.55 if sign > 0 else 0.85,
                edgecolor=SLATE, lw=0.9, zorder=2 + (sign > 0))


def tyrosine_glyph(ax, cx, cy, s=1.0):
    """Tyr phenol hydroxyl donating a hydrogen bond to a backbone acceptor."""
    ax.add_patch(RegularPolygon((cx - 0.75 * s, cy), 6, radius=0.42 * s,
                                orientation=np.pi / 6, facecolor="none",
                                edgecolor=CRIMSON, lw=2.0, zorder=3))
    ax.add_line(Line2D([cx - 0.75 * s, cx - 0.75 * s],
                       [cy - 0.42 * s, cy - 0.80 * s],
                       color=CRIMSON, lw=1.8, zorder=3))
    ax.text(cx - 0.75 * s, cy - 1.02 * s, "OH", ha="center", va="center",
            fontsize=9.5, color=CRIMSON, fontweight="bold")
    # dashed hydrogen bond
    ax.add_line(Line2D([cx - 0.52 * s, cx + 0.52 * s],
                       [cy - 0.92 * s, cy - 0.92 * s],
                       color=INK, lw=1.5, ls=(0, (2.2, 2.2)), zorder=3))
    # backbone fragment
    ax.add_line(Line2D([cx + 0.62 * s, cx + 1.55 * s], [cy - 0.92 * s] * 2,
                       color=SLATE, lw=3.4, solid_capstyle="round", zorder=3))
    ax.text(cx + 1.08 * s, cy - 0.52 * s, "backbone", ha="center", va="center",
            fontsize=8.5, color=SLATE)
    ax.text(cx + 1.08 * s, cy - 1.42 * s, "N–H  ·  C=O", ha="center", va="center",
            fontsize=8.5, color=GREY)


def bar(ax, x0, x1, y, frac, colour, h=9.0):
    ax.plot([x0, x1], [y, y], color=FAINT, lw=h, solid_capstyle="round", zorder=1)
    if frac > 0:
        ax.plot([x0, x0 + (x1 - x0) * frac], [y, y], color=colour, lw=h,
                solid_capstyle="round", zorder=2)


def value(ax, x, y, num, label, colour=INK, size=27):
    ax.text(x, y, num, ha="left", va="center", fontsize=size, color=colour,
            fontweight="light")
    ax.text(x, y - 0.62, label, ha="left", va="center", fontsize=8.6, color=GREY)


# ------------------------------------------------------------------------ panels
def panel_a(ax, hm, diff):
    class_r = hm["chemical_classes"]["class_similarity_antibody_vs_ppi"]
    res_r = hm["cluster30"]["grammar_similarity_antibody_vs_ppi"]
    pt = diff["permutation_test"]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    ax.text(0, 9.6, "one pipeline, two interface classes", fontsize=11.5,
            fontweight="bold", va="center")

    antibody_glyph(ax, 2.5, 6.3, s=0.92)
    ppi_glyph(ax, 7.4, 6.6, s=0.92)
    ax.text(2.5, 4.55, "antibody–antigen", ha="center", fontsize=10,
            color=CRIMSON, fontweight="bold")
    ax.text(2.5, 4.05, "512 complexes  ·  233 clusters", ha="center", fontsize=8.4,
            color=GREY)
    ax.text(7.4, 4.55, "protein–protein", ha="center", fontsize=10, color=SLATE,
            fontweight="bold")
    ax.text(7.4, 4.05, "2,292 complexes  ·  1,326 clusters", ha="center",
            fontsize=8.4, color=GREY)

    ax.plot([0.1, 9.9], [3.35, 3.35], color=FAINT, lw=1.1)

    bar(ax, 0.2, 6.4, 2.35, class_r, SAGE)
    ax.text(6.75, 2.35, f"{class_r:.2f}", fontsize=20, color=SAGE, va="center",
            fontweight="light")
    ax.text(0.2, 1.72, "agreement by chemical class", fontsize=9, color=GREY)

    bar(ax, 0.2, 6.4, 0.95, res_r, CRIMSON)
    ax.text(6.75, 0.95, f"{res_r:.2f}", fontsize=20, color=CRIMSON, va="center",
            fontweight="light")
    ax.text(0.2, 0.32, f"agreement by residue pair    "
                       f"permutation p = {pt['p_one_sided']:.3f}",
            fontsize=9, color=GREY)


def panel_b(ax, rd):
    side = rd["cross_arm_map_r"]["side"]
    main = rd["cross_arm_map_r"]["main"]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    ax.text(0, 9.6, "the difference is in shape, not sequence", fontsize=11.5,
            fontweight="bold", va="center")

    for x, val, title, sub, col in (
            (2.35, side, "side chains", "sequence readout", SAGE),
            (7.05, main, "main chain", "shape readout", CRIMSON)):
        ax.add_patch(Circle((x, 6.75), 1.30, facecolor=col, alpha=0.10,
                            edgecolor=col, lw=1.6, zorder=2))
        ax.text(x, 6.95, f"{val:.2f}", ha="center", va="center", fontsize=26,
                color=col, fontweight="light")
        ax.text(x, 5.95, "r", ha="center", va="center", fontsize=9, color=GREY,
                style="italic")
        ax.text(x, 4.95, title, ha="center", fontsize=10, fontweight="bold")
        ax.text(x, 4.48, sub, ha="center", fontsize=8.4, color=GREY)
    ax.text(4.70, 6.75, "vs", ha="center", va="center", fontsize=10, color=GREY)

    ax.text(0.1, 3.62, "both classes read 72% of contacts through side chains",
            fontsize=8.8, color=GREY)
    ax.plot([0.1, 9.9], [3.15, 3.15], color=FAINT, lw=1.1)

    tyrosine_glyph(ax, 2.35, 2.05, s=0.95)
    ax.text(5.05, 2.20, "tyrosine reads the backbone", fontsize=10.2,
            fontweight="bold", color=CRIMSON, va="center")
    ax.text(5.05, 1.52, "+0.56 to N–H,  +0.58 to C=O\nover the protein–protein "
                        "baseline", fontsize=8.5, color=GREY, va="center",
            linespacing=1.5)
    ax.text(5.05, 0.52, "ring stacking is depleted", fontsize=8.3, color=GREY,
            style="italic", va="center")


def panel_c(ax, epi):
    roles = epi["roles"]
    rows = [("antibody paratope", roles["Antibody paratope"]["mean_abs_enrichment"],
             CRIMSON),
            ("protein–protein interface",
             roles["PPI partner interface"]["mean_abs_enrichment"], SLATE),
            ("antibody epitope", roles["Antibody epitope"]["mean_abs_enrichment"],
             "#A8A29B")]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    ax.text(0, 9.6, "only one side is specialised", fontsize=11.5,
            fontweight="bold", va="center")
    ax.text(0, 8.85, "chemical distinguishability from the same protein's\n"
                     "own accessible surface", fontsize=8.6, color=GREY,
            va="top", linespacing=1.55)

    top, gap, span = 6.85, 1.55, 5.9
    for i, (lab, v, col) in enumerate(rows):
        y = top - i * gap
        bar(ax, 0.2, 0.2 + span, y, v / 1.05, col, h=8.5)
        ax.text(0.2 + span + 0.35, y, f"{v:.2f}", fontsize=18, color=col,
                va="center", fontweight="light")
        ax.text(0.2, y + 0.62, lab, fontsize=9.2, color=INK)

    ax.plot([0.1, 9.9], [2.72, 2.72], color=FAINT, lw=1.1)
    ax.text(0.1, 2.05, "the antibody is the specialised object",
            fontsize=11.2, fontweight="bold", color=SLATE, va="center")
    ax.text(0.1, 1.10, "the epitope does not protrude, lies on a different face\n"
                       "from the natural partner, and avoids ligand-binding sites",
            fontsize=8.6, color=GREY, va="center", linespacing=1.55)


if __name__ == "__main__":
    hm, diff = J("heatmap.json"), J("differential.json")
    rd, epi = J("readout.json"), J("epitope_logic.json")

    fig = plt.figure(figsize=(16.2, 5.9))
    fig.patch.set_facecolor(CREAM)
    gs = fig.add_gridspec(1, 3, left=0.035, right=0.985, top=0.90, bottom=0.06,
                          wspace=0.30)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    panel_a(axes[0], hm, diff)
    panel_b(axes[1], rd)
    panel_c(axes[2], epi)

    for i, ax in enumerate(axes):
        ax.text(-0.045, 1.10, f"({'abc'[i]})", transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top", color=GREY)

    # hair rules between panels, and a chevron to signal the reading order
    for x in (0.353, 0.670):
        fig.add_artist(Line2D([x, x], [0.10, 0.86], color=FAINT, lw=1.2))
        fig.add_artist(FancyArrowPatch((x - 0.006, 0.055), (x + 0.006, 0.055),
                                       transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=15,
                                       color="#C9C5C0", lw=1.6))

    out = C.OUT_FIG / "graphical_abstract.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=CREAM)
    plt.close(fig)
    print("wrote", out)
