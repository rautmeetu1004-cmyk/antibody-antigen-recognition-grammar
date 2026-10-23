"""Recreate Fig. 11 (comparator decomposition) from the committed result JSONs.

Panel A: cross-class recognition-grammar correlation within biologically defined subclasses
of the protein-protein reference set, each with its 95% bootstrap interval and a size-matched
null band. Panel B: the antibody-vs-comparator effect (rank-biserial) for apolar-apolar
burial, within the same subclasses. Reads comparator_grammar.json and comparator_physics.json;
writes figures/comparator_combined.png. (Replaces an earlier ad-hoc figure with no generator;
legend repositioned so it no longer overlaps the bars, rendered at 400 dpi.)
"""
from __future__ import annotations
import json, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = pathlib.Path(__file__).resolve().parent.parent
G = json.load(open(ROOT / "metadata/analysis/comparator_grammar.json"))["results"]
P = json.load(open(ROOT / "metadata/analysis/comparator_physics.json"))["results"]["f_apolar_apolar"]["subclasses"]

# (json key, label, colour)
ROWS = [
    ("all_ppi", "full protein-protein", "#000000"),
    ("heteromeric", "heteromeric", "#2166ac"),
    ("no_disulfide", "no disulfide", "#2166ac"),
    ("homomeric", "homomeric", "#2166ac"),
    ("disulfide", "disulfide-bonded", "#762a83"),
    ("enzyme_inhibitor", "enzyme-inhibitor", "#b2182b"),
    ("receptor_ligand", "receptor-ligand", "#b2182b"),
]
GROUPS = [Patch(fc="#2166ac", label="clean co-evolved"),
          Patch(fc="#762a83", label="hydrophobic/obligate"),
          Patch(fc="#b2182b", label="transient/polar")]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.4))
y = np.arange(len(ROWS))[::-1]

# ---- Panel A: grammar correlation ----
for yi, (k, lab, col) in zip(y, ROWS):
    g = G[k]; r = g["r_obs"]; lo, hi = g["boot_ci"]; nl, nm, nh = g["size_matched_null"]
    if k != "all_ppi":
        axA.barh(yi, nh - nl, left=nl, height=0.62, color="0.85", zorder=0)   # null band
        axA.plot([nm, nm], [yi - 0.31, yi + 0.31], color="0.5", lw=1, zorder=1)
    axA.errorbar(r, yi, xerr=[[r - lo], [hi - r]], fmt="o", color=col, ms=7, capsize=3, zorder=3)
    axA.text(r, yi + 0.33, f"{r:.2f}", ha="center", va="bottom", fontsize=9, color=col, fontweight="bold")
    axA.text(-0.02, yi, f"{lab}\n(u={g['n_units']})", ha="right", va="center", fontsize=8.5, color=col)
axA.axvline(0.42, ls="--", color="0.4", lw=1); axA.text(0.42, len(ROWS)-0.3, " full-PPI 0.42", fontsize=7.5, color="0.4")
axA.axvline(0.80, ls=":", color="#2166ac", lw=1.4); axA.text(0.80, len(ROWS)-1.2, " one-grammar\n null ≈ 0.80", fontsize=8, color="#2166ac")
axA.set_xlim(-0.55, 0.95); axA.set_ylim(-0.6, len(ROWS)-0.3); axA.set_yticks([])
axA.set_xlabel("cross-class grammar correlation  r\n(lower = more divergent)", fontsize=9)
axA.set_title("A   Recognition grammar", fontsize=11, fontweight="bold", loc="left")
axA.legend(handles=GROUPS, fontsize=8, frameon=True, loc="lower right")   # legend in the empty lower-right of A

# ---- Panel B: apolar-burial effect ----
for yi, (k, lab, col) in zip(y, ROWS):
    rb = P.get(k, {}).get("rank_biserial")
    if rb is None:
        continue
    axB.barh(yi, rb, height=0.62, color=col, alpha=0.85)
    axB.text(rb - 0.015, yi + 0.33, f"{rb:.2f}", ha="right", va="bottom", fontsize=9, color=col, fontweight="bold")
axB.axvline(0, color="0.3", lw=1)
axB.set_ylim(-0.6, len(ROWS)-0.3); axB.set_yticks([])
axB.set_xlabel("antibody vs comparator effect\n(rank-biserial, apolar-apolar burial)", fontsize=9)
axB.set_title("B   Interface physics", fontsize=11, fontweight="bold", loc="left")

fig.suptitle("Decomposing the protein-protein reference set: physics stays stable, grammar divergence is sharpened not explained",
             fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(ROOT / "figures/comparator_combined.png", dpi=400, bbox_inches="tight")
print("wrote figures/comparator_combined.png")
