"""Regenerate every figure in publication form: no baked-in titles, panel letters.

Journals expect the caption to carry the description and the artwork to carry only
data and panel identifiers. The analysis scripts were written with titles on the
axes, which is useful while working and wrong for submission.

Rather than editing twenty scripts and risking divergence between the working and
published versions, this patches matplotlib for the duration of the run:

  * Axes.set_title and Figure.suptitle become no-ops, so every title disappears
  * Figure.savefig inserts a bold (a), (b), (c) ... into the upper left of each
    visible axes, in creation order, whenever a figure has more than one panel
  * fig.text annotations are left alone, since those carry legend keys rather than
    titles

Each script is then executed unchanged. Two scripts whose figures depend on an
expensive structural pass are replotted from their cached result tables instead.

Run:  python3 code/make_publication_figures.py
"""
from __future__ import annotations

import pathlib
import runpy
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aa_common as C

PANEL_LETTERS = "abcdefghijklmnopqrstuvwxyz"

# scripts that draw the figures, cheapest first; the expensive structural passes are
# handled separately below
SCRIPTS = [
    "phase2_inventory.py",
    "aa_heatmap.py",
    "phase4_conservation.py",
    "phase5_readout.py",
    "phase8_differential.py",
    "case_thermodynamics.py",
    "phase6b_functional_groups.py",
    "phase6c_conformation.py",
    "phase11_controls.py",
    "add1_transfer.py",
    "immunodominance.py",
    "phase6d_folds.py",
    "phase7_paratope.py",
    "add6_functional_sites.py",
    "add4_pca_geometry.py",
    "add5_epitope_logic.py",
    "add7_review_response.py",
    "phase9_modified.py",
    "case_lysozyme.py",
    "case_matched_targets.py",
]


_ORIGINAL_SAVEFIG = matplotlib.figure.Figure.savefig


def install_patches():
    """Suppress titles and stamp panel letters at save time."""
    matplotlib.axes.Axes.set_title = lambda self, *a, **k: None
    matplotlib.figure.Figure.suptitle = lambda self, *a, **k: None

    original_savefig = _ORIGINAL_SAVEFIG

    def savefig(self, *args, **kwargs):
        # a colorbar is an Axes with data, but it is not a panel and must not be
        # lettered; matplotlib tags the axes it creates for a Colorbar
        axes = [ax for ax in self.get_axes()
                if ax.get_visible() and ax.has_data()
                and getattr(ax, "_colorbar", None) is None
                and not getattr(ax, "_colorbar_info", None)]
        if len(axes) > 1:
            for i, ax in enumerate(axes[:len(PANEL_LETTERS)]):
                ax.text(-0.08, 1.06, f"({PANEL_LETTERS[i]})",
                        transform=ax.transAxes, fontsize=13, fontweight="bold",
                        va="top", ha="left", zorder=1000)
        return original_savefig(self, *args, **kwargs)

    matplotlib.figure.Figure.savefig = savefig


def replot_thermodynamics():
    """Addition 2 figure from cached per-interface values, no SASA recomputation."""
    import csv
    from collections import defaultdict
    from scipy.stats import mannwhitneyu

    rows = list(csv.DictReader((C.OUT_INT / "thermodynamics.csv").open()))
    units = {}
    for arm in C.ARMS:
        for p, v in C.units(arm).items():
            units[(arm, p)] = (v.get("cluster30") if arm == "antibody"
                               else (v.get("a_to_b") or {}).get("cluster30"))
    per = {arm: defaultdict(lambda: defaultdict(list)) for arm in C.ARMS}
    for r in rows:
        u = units.get((r["arm"], r["pdb_id"]))
        if not u:
            continue
        ap, po = float(r["dA_apolar"]), float(r["dA_polar"])
        tot = ap + po
        if tot <= 0:
            continue
        dcp = 0.45 * (-ap) + (-0.26) * (-po)
        per[r["arm"]][u]["frac_apolar_buried"].append(ap / tot)
        per[r["arm"]][u]["dCp_per_1000A2"].append(dcp / (tot / 1000.0))
        per[r["arm"]][u]["dA_total"].append(tot)

    fields = [("frac_apolar_buried", "fraction of buried surface that is apolar"),
              ("dCp_per_1000A2",
               "$\\Delta C_p$ per 1000 $\\AA^2$ buried (cal mol$^{-1}$ K$^{-1}$)"),
              ("dA_total", "total buried surface ($\\AA^2$)")]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    for ax, (f, lab) in zip(axes, fields):
        vals = {}
        for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
            v = np.array([np.mean(d[f]) for d in per[arm].values() if d[f]])
            vals[arm] = v
            ax.hist(v, bins=35, density=True, alpha=0.55, color=col,
                    label=C.ARM_LABEL[arm])
            ax.axvline(np.median(v), color=col, ls="--", lw=1.3)
        ax.set_xlabel(lab, fontsize=9)
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "add2_thermodynamics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  replotted add2_thermodynamics.png")


def replot_shape_complementarity():
    """Addition 3 figure from cached per-interface Sc values."""
    import csv
    from collections import defaultdict

    rows = list(csv.DictReader((C.OUT_INT / "shape_complementarity.csv").open()))
    units = {}
    for arm in C.ARMS:
        for p, v in C.units(arm).items():
            units[(arm, p)] = (v.get("cluster30") if arm == "antibody"
                               else (v.get("a_to_b") or {}).get("cluster30"))
    per = {arm: defaultdict(list) for arm in C.ARMS}
    for r in rows:
        u = units.get((r["arm"], r["pdb_id"]))
        if u:
            per[r["arm"]][u].append(float(r["sc"]))
    fig, ax = plt.subplots(figsize=(7.2, 5))
    for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
        v = np.array([np.mean(x) for x in per[arm].values() if x])
        ax.hist(v, bins=35, density=True, alpha=0.55, color=col,
                label=f"{C.ARM_LABEL[arm]} (median {np.median(v):.3f})")
        ax.axvline(np.median(v), color=col, ls="--", lw=1.3)
    ax.set_xlabel("shape complementarity $S_c$")
    ax.set_ylabel("density")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "add3_shape_complementarity.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("  replotted add3_shape_complementarity.png")



def remontage(panels, out_name, legend=None, figsize_per=5.6):
    """Rebuild a PyMOL montage from rendered panels: no title, panel letters only."""
    import matplotlib.image as mpimg
    from matplotlib.patches import Patch
    present = [(p, lab) for p, lab in panels if (C.OUT_FIG / p).exists()]
    if not present:
        print("  no panels for " + out_name)
        return
    fig, axes = plt.subplots(1, len(present),
                             figsize=(figsize_per * len(present), figsize_per + 0.6))
    if len(present) == 1:
        axes = [axes]
    for i, (ax, (png, lab)) in enumerate(zip(axes, present)):
        ax.imshow(mpimg.imread(str(C.OUT_FIG / png)))
        ax.axis("off")
        ax.text(0.01, 0.99, f"({PANEL_LETTERS[i]})", transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top", ha="left")
        if lab:
            ax.text(0.5, -0.02, lab, transform=ax.transAxes, fontsize=9,
                    ha="center", va="top")
    if legend:
        fig.legend(handles=[Patch(facecolor=c, label=l) for c, l in legend],
                   loc="lower center", ncol=len(legend), fontsize=9, frameon=False)
        fig.tight_layout(rect=[0, 0.07, 1, 1])
    else:
        fig.tight_layout()
    # this function letters its own panels, so save through the unpatched method;
    # going through the patch would stamp a second set of letters on every panel
    _ORIGINAL_SAVEFIG(fig, C.OUT_FIG / out_name, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  rebuilt " + out_name)


def remontage_all():
    remontage([("pymol_mechanism_3D9A.png", "3D9A, 1.20 A, lysozyme C"),
               ("pymol_mechanism_2P45.png", "2P45, 1.10 A, ribonuclease"),
               ("pymol_mechanism_7AQZ.png", "7AQZ, 1.30 A, nanobody complex")],
              "pymol_mechanism_montage.png")
    remontage([("pymol_immunodominance.png", "front"),
               ("pymol_immunodominance_back.png", "rotated 180 degrees")],
              "pymol_immunodominance_montage.png",
              legend=[("#FFFFCC", "rarely contacted"),
                      ("#FF8C00", "frequently contacted"),
                      ("#4169E1", "catalytic residues")])


if __name__ == "__main__":
    install_patches()
    t0 = time.time()
    ok, failed = [], []
    for name in SCRIPTS:
        path = HERE / name
        if not path.exists():
            failed.append((name, "missing"))
            continue
        print(f"[{time.time()-t0:6.0f}s] {name}", flush=True)
        try:
            runpy.run_path(str(path), run_name="__main__")
            ok.append(name)
            plt.close("all")
        except SystemExit:
            ok.append(name)
        except Exception as e:                                   # noqa: BLE001
            failed.append((name, f"{type(e).__name__}: {e}"))
            print(f"    FAILED {type(e).__name__}: {e}", flush=True)
            plt.close("all")

    print("\nreplotting from cached tables:")
    for fn in (replot_thermodynamics, replot_shape_complementarity,
               remontage_all):
        try:
            fn()
            ok.append(fn.__name__)
        except Exception as e:                                   # noqa: BLE001
            failed.append((fn.__name__, f"{type(e).__name__}: {e}"))
            print(f"    FAILED {fn.__name__}: {e}")

    print(f"\n{len(ok)} succeeded, {len(failed)} failed, "
          f"{time.time()-t0:.0f}s total")
    for n, why in failed:
        print(f"  FAILED {n}: {why}")
