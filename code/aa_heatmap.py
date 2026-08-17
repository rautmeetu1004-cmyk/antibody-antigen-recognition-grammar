"""The core 20x20 amino-acid recognition matrices for both arms.

Outputs
  figures/heatmap_20x20_{cluster30,uniprot}.png   the two arms side by side
  figures/heatmap_chemical_classes.png            6x6 collapse + class marginals
  figures/heatmap_bond_decomposition.png          20x20 per bond type, both arms
  interactions/analysis/matrix_{arm}_{level}.csv  the numbers behind every figure
  metadata/analysis/heatmap.json                  stats roll-up

Reading the figure: colour is log2(observed/expected) where expected comes from
the two residues' own interface abundances. Red = the pair is enriched over what
its two residues' commonness predicts; blue = depleted. A cell is stippled when
it fails BH q<0.05, and outlined when the unit bootstrap keeps its sign in >=95%
of draws. Effect size and evidence are deliberately shown as separate channels.
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import aa_common as C


def _axis_labels(ax, xlab, ylab):
    ax.set_xticks(range(20))
    ax.set_yticks(range(20))
    ax.set_xticklabels(C.AAS, fontsize=7)
    ax.set_yticklabels(C.AAS, fontsize=7)
    for t, a in zip(ax.get_xticklabels(), C.AAS):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    for t, a in zip(ax.get_yticklabels(), C.AAS):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel(ylab, fontsize=9)


def _draw(ax, lo2, q, same, title, xlab, ylab, vmax=1.5):
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.imshow(np.ma.masked_invalid(lo2), cmap="RdBu_r", norm=norm)
    for i in range(20):
        for j in range(20):
            if not np.isfinite(lo2[i, j]) or q[i, j] >= 0.05:
                # stipple = not significant after BH over all 400 cells
                ax.plot(j, i, marker="x", ms=2.2, color="0.45", mew=0.5)
            elif same[i, j] >= 0.95:
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                       ec="k", lw=0.45))
    _axis_labels(ax, xlab, ylab)
    ax.set_title(title, fontsize=10)
    return im


def run(level: str = "cluster30") -> dict:
    out: dict = {"unit_level": level, "arms": {}}
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.2))

    cache = {}
    for ax, arm in zip(axes, C.ARMS):
        rows = C.load_contacts(arm, level)
        sym = (arm == "ppi")     # co-evolved partners are exchangeable; Ab/Ag are not
        Mw = C.contact_matrix(rows, weighted=True, symmetrise=sym)
        # Fisher needs the TRUE number of physical contacts, not the doubled
        # both-directions row count -- see aa_common.dedup_rows
        Mraw = C.contact_matrix(C.dedup_rows(arm, rows), weighted=False,
                                symmetrise=sym)
        lo2 = C.log2_odds(Mw)
        _, q = C.fisher_bh(Mraw)
        point, blo, bhi, same = C.bootstrap_ci(rows, symmetrise=sym)

        xlab = "Epitope (antigen) residue" if arm == "antibody" else "Partner B residue"
        ylab = "Paratope (antibody) residue" if arm == "antibody" else "Partner A residue"
        sm = C.corpus_summary(arm, level)
        title = (f"{C.ARM_LABEL[arm]}\n{sm['n_structures']} structures, "
                 f"{sm['n_units']} {level} units, {sm['n_contacts']:,} contacts"
                 + ("  (symmetrised)" if sym else ""))
        im = _draw(ax, lo2, q, same, title, xlab, ylab)

        cache[arm] = dict(lo2=lo2, q=q, same=same, Mw=Mw, Mraw=Mraw,
                          blo=blo, bhi=bhi)

        sig = np.isfinite(lo2) & (q < 0.05)
        out["arms"][arm] = {
            **sm,
            "n_cells_finite": int(np.isfinite(lo2).sum()),
            "n_cells_bh_sig": int(sig.sum()),
            "n_cells_sign_stable": int((same >= 0.95).sum()),
            "log2_or_range": [float(np.nanmin(lo2)), float(np.nanmax(lo2))],
            "top_enriched": _top(lo2, q, same, 12, True),
            "top_depleted": _top(lo2, q, same, 12, False),
        }
        _write_matrix(arm, level, lo2, q, same, Mw, Mraw)

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("log$_2$(observed / expected)", fontsize=9)
    fig.suptitle("Amino-acid recognition grammar: antibody-antigen vs general "
                 "protein-protein interfaces", fontsize=12.5)
    fig.text(0.5, 0.015, "x = fails BH q<0.05 (400 tests/panel)   |   black outline = "
             "bootstrap sign stable in >=95% of 2,000 unit resamples",
             ha="center", fontsize=8, color="0.3")
    fig.savefig(C.OUT_FIG / f"heatmap_20x20_{level}.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    out["grammar_similarity_antibody_vs_ppi"] = C.matrix_similarity(
        cache["antibody"]["lo2"], cache["ppi"]["lo2"])
    return out, cache


def _top(lo2, q, same, n, enriched):
    cells = []
    for i in range(20):
        for j in range(20):
            if np.isfinite(lo2[i, j]) and q[i, j] < 0.05 and same[i, j] >= 0.95:
                cells.append((lo2[i, j], C.AAS[i], C.AAS[j], q[i, j]))
    cells.sort(reverse=enriched)
    return [{"pair": f"{a}-{b}", "log2_or": round(float(v), 3), "q": float(qq)}
            for v, a, b, qq in cells[:n]]


def _write_matrix(arm, level, lo2, q, same, Mw, Mraw):
    path = C.OUT_INT / f"matrix_{arm}_{level}.csv"
    with path.open("w") as fh:
        fh.write("probe_aa,target_aa,log2_or,q_bh,boot_sign_frac,"
                 "weighted_count,raw_count\n")
        for i, a in enumerate(C.AAS):
            for j, b in enumerate(C.AAS):
                fh.write(f"{a},{b},{lo2[i,j]:.4f},{q[i,j]:.3e},{same[i,j]:.3f},"
                         f"{Mw[i,j]:.4f},{int(Mraw[i,j])}\n")


def chemical_classes(cache) -> dict:
    """6x6 chemical-class view -- where the shared-chemistry question is powered."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    res = {}
    for ax, arm in zip(axes, C.ARMS):
        G = C.group_matrix(cache[arm]["Mw"])
        lo2 = C.log2_odds(G)
        norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
        im = ax.imshow(lo2, cmap="RdBu_r", norm=norm)
        for i in range(6):
            for j in range(6):
                ax.text(j, i, f"{lo2[i,j]:+.2f}", ha="center", va="center",
                        fontsize=7.5,
                        color="w" if abs(lo2[i, j]) > 0.6 else "k")
        ax.set_xticks(range(6)); ax.set_yticks(range(6))
        ax.set_xticklabels(C.GROUPS, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(C.GROUPS, fontsize=8)
        for t, g in zip(ax.get_xticklabels(), C.GROUPS):
            t.set_color(C.GROUP_COLOR[g])
        for t, g in zip(ax.get_yticklabels(), C.GROUPS):
            t.set_color(C.GROUP_COLOR[g])
        ax.set_title(C.ARM_LABEL[arm], fontsize=10)
        res[arm] = {f"{C.GROUPS[i]}|{C.GROUPS[j]}": round(float(lo2[i, j]), 3)
                    for i in range(6) for j in range(6)}
        res[f"{arm}_matrix"] = lo2.tolist()
    fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02).set_label(
        "log$_2$(obs/exp)", fontsize=9)
    fig.suptitle("Chemical-class recognition (6x6 collapse of the 20x20)", fontsize=12)
    fig.savefig(C.OUT_FIG / "heatmap_chemical_classes.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    res["class_similarity_antibody_vs_ppi"] = C.matrix_similarity(
        np.array(res["antibody_matrix"]), np.array(res["ppi_matrix"]))
    return res


def bond_decomposition(level="cluster30") -> dict:
    """The same matrix split by bond chemistry -- salt / H-bond / vdW / stacking."""
    bonds = ["salt", "hbond", "vdw", "stack"]
    fig, axes = plt.subplots(2, 4, figsize=(20, 9.5))
    res = {}
    for r, arm in enumerate(C.ARMS):
        rows = C.load_contacts(arm, level)
        sym = (arm == "ppi")
        full = C.log2_odds(C.contact_matrix(rows, symmetrise=sym))
        for c, bond in enumerate(bonds):
            ax = axes[r, c]
            Mw = C.contact_matrix(rows, symmetrise=sym, bond=bond)
            Mraw = C.contact_matrix(C.dedup_rows(arm, rows), weighted=False,
                                    symmetrise=sym, bond=bond)
            lo2 = C.log2_odds(Mw)
            _, q = C.fisher_bh(Mraw)
            same = np.ones_like(lo2)          # bootstrap is run only on the full map
            _draw(ax, lo2, q, same, f"{C.ARM_LABEL[arm]} - {bond}", "", "", vmax=2.0)
            n = int(Mraw.sum())
            res[f"{arm}_{bond}"] = {
                "n_contacts_raw": n,
                "frac_of_all_contacts": round(n / max(1, len(rows)), 4),
                "similarity_to_full_map": C.matrix_similarity(lo2, full),
                "top": _top(lo2, q, same, 8, True),
            }
    fig.suptitle("Recognition grammar decomposed by bond chemistry", fontsize=13)
    fig.savefig(C.OUT_FIG / "heatmap_bond_decomposition.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    return res


if __name__ == "__main__":
    summary = {}
    primary, cache = run("cluster30")
    summary["cluster30"] = primary
    sens, _ = run("uniprot")
    summary["uniprot"] = sens
    summary["sensitivity_note"] = (
        "cluster30 is primary; uniprot is the sensitivity analysis. If the two "
        "maps disagree, the finding is a clustering artefact.")
    summary["chemical_classes"] = chemical_classes(cache)
    summary["bond_decomposition"] = bond_decomposition("cluster30")

    (C.OUT_MET / "heatmap.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({
        "antibody": summary["cluster30"]["arms"]["antibody"],
        "ppi": summary["cluster30"]["arms"]["ppi"],
        "grammar_r": summary["cluster30"]["grammar_similarity_antibody_vs_ppi"],
        "class_r": summary["chemical_classes"]["class_similarity_antibody_vs_ppi"],
    }, indent=1))
