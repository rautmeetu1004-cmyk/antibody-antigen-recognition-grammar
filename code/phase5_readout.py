"""Phase 5 - does the side partner matter? Side-chain vs main-chain readout.

A protein surface can be read two ways, and they mean different things:

  SIDE-CHAIN readout   identity-specific. Reading Arg's guanidinium is reading
                       "this is an arginine". This is sequence readout.
  MAIN-CHAIN readout   identity-agnostic. The backbone N-H and C=O of a leucine are
                       chemically identical to those of a lysine; only the local
                       conformation differs. This is shape readout.

This is the direct protein-protein analogue of the nucleic-acid arm's base-versus-
backbone decomposition, and it bears on why antibodies work: a binder that reads
mostly backbone is reading fold and loop conformation, and would generalise across
sequences; a binder that reads mostly side chains is reading sequence identity, and
must discriminate between similar surfaces.

The contact tables carry the antigen-side moiety breakdown per residue pair
(`ag_side_contacts` / `ag_main_contacts`, and the same split for H-bonds and vdW),
so each contact is apportioned between the two channels rather than assigned
whole to one.

Outputs
  figures/phase5_readout.png
  metadata/analysis/readout.json
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C


def split_matrices(arm: str, level: str = "cluster30"):
    """Two 20x20 maps: contacts apportioned to target side-chain vs main-chain.

    A pair contributing 5 side-chain and 3 main-chain atom contacts adds 5/8 of its
    unit weight to the side map and 3/8 to the main map. Apportioning rather than
    hard-assigning avoids inventing a threshold, and keeps the two maps summing to
    the full map.
    """
    rows = C.load_contacts(arm, level)
    sym = (arm == "ppi")
    idx = {a: i for i, a in enumerate(C.AAS)}
    side = np.zeros((20, 20))
    main = np.zeros((20, 20))
    per_res_side = defaultdict(float)
    per_res_main = defaultdict(float)
    unit_side = defaultdict(float)
    unit_tot = defaultdict(float)

    for r in rows:
        s = float(r.get("ag_side_contacts") or 0)
        m = float(r.get("ag_main_contacts") or 0)
        tot = s + m
        if tot <= 0:
            continue
        w = r["weight"]
        i, j = idx[r["aa_probe"]], idx[r["aa_target"]]
        side[i, j] += w * s / tot
        main[i, j] += w * m / tot
        per_res_side[r["aa_target"]] += w * s / tot
        per_res_main[r["aa_target"]] += w * m / tot
        unit_side[r["unit"]] += w * s / tot
        unit_tot[r["unit"]] += w

    if sym:
        side = (side + side.T) / 2.0
        main = (main + main.T) / 2.0

    per_res = {a: {"side": per_res_side.get(a, 0.0), "main": per_res_main.get(a, 0.0)}
               for a in C.AAS}
    for a, d in per_res.items():
        t = d["side"] + d["main"]
        d["f_side"] = d["side"] / t if t else np.nan
    unit_f = np.array([unit_side[u] / unit_tot[u] for u in unit_tot if unit_tot[u] > 0])
    return side, main, per_res, unit_f


def figure(data):
    fig = plt.figure(figsize=(16, 8.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.15], hspace=0.35, wspace=0.28)

    for r, arm in enumerate(C.ARMS):
        side, main, per_res, _ = data[arm]
        for c, (M, name) in enumerate(((side, "side-chain-mediated"),
                                       (main, "main-chain-mediated"))):
            ax = fig.add_subplot(gs[r, c])
            lo2 = C.log2_odds(M)
            im = ax.imshow(np.ma.masked_invalid(lo2), cmap="RdBu_r",
                           norm=TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5))
            ax.set_xticks(range(20)); ax.set_yticks(range(20))
            ax.set_xticklabels(C.AAS, fontsize=6)
            ax.set_yticklabels(C.AAS, fontsize=6)
            for t, a in zip(ax.get_xticklabels(), C.AAS):
                t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
            for t, a in zip(ax.get_yticklabels(), C.AAS):
                t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
            ax.set_title(f"{C.ARM_LABEL[arm]}\n{name}", fontsize=9)
            if c == 0:
                ax.set_ylabel("probe residue", fontsize=8)
            ax.set_xlabel("target residue", fontsize=8)

    # per-target-residue side-chain fraction, both arms
    ax = fig.add_subplot(gs[:, 2])
    ys = np.arange(20)
    for arm, col, off in (("antibody", "#b2182b", -0.19), ("ppi", "#2166ac", 0.19)):
        per_res = data[arm][2]
        vals = [per_res[a]["f_side"] for a in C.AAS]
        ax.barh(ys + off, vals, height=0.36, color=col, label=C.ARM_LABEL[arm])
    ax.set_yticks(ys); ax.set_yticklabels(C.AAS, fontsize=8)
    for t, a in zip(ax.get_yticklabels(), C.AAS):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    ax.invert_yaxis()
    ax.axvline(0.5, c="0.4", ls="--", lw=1)
    ax.set_xlabel("fraction of contacts made through the residue's SIDE CHAIN")
    ax.set_title("How each target residue is read\n"
                 "(left of dashed line = read mostly by backbone)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower right")

    fig.colorbar(im, ax=fig.axes[:4], fraction=0.02, pad=0.02).set_label(
        "log$_2$(obs/exp)", fontsize=8)
    fig.suptitle("Phase 5 - sequence readout (side chain) vs shape readout "
                 "(main chain)", fontsize=13)
    fig.savefig(C.OUT_FIG / "phase5_readout.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    data = {arm: split_matrices(arm) for arm in C.ARMS}

    out = {"per_residue": {}, "maps": {}}
    for arm in C.ARMS:
        side, main, per_res, unit_f = data[arm]
        st, mt = side.sum(), main.sum()
        out["maps"][arm] = {
            "f_side_overall": float(st / (st + mt)),
            "side_vs_main_map_r": C.matrix_similarity(C.log2_odds(side),
                                                      C.log2_odds(main)),
            "n_units": int(len(unit_f)),
            "unit_f_side_median": float(np.median(unit_f)),
        }
        out["per_residue"][arm] = {a: round(per_res[a]["f_side"], 4) for a in C.AAS}

    # do the two arms read the same residues the same way?
    ab = np.array([out["per_residue"]["antibody"][a] for a in C.AAS])
    pp = np.array([out["per_residue"]["ppi"][a] for a in C.AAS])
    ok = np.isfinite(ab) & np.isfinite(pp)
    out["per_residue_correlation"] = float(np.corrcoef(ab[ok], pp[ok])[0, 1])
    out["per_residue_delta"] = {a: round(float(ab[i] - pp[i]), 4)
                                for i, a in enumerate(C.AAS)}

    u_ab, u_pp = data["antibody"][3], data["ppi"][3]
    u, p = mannwhitneyu(u_ab, u_pp, alternative="two-sided")
    out["unit_level_test"] = {
        "median_antibody": float(np.median(u_ab)),
        "median_ppi": float(np.median(u_pp)),
        "p": float(p),
        "effect_rank_biserial": float(2 * u / (len(u_ab) * len(u_pp)) - 1),
        "n_units_antibody": len(u_ab), "n_units_ppi": len(u_pp),
    }
    out["cross_arm_map_r"] = {
        "side": C.matrix_similarity(C.log2_odds(data["antibody"][0]),
                                    C.log2_odds(data["ppi"][0])),
        "main": C.matrix_similarity(C.log2_odds(data["antibody"][1]),
                                    C.log2_odds(data["ppi"][1])),
    }
    out["note"] = (
        "main-chain readout is identity-agnostic by construction: the backbone of "
        "every residue except proline and glycine is chemically identical, so a "
        "main-chain map's residue axis reports WHICH residues present their "
        "backbone, not which chemistry is being read.")

    (C.OUT_MET / "readout.json").write_text(json.dumps(out, indent=1))
    figure(data)

    for arm in C.ARMS:
        m = out["maps"][arm]
        print(f"{arm:9s} side-chain-mediated fraction = {m['f_side_overall']:.3f}  "
              f"(unit median {m['unit_f_side_median']:.3f})  "
              f"side-vs-main map r = {m['side_vs_main_map_r']:.3f}")
    t = out["unit_level_test"]
    print(f"\nunit-level side fraction: antibody {t['median_antibody']:.3f} vs "
          f"PPI {t['median_ppi']:.3f}  effect={t['effect_rank_biserial']:+.3f}  "
          f"p={t['p']:.2e}")
    print(f"cross-arm map r: side={out['cross_arm_map_r']['side']:.3f}  "
          f"main={out['cross_arm_map_r']['main']:.3f}")
    print(f"per-residue side-fraction correlation between arms = "
          f"{out['per_residue_correlation']:.3f}")
    d = out["per_residue_delta"]
    top = sorted(d.items(), key=lambda kv: -abs(kv[1]))[:8]
    print("largest per-residue differences (antibody - PPI):",
          ", ".join(f"{a}{v:+.3f}" for a, v in top))
