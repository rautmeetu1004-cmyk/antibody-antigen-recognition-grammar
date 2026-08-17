"""Phase 4 - is the recognition chemistry conserved between the two groups?

Phase 3 reported r = 0.42 between the antibody and PPI grammars, and r = 0.84 once
collapsed to chemical classes. A bare correlation is not an answer, so this phase
asks the question four ways:

  1. WITH A CONFIDENCE INTERVAL. r is a statistic like any other. Both arms'
     units are resampled independently and r recomputed, giving a CI on the
     similarity itself. Without this, "r = 0.42" cannot be distinguished from
     "r = 0.42 +- 0.30".
  2. PER BOND TYPE. If the two groups share an H-bond vocabulary but differ in
     packing, the aggregate r hides it. The nucleic-acid arm found H-bond maps
     conserved at r = 0.87-0.99 while whole-map grammar diverged; the same
     dissection is applied here.
  3. UNIVERSAL vs ARM-SPECIFIC PAIRS. A cell is UNIVERSAL when it is BH-significant
     AND bootstrap-sign-stable AND same-signed in both arms -- the protein-protein
     analogue of the nucleic-acid arm's 18 universal bonds.
  4. WITH CYSTEINE MASKED. Inter-chain disulfides are covalent bonds, not
     recognition preferences, and C-C is the single most extreme PPI cell
     (log2 OR +2.21). Every headline number is reported with and without Cys so
     the reader can see how much of the "conservation" is one covalent artefact.

Outputs
  figures/phase4_conservation.png
  metadata/analysis/conservation.json
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

N_BOOT_R = 1000
CYS = C.AAS.index("C")


def _mask_cys(M: np.ndarray) -> np.ndarray:
    out = M.copy()
    out[CYS, :] = np.nan
    out[:, CYS] = np.nan
    return out


def arm_maps(arm: str, level: str = "cluster30"):
    rows = C.load_contacts(arm, level)
    sym = (arm == "ppi")
    out = {"rows": rows, "sym": sym}
    out["full"] = C.log2_odds(C.contact_matrix(rows, symmetrise=sym))
    _, q = C.fisher_bh(C.contact_matrix(C.dedup_rows(arm, rows), weighted=False,
                                        symmetrise=sym))
    _, _, _, same = C.bootstrap_ci(rows, symmetrise=sym)
    out["q"], out["same"] = q, same
    for b in C.BOND_COLS:
        out[b] = C.log2_odds(C.contact_matrix(rows, symmetrise=sym, bond=b))
        out[f"{b}_n"] = int(C.contact_matrix(rows, weighted=False, symmetrise=sym,
                                             bond=b).sum())
    return out


def r_with_ci(rows_a, rows_b, bond=None, mask_cys=False, n_boot=N_BOOT_R,
              seed=C.SEED):
    """Pearson r between the two arms' maps, with a unit-bootstrap CI.

    Both arms are resampled independently because they are independent corpora;
    the CI therefore propagates sampling error from BOTH sides, which a naive
    correlation test would ignore entirely.
    """
    rng = np.random.default_rng(seed)
    ids_a, Ua = C.unit_matrices(rows_a, symmetrise=False, bond=bond)
    ids_b, Ub = C.unit_matrices(rows_b, symmetrise=True, bond=bond)
    na, nb = len(ids_a), len(ids_b)
    if na < 2 or nb < 2:
        return {"r": float("nan"), "ci": [float("nan")] * 2, "n_cells": 0}

    def _lo2(stack):
        M = stack.sum(0)
        return _mask_cys(C.log2_odds(M)) if mask_cys else C.log2_odds(M)

    point = C.matrix_similarity(_lo2(Ua), _lo2(Ub))

    fa, fb = Ua.reshape(na, 400), Ub.reshape(nb, 400)
    ca = rng.multinomial(na, np.full(na, 1 / na), size=n_boot).astype(float)
    cb = rng.multinomial(nb, np.full(nb, 1 / nb), size=n_boot).astype(float)
    A = (ca @ fa).reshape(n_boot, 20, 20)
    B = (cb @ fb).reshape(n_boot, 20, 20)

    rs = np.empty(n_boot)
    for i in range(n_boot):
        a = C.log2_odds(A[i]); b = C.log2_odds(B[i])
        if mask_cys:
            a, b = _mask_cys(a), _mask_cys(b)
        rs[i] = C.matrix_similarity(a, b)
    ok = np.isfinite(rs)
    n_cells = int((np.isfinite(_lo2(Ua)) & np.isfinite(_lo2(Ub))).sum())
    return {"r": float(point),
            "ci": [float(np.percentile(rs[ok], 2.5)),
                   float(np.percentile(rs[ok], 97.5))],
            "n_cells": n_cells}


def classify_cells(A, B):
    """Universal / arm-specific / divergent, using significance AND sign stability."""
    ab, pp = A["full"], B["full"]
    sig_a = np.isfinite(ab) & (A["q"] < 0.05) & (A["same"] >= 0.95)
    sig_p = np.isfinite(pp) & (B["q"] < 0.05) & (B["same"] >= 0.95)
    same_sign = np.sign(ab) == np.sign(pp)

    universal, ab_only, ppi_only, opposed = [], [], [], []
    for i in range(20):
        for j in range(20):
            rec = {"pair": f"{C.AAS[i]}-{C.AAS[j]}",
                   "log2_or_antibody": round(float(ab[i, j]), 3),
                   "log2_or_ppi": round(float(pp[i, j]), 3)}
            if sig_a[i, j] and sig_p[i, j]:
                (universal if same_sign[i, j] else opposed).append(rec)
            elif sig_a[i, j]:
                ab_only.append(rec)
            elif sig_p[i, j]:
                ppi_only.append(rec)
    key = lambda r: -abs(r["log2_or_antibody"])          # noqa: E731
    return {
        "n_universal": len(universal),
        "n_antibody_specific": len(ab_only),
        "n_ppi_specific": len(ppi_only),
        "n_opposed": len(opposed),
        "universal": sorted(universal, key=key)[:25],
        "antibody_specific": sorted(ab_only, key=key)[:20],
        "opposed": sorted(opposed, key=key)[:20],
    }


def figure(A, B, res):
    fig = plt.figure(figsize=(15.5, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1.1], wspace=0.3)

    # -- scatter: every cell, antibody vs PPI, coloured by probe chemical class
    ax = fig.add_subplot(gs[0, 0])
    ab, pp = A["full"], B["full"]
    for gi, g in enumerate(C.GROUPS):
        xs, ys = [], []
        for i, a in enumerate(C.AAS):
            if C.AA_GROUP[a] != g:
                continue
            for j in range(20):
                if np.isfinite(ab[i, j]) and np.isfinite(pp[i, j]):
                    xs.append(pp[i, j]); ys.append(ab[i, j])
        ax.scatter(xs, ys, s=14, alpha=0.7, color=C.GROUP_COLOR[g], label=g,
                   edgecolors="none")
    lim = [-2.2, 2.6]
    ax.plot(lim, lim, ls="--", c="0.5", lw=1)
    ax.axhline(0, c="0.8", lw=0.7); ax.axvline(0, c="0.8", lw=0.7)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("PPI log$_2$(obs/exp)"); ax.set_ylabel("Antibody log$_2$(obs/exp)")
    r = res["overall"]["all_contacts"]
    ax.set_title(f"Cell-by-cell agreement\nr = {r['r']:.2f} "
                 f"[{r['ci'][0]:.2f}, {r['ci'][1]:.2f}]", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="upper left")

    # -- per-bond correlation with CI
    ax = fig.add_subplot(gs[0, 1])
    keys = ["all_contacts", "vdw", "hbond", "salt"]
    labels = ["All contacts", "van der Waals", "H-bond", "Salt bridge"]
    ys = np.arange(len(keys))
    for k, y in zip(keys, ys):
        d = res["overall"][k]
        ax.plot([d["ci"][0], d["ci"][1]], [y, y], c="0.35", lw=2)
        ax.plot(d["r"], y, "o", ms=8, c="#b2182b")
        dm = res["cys_masked"][k]
        ax.plot(dm["r"], y + 0.22, "s", ms=6, c="#2166ac")
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(0, c="0.8", lw=0.8)
    ax.set_xlabel("Antibody vs PPI grammar correlation (r)")
    ax.set_title("Conservation by bond chemistry\n"
                 "red = all cells, blue = Cys masked", fontsize=10)
    ax.set_xlim(-0.2, 1.0)
    ax.invert_yaxis()

    # -- cell classification
    ax = fig.add_subplot(gs[0, 2])
    cats = ["n_universal", "n_antibody_specific", "n_ppi_specific", "n_opposed"]
    names = ["Universal\n(both, same sign)", "Antibody-\nspecific",
             "PPI-\nspecific", "Opposed\n(both, opposite sign)"]
    vals = [res["cells"][c] for c in cats]
    ax.bar(range(4), vals, color=["#1b7837", "#b2182b", "#2166ac", "#762a83"])
    for i, v in enumerate(vals):
        ax.text(i, v + 3, str(v), ha="center", fontsize=9)
    ax.set_xticks(range(4)); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("cells (of 400)")
    ax.set_title("Which cells are shared?\nBH q<0.05 AND sign-stable", fontsize=10)

    fig.suptitle("Phase 4 - conservation of recognition chemistry between "
                 "antibody-antigen and general protein-protein interfaces",
                 fontsize=12.5)
    fig.savefig(C.OUT_FIG / "phase4_conservation.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    A = arm_maps("antibody")
    B = arm_maps("ppi")

    res = {"overall": {}, "cys_masked": {}, "class_level": {}}
    for key, bond in (("all_contacts", None), ("vdw", "vdw"),
                      ("hbond", "hbond"), ("salt", "salt")):
        res["overall"][key] = r_with_ci(A["rows"], B["rows"], bond=bond)
        res["cys_masked"][key] = r_with_ci(A["rows"], B["rows"], bond=bond,
                                           mask_cys=True)
        res["overall"][key]["n_contacts_antibody"] = (
            len(A["rows"]) if bond is None else A[f"{bond}_n"])
        res["overall"][key]["n_contacts_ppi"] = (
            len(B["rows"]) if bond is None else B[f"{bond}_n"])

    res["class_level"]["r"] = C.matrix_similarity(
        C.log2_odds(C.group_matrix(C.contact_matrix(A["rows"]))),
        C.log2_odds(C.group_matrix(C.contact_matrix(B["rows"], symmetrise=True))))

    res["cells"] = classify_cells(A, B)
    res["interpretation_note"] = (
        "r is computed over the 400 cells of the log2-odds maps. The CI resamples "
        "cluster30 units in BOTH arms independently. A high per-bond r with a low "
        "overall r would mean the groups share a bond vocabulary but deploy it "
        "differently; the reverse would mean the aggregate agreement is carried by "
        "one dominant bond type.")

    (C.OUT_MET / "conservation.json").write_text(json.dumps(res, indent=1))
    figure(A, B, res)

    print(f"{'map':14s} {'r':>7s} {'95% CI':>18s} {'r (Cys masked)':>16s}")
    for k in ("all_contacts", "vdw", "hbond", "salt"):
        d, m = res["overall"][k], res["cys_masked"][k]
        print(f"{k:14s} {d['r']:7.3f}  [{d['ci'][0]:6.3f},{d['ci'][1]:6.3f}]  "
              f"{m['r']:16.3f}")
    print(f"\nchemical-class r = {res['class_level']['r']:.3f}")
    c = res["cells"]
    print(f"universal={c['n_universal']}  antibody-specific={c['n_antibody_specific']}"
          f"  ppi-specific={c['n_ppi_specific']}  opposed={c['n_opposed']}")
    print("\ntop universal pairs:",
          ", ".join(f"{r['pair']}({r['log2_or_antibody']:+.2f}/"
                    f"{r['log2_or_ppi']:+.2f})" for r in c["universal"][:10]))
    print("\nopposed pairs:",
          ", ".join(f"{r['pair']}({r['log2_or_antibody']:+.2f}/"
                    f"{r['log2_or_ppi']:+.2f})" for r in c["opposed"][:10]))
