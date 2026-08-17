"""Phase 8 - the formal comparison: differential map and inferential statistics.

Phases 3-6 established that the two grammars correlate at r ~ 0.4-0.6 and that the
gap survives matching on secondary structure and fold. This phase asks the question
those descriptive results cannot answer on their own:

  IS THE DIFFERENCE SIGNIFICANT, AND WHERE DOES IT LIVE?

Four tests:

  1. PERMUTATION TEST OF GRAMMAR IDENTITY. Pool the cluster30 units from both arms,
     reassign them at random to two groups of the original sizes, and recompute the
     cross-group correlation. This is the null "both arms are samples of ONE
     grammar". If the observed r sits far below that null, the grammars genuinely
     differ -- and unlike the Phase 4 bootstrap CI, this test is not attenuated by
     resampling noise, because the null is built the same way as the observation.

  2. THE DIFFERENTIAL MAP with a CI on the DIFFERENCE. Delta = log2OR(antibody) -
     log2OR(PPI), with both arms' units resampled independently, then BH across the
     400 cells. A cell only counts as differing if its CI on Delta excludes zero.

  3. ASYMMETRY TEST. The antibody matrix should be asymmetric (paratope and epitope
     are not exchangeable); the PPI matrix should not be. Observed asymmetry is
     compared against a within-arm bootstrap null, turning a design assumption into
     a measured result.

  4. DYNAMIC-RANGE TEST. Phase 4 noted antibody effect sizes are ~3x larger than PPI
     for the same residue pair. Two explanations: PPI is a heterogeneous mixture of
     binding modes whose average washes out preferences, or the antibody arm's
     smaller sample inflates extremes. Subsampling PPI down to the antibody's unit
     count separates them -- if PPI's spread rises to match, it was sample size.

Outputs
  figures/phase8_differential.png
  metadata/analysis/differential.json
  interactions/analysis/differential_cells.csv
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

N_PERM = 1000
N_BOOT = 1000


def arm_units(arm, symmetrise=True):
    rows = C.load_contacts(arm)
    return C.unit_matrices(rows, symmetrise=symmetrise)


def permutation_test(Ua, Up, n_perm=N_PERM, seed=C.SEED):
    """Null: both arms are samples of a single grammar.

    Units are the exchangeable objects, so units are what gets shuffled. Both arms
    are symmetrised here so the pooled objects are commensurable; antibody
    asymmetry is tested separately in asymmetry_test().
    """
    rng = np.random.default_rng(seed)
    na, np_ = len(Ua), len(Up)
    pool = np.concatenate([Ua, Up], axis=0)
    obs = C.matrix_similarity(C.log2_odds(Ua.sum(0)), C.log2_odds(Up.sum(0)))

    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(na + np_)
        A = pool[idx[:na]].sum(0)
        B = pool[idx[na:]].sum(0)
        null[i] = C.matrix_similarity(C.log2_odds(A), C.log2_odds(B))
    p = (np.sum(null <= obs) + 1) / (n_perm + 1)
    return {"observed_r": float(obs),
            "null_median": float(np.median(null)),
            "null_ci": [float(np.percentile(null, 2.5)),
                        float(np.percentile(null, 97.5))],
            "p_one_sided": float(p),
            "n_perm": n_perm,
            "interpretation": (
                "p is the probability of seeing a cross-group correlation this LOW "
                "if both arms were samples of one grammar")}


def differential_map(Ua, Up, n_boot=N_BOOT, seed=C.SEED):
    """Delta = log2OR(antibody) - log2OR(PPI), with a CI on the difference."""
    rng = np.random.default_rng(seed)
    na, np_ = len(Ua), len(Up)
    lo_a = C.log2_odds(Ua.sum(0))
    lo_p = C.log2_odds(Up.sum(0))
    delta = lo_a - lo_p

    fa, fp = Ua.reshape(na, 400), Up.reshape(np_, 400)
    ca = rng.multinomial(na, np.full(na, 1 / na), size=n_boot).astype(float)
    cp = rng.multinomial(np_, np.full(np_, 1 / np_), size=n_boot).astype(float)
    A = (ca @ fa).reshape(n_boot, 20, 20)
    B = (cp @ fp).reshape(n_boot, 20, 20)

    draws = np.empty((n_boot, 20, 20))
    for i in range(n_boot):
        draws[i] = C.log2_odds(A[i]) - C.log2_odds(B[i])
    with np.errstate(invalid="ignore"):
        lo = np.nanpercentile(draws, 2.5, axis=0)
        hi = np.nanpercentile(draws, 97.5, axis=0)
        # two-sided bootstrap p: how often does the difference cross zero
        frac_cross = np.nanmean(np.sign(draws) != np.sign(delta)[None, :, :], axis=0)
    p = np.clip(2 * frac_cross, 1e-4, 1.0)

    flat = p.ravel()
    order = np.argsort(flat)
    m = flat.size
    q = np.empty(m)
    prev = 1.0
    for rank, k in enumerate(order[::-1], start=1):
        prev = min(prev, flat[k] * m / (m - rank + 1))
        q[k] = prev
    return delta, lo, hi, q.reshape(20, 20), lo_a, lo_p


def asymmetry_test(arm, n_boot=500, seed=C.SEED):
    """|M - M'| / sum(M) against a within-arm bootstrap null."""
    rng = np.random.default_rng(seed)
    ids, U = C.unit_matrices(C.load_contacts(arm), symmetrise=False)
    n = len(ids)
    M = U.sum(0)
    obs = float(np.abs(M - M.T).sum() / M.sum())
    flat = U.reshape(n, 400)
    counts = rng.multinomial(n, np.full(n, 1 / n), size=n_boot).astype(float)
    draws = (counts @ flat).reshape(n_boot, 20, 20)
    vals = np.abs(draws - draws.transpose(0, 2, 1)).sum((1, 2)) / draws.sum((1, 2))
    return {"asymmetry": obs,
            "bootstrap_ci": [float(np.percentile(vals, 2.5)),
                             float(np.percentile(vals, 97.5))],
            "n_units": n}


def dynamic_range_test(Ua, Up, n_draw=200, seed=C.SEED):
    """Is PPI's compressed effect-size range real, or a sample-size artefact?"""
    rng = np.random.default_rng(seed)
    na = len(Ua)
    sd_a = float(np.nanstd(C.log2_odds(Ua.sum(0))))
    sd_p = float(np.nanstd(C.log2_odds(Up.sum(0))))
    sub = []
    for _ in range(n_draw):
        idx = rng.choice(len(Up), size=na, replace=False)
        sub.append(np.nanstd(C.log2_odds(Up[idx].sum(0))))
    sub = np.array(sub)
    return {"sd_antibody": sd_a, "sd_ppi_full": sd_p,
            "sd_ppi_subsampled_to_antibody_units_median": float(np.median(sub)),
            "sd_ppi_subsampled_ci": [float(np.percentile(sub, 2.5)),
                                     float(np.percentile(sub, 97.5))],
            "verdict": ("sample size explains it" if sd_a <= np.percentile(sub, 97.5)
                        else "antibody grammar is genuinely sharper than PPI at "
                             "equal sample size")}


if __name__ == "__main__":
    _, Ua = arm_units("antibody", symmetrise=True)
    _, Up = arm_units("ppi", symmetrise=True)
    print(f"units: antibody={len(Ua)}, ppi={len(Up)}")

    perm = permutation_test(Ua, Up)
    delta, lo, hi, q, lo_a, lo_p = differential_map(Ua, Up)
    asym = {arm: asymmetry_test(arm) for arm in C.ARMS}
    dyn = dynamic_range_test(Ua, Up)

    sig = np.isfinite(delta) & (q < 0.05) & ((lo > 0) | (hi < 0))
    cells = []
    for i in range(20):
        for j in range(20):
            if np.isfinite(delta[i, j]):
                cells.append({"pair": f"{C.AAS[i]}-{C.AAS[j]}",
                              "delta": round(float(delta[i, j]), 3),
                              "ci": [round(float(lo[i, j]), 3),
                                     round(float(hi[i, j]), 3)],
                              "q": float(q[i, j]),
                              "antibody": round(float(lo_a[i, j]), 3),
                              "ppi": round(float(lo_p[i, j]), 3),
                              "significant": bool(sig[i, j])})
    cells.sort(key=lambda c: -abs(c["delta"]))

    out = {"permutation_test": perm,
           "n_cells_significantly_different": int(sig.sum()),
           "top_differential_cells": cells[:25],
           "asymmetry": asym,
           "dynamic_range": dyn,
           "power_note": (
               f"400 cells are tested per map. The antibody arm has {len(Ua)} "
               f"cluster30 units against PPI's {len(Up)}, so the antibody map is "
               f"the power-limiting side throughout; a cell failing significance "
               f"here is not evidence of no difference.")}
    (C.OUT_MET / "differential.json").write_text(json.dumps(out, indent=1))

    with (C.OUT_INT / "differential_cells.csv").open("w") as fh:
        fh.write("pair,delta,ci_lo,ci_hi,q_bh,log2or_antibody,log2or_ppi,significant\n")
        for c in cells:
            fh.write(f"{c['pair']},{c['delta']},{c['ci'][0]},{c['ci'][1]},"
                     f"{c['q']:.3e},{c['antibody']},{c['ppi']},{int(c['significant'])}\n")

    # ------------------------------------------------------------------- figure
    fig = plt.figure(figsize=(16.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1])

    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(np.ma.masked_invalid(delta), cmap="PuOr_r",
                   norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2))
    for i in range(20):
        for j in range(20):
            if sig[i, j]:
                ax.plot(j, i, marker="o", ms=1.8, color="k")
    ax.set_xticks(range(20)); ax.set_yticks(range(20))
    ax.set_xticklabels(C.AAS, fontsize=6.5); ax.set_yticklabels(C.AAS, fontsize=6.5)
    for t, a in zip(ax.get_xticklabels(), C.AAS):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    for t, a in zip(ax.get_yticklabels(), C.AAS):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    ax.set_title(f"Differential map (antibody - PPI)\ndot = CI excludes 0 and "
                 f"BH q<0.05 ({int(sig.sum())} cells)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046).set_label("$\\Delta$ log$_2$(obs/exp)",
                                                      fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(C.SEED)
    na, np_ = len(Ua), len(Up)
    pool = np.concatenate([Ua, Up], axis=0)
    null = []
    for i in range(300):
        idx = rng.permutation(na + np_)
        null.append(C.matrix_similarity(C.log2_odds(pool[idx[:na]].sum(0)),
                                        C.log2_odds(pool[idx[na:]].sum(0))))
    ax.hist(null, bins=30, color="0.7", label="null: one shared grammar")
    ax.axvline(perm["observed_r"], color="#b2182b", lw=2.5,
               label=f"observed r={perm['observed_r']:.3f}")
    ax.set_xlabel("cross-group grammar correlation")
    ax.set_ylabel("permutations")
    ax.set_title(f"Permutation test\np = {perm['p_one_sided']:.4f}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = fig.add_subplot(gs[0, 2])
    ax.hist(np.random.default_rng(1).normal(0, 0, 1), alpha=0)   # keep axes
    labels = ["antibody", "PPI"]
    vals = [asym[a]["asymmetry"] for a in C.ARMS]
    cis = [asym[a]["bootstrap_ci"] for a in C.ARMS]
    ax.bar(range(2), vals, color=["#b2182b", "#2166ac"], width=0.5)
    for i, (v, ci) in enumerate(zip(vals, cis)):
        ax.plot([i, i], ci, color="k", lw=2)
        ax.text(i, v + 0.03, f"{v:.3f}", ha="center", fontsize=9)
    ax.set_xticks(range(2)); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("|M - M'| / sum(M)")
    ax.set_title("Matrix asymmetry\n(paratope/epitope are not exchangeable;\n"
                 "PPI partners are)", fontsize=10)

    fig.suptitle("Phase 8 - differential map and inferential statistics", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(C.OUT_FIG / "phase8_differential.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\npermutation test: observed r={perm['observed_r']:.3f}, "
          f"null median={perm['null_median']:.3f} "
          f"[{perm['null_ci'][0]:.3f}, {perm['null_ci'][1]:.3f}], "
          f"p={perm['p_one_sided']:.4f}")
    print(f"cells significantly different: {int(sig.sum())}/400")
    print("\ntop differential cells (antibody - PPI):")
    for c in cells[:10]:
        print(f"   {c['pair']:5s} delta={c['delta']:+.2f} "
              f"[{c['ci'][0]:+.2f},{c['ci'][1]:+.2f}] "
              f"ab={c['antibody']:+.2f} ppi={c['ppi']:+.2f} "
              f"{'*' if c['significant'] else ''}")
    print("\nasymmetry:")
    for a in C.ARMS:
        d = asym[a]
        print(f"   {a:9s} {d['asymmetry']:.3f} "
              f"[{d['bootstrap_ci'][0]:.3f}, {d['bootstrap_ci'][1]:.3f}]")
    print(f"\ndynamic range: SD_antibody={dyn['sd_antibody']:.3f}, "
          f"SD_ppi_full={dyn['sd_ppi_full']:.3f}, "
          f"SD_ppi_subsampled={dyn['sd_ppi_subsampled_to_antibody_units_median']:.3f} "
          f"{dyn['sd_ppi_subsampled_ci']}")
    print(f"   -> {dyn['verdict']}")
