"""Addition 7 - three analyses added in response to peer review.

Each addresses a specific gap identified in review, where a claim in the manuscript
was asserted at a lower standard of evidence than the rest of the study applies.

  1. EQUIVALENCE TESTING for the null results. The manuscript states that salt-bridge
     frequency, charge pairing, side-chain readout fraction and shape complementarity
     are "identical" or "show no difference" between the two interface classes. Those
     statements rest on non-significant p values, which is absence of evidence rather
     than evidence of absence. Two one-sided tests (TOST) invert the logic: the null
     becomes "the difference is at least as large as delta" and rejecting it licenses
     the claim of practical equivalence. Equivalence margins are set from the study's
     own effect-size distribution rather than chosen arbitrarily.

  2. PERMUTATION NULL FOR CHANNEL INDEPENDENCE. The manuscript reports that antibody
     side-chain and main-chain recognition maps correlate at 0.28 while the
     protein-protein equivalents correlate at 0.71, and builds the mechanistic
     argument on that contrast. Unlike the cross-group grammar test, this comparison
     had no null. Here the moiety label of each contact is permuted within interface,
     preserving every other property, and the two maps rebuilt. This gives the
     correlation expected if there were no channel-specific grammar at all.

  3. CRYSTALLISABILITY BIAS ON INTERFACE GEOMETRY. The antibody corpus is entirely
     X-ray, and well-diffracting antigens are plausibly smaller and more rigid than
     the general run of proteins. That could inflate the reported compactness and
     flatness of antibody interfaces without any antibody-specific cause. The test
     stratifies both groups by the size of the bound protein and by resolution, and
     asks whether the geometry difference survives within matched bins.

Outputs
  metadata/analysis/review_response.json
  figures/add7_equivalence.png
  figures/add7_channel_null.png
  figures/add7_crystallisability.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

N_PERM = 1000
IDX = {a: i for i, a in enumerate(C.AAS)}


# ------------------------------------------------------------ 1. equivalence tests
def bootstrap_diff_ci(a, b, n_boot=5000, seed=C.SEED):
    """Percentile CI on the difference of means between two unit-level samples."""
    rng = np.random.default_rng(seed)
    d = np.empty(n_boot)
    for i in range(n_boot):
        d[i] = (rng.choice(a, size=len(a), replace=True).mean()
                - rng.choice(b, size=len(b), replace=True).mean())
    return float(np.mean(a) - np.mean(b)), [float(np.percentile(d, 2.5)),
                                            float(np.percentile(d, 97.5))]


def tost(a, b, margin):
    """Two one-sided tests. Returns the larger of the two p values.

    H0 for the upper test is (mean_a - mean_b) >= +margin, and for the lower test
    (mean_a - mean_b) <= -margin. Rejecting BOTH establishes that the true difference
    lies inside (-margin, +margin), which is what a claim of equivalence requires.
    Implemented with Welch t statistics on unit-level means.
    """
    from scipy.stats import t as tdist
    na, nb = len(a), len(b)
    ma, mb = np.mean(a), np.mean(b)
    va, vb = np.var(a, ddof=1) / na, np.var(b, ddof=1) / nb
    se = np.sqrt(va + vb)
    if se == 0:
        return {"p_tost": 1.0, "note": "zero variance"}
    df = (va + vb) ** 2 / (va ** 2 / (na - 1) + vb ** 2 / (nb - 1))
    t_low = (ma - mb + margin) / se        # H0: diff <= -margin
    t_up = (ma - mb - margin) / se         # H0: diff >= +margin
    p_low = 1 - tdist.cdf(t_low, df)
    p_up = tdist.cdf(t_up, df)
    return {"p_tost": float(max(p_low, p_up)), "margin": float(margin),
            "observed_diff": float(ma - mb), "se": float(se), "df": float(df)}


def equivalence_block():
    """TOST for every quantity the manuscript describes as showing no difference."""
    import importlib
    p2 = importlib.import_module("phase2_inventory")
    recs = p2.per_interface("antibody") + p2.per_interface("ppi")
    ab = [r for r in recs if r["arm"] == "antibody"]
    pp = [r for r in recs if r["arm"] == "ppi"]

    # margin = 20% of the LARGEST effect the study reports on the same scale, so the
    # threshold for "negligible" is set by this corpus rather than by convention
    ref_field = "f_apolar_apolar"
    ref_a = p2.unit_level_values(ab, ref_field)
    ref_b = p2.unit_level_values(pp, ref_field)
    ref_diff = abs(np.mean(ref_a) - np.mean(ref_b))

    out = {}
    for field, label in (("f_salt", "Salt-bridge fraction"),
                         ("f_opp_charge", "Opposite-charge pair fraction"),
                         ("f_target_sidechain_read", "Side-chain readout fraction"),
                         ("f_stack", "Aromatic stacking fraction")):
        a = p2.unit_level_values(ab, field)
        b = p2.unit_level_values(pp, field)
        if len(a) < 10 or len(b) < 10:
            continue
        margin = 0.20 * ref_diff
        diff, ci = bootstrap_diff_ci(a, b)
        u, p_diff = mannwhitneyu(a, b, alternative="two-sided")
        res = tost(a, b, margin)
        out[field] = {
            "label": label,
            "median_antibody": float(np.median(a)),
            "median_ppi": float(np.median(b)),
            "difference_of_means": diff,
            "bootstrap_ci_on_difference": ci,
            "difference_test_p": float(p_diff),
            "equivalence_margin": float(margin),
            "tost_p": res["p_tost"],
            "equivalent_at_margin": bool(res["p_tost"] < 0.05),
        }
    out["_margin_note"] = (
        f"Equivalence margin is 20 percent of the largest composition difference the "
        f"study measures on the same scale, the apolar-apolar contact fraction "
        f"(absolute difference {ref_diff:.4f}), giving a margin of "
        f"{0.20 * ref_diff:.4f}. A significant TOST licenses the statement that the "
        f"two classes do not differ by an amount comparable to the effects the study "
        f"reports elsewhere.")
    return out


# --------------------------------------------------- 2. channel-independence null
def channel_maps(rows, symmetrise, moiety_of):
    """Two 20x20 maps, apportioning each contact by the given moiety function."""
    side = np.zeros((20, 20))
    main = np.zeros((20, 20))
    for r in rows:
        s, m = moiety_of(r)
        tot = s + m
        if tot <= 0:
            continue
        i, j = IDX[r["aa_probe"]], IDX[r["aa_target"]]
        w = r["weight"]
        side[i, j] += w * s / tot
        main[i, j] += w * m / tot
    if symmetrise:
        side = (side + side.T) / 2
        main = (main + main.T) / 2
    return C.log2_odds(side), C.log2_odds(main)


def channel_null(arm, n_perm=200, seed=C.SEED):
    """Observed side-vs-main map correlation against a moiety-permuted null.

    The permutation shuffles which contacts are side-chain mediated and which are
    main-chain mediated, within each interface, preserving interface identity,
    residue composition and the overall side/main ratio. If there is no
    channel-specific grammar, the observed correlation should match this null.
    """
    rng = np.random.default_rng(seed)
    rows = C.dedup_rows(arm, C.load_contacts(arm))
    sym = (arm == "ppi")

    def real(r):
        return (float(r.get("ag_side_contacts") or 0),
                float(r.get("ag_main_contacts") or 0))

    a, b = channel_maps(rows, sym, real)
    observed = C.matrix_similarity(a, b)

    by_iface = defaultdict(list)
    for r in rows:
        by_iface[(r["pdb_id"], r.get("direction", ""))].append(r)

    null = np.empty(n_perm)
    for k in range(n_perm):
        shuffled = {}
        for key, rs in by_iface.items():
            pairs = [real(r) for r in rs]
            order = rng.permutation(len(pairs))
            for r, idx in zip(rs, order):
                shuffled[id(r)] = pairs[idx]
        a2, b2 = channel_maps(rows, sym, lambda r: shuffled[id(r)])
        null[k] = C.matrix_similarity(a2, b2)
    p = (np.sum(null <= observed) + 1) / (n_perm + 1)
    return {"observed": float(observed),
            "null_median": float(np.median(null)),
            "null_ci": [float(np.percentile(null, 2.5)),
                        float(np.percentile(null, 97.5))],
            "p_one_sided_observed_below_null": float(p),
            "n_perm": n_perm,
            "interpretation": (
                "p is the probability of a side-vs-main correlation this LOW if the "
                "moiety of each contact carried no grammar information. A low p means "
                "the two channels encode genuinely different residue preferences.")}


# ------------------------------------------- 3. crystallisability bias on geometry
def crystallisability_block():
    """Does the compactness and flatness difference survive size and resolution matching?"""
    geo = list(csv.DictReader(
        (C.OUT_INT / "interface_geometry.csv").open()))
    man = {(r["group"], r["pdb_id"]): r
           for r in csv.DictReader((C.MET / "corpus_manifest.csv").open())}

    recs = []
    for g in geo:
        key = ("antibody" if g["arm"] == "antibody" else "ppi", g["pdb_id"])
        m = man.get(key)
        if not m:
            continue
        try:
            res = float(m["resolution"])
        except (ValueError, KeyError):
            continue
        # size of the protein being READ, which is the quantity crystallisability
        # would bias: well-diffracting antigens skew small
        try:
            length = float(m["side_b_length"] or 0)
        except ValueError:
            length = 0.0
        if length <= 0:
            continue
        for f in ("probe_rg", "probe_planarity", "target_rg", "target_planarity"):
            if not g.get(f):
                break
        else:
            recs.append({"arm": g["arm"], "resolution": res, "length": length,
                         **{f: float(g[f]) for f in ("probe_rg", "probe_planarity",
                                                     "target_rg",
                                                     "target_planarity")}})

    out = {"n_with_metadata": len(recs), "strata": {}}
    ab_len = np.array([r["length"] for r in recs if r["arm"] == "antibody"])
    pp_len = np.array([r["length"] for r in recs if r["arm"] == "ppi"])
    u, p = mannwhitneyu(ab_len, pp_len, alternative="two-sided")
    out["target_length"] = {
        "median_antibody": float(np.median(ab_len)),
        "median_ppi": float(np.median(pp_len)),
        "effect_rank_biserial": float(2 * u / (len(ab_len) * len(pp_len)) - 1),
        "p": float(p),
        "note": ("If antibody antigens were systematically smaller, the compactness "
                 "difference could follow from protein size rather than from "
                 "antibody chemistry. This tests the premise of the objection.")}

    # length-matched strata, quartiles of the pooled distribution
    edges = np.percentile(np.concatenate([ab_len, pp_len]), [0, 25, 50, 75, 100])
    for i in range(4):
        lo, hi = edges[i], edges[i + 1]
        sub = [r for r in recs if lo <= r["length"] <= hi]
        a = [r for r in sub if r["arm"] == "antibody"]
        b = [r for r in sub if r["arm"] == "ppi"]
        if len(a) < 15 or len(b) < 15:
            continue
        block = {"n_antibody": len(a), "n_ppi": len(b),
                 "length_range": [float(lo), float(hi)]}
        for f in ("probe_rg", "probe_planarity"):
            va = np.array([r[f] for r in a])
            vb = np.array([r[f] for r in b])
            u_, p_ = mannwhitneyu(va, vb, alternative="two-sided")
            block[f] = {"median_antibody": float(np.median(va)),
                        "median_ppi": float(np.median(vb)),
                        "effect_rank_biserial": float(
                            2 * u_ / (len(va) * len(vb)) - 1),
                        "p": float(p_)}
        out["strata"][f"length_{int(lo)}_{int(hi)}"] = block

    # high-resolution-only stratum, both arms
    a = [r for r in recs if r["arm"] == "antibody" and r["resolution"] < 2.0]
    b = [r for r in recs if r["arm"] == "ppi" and r["resolution"] < 2.0]
    if len(a) > 15 and len(b) > 15:
        block = {"n_antibody": len(a), "n_ppi": len(b)}
        for f in ("probe_rg", "probe_planarity"):
            va = np.array([r[f] for r in a]); vb = np.array([r[f] for r in b])
            u_, p_ = mannwhitneyu(va, vb, alternative="two-sided")
            block[f] = {"median_antibody": float(np.median(va)),
                        "median_ppi": float(np.median(vb)),
                        "effect_rank_biserial": float(
                            2 * u_ / (len(va) * len(vb)) - 1),
                        "p": float(p_)}
        out["strata"]["resolution_better_than_2A"] = block
    return out


if __name__ == "__main__":
    print("1. equivalence tests", flush=True)
    eq = equivalence_block()
    print("2. channel-independence null", flush=True)
    ch = {arm: channel_null(arm) for arm in C.ARMS}
    print("3. crystallisability bias", flush=True)
    cb = crystallisability_block()

    out = {"equivalence": eq, "channel_null": ch, "crystallisability": cb}
    (C.OUT_MET / "review_response.json").write_text(json.dumps(out, indent=1))

    # ---------------------------------------------------------------- figures
    fig, ax = plt.subplots(figsize=(9, 4.6))
    keys = [k for k in eq if not k.startswith("_")]
    ys = np.arange(len(keys))
    margin = eq[keys[0]]["equivalence_margin"]
    ax.axvspan(-margin, margin, color="0.88",
               label=f"equivalence region (+/-{margin:.4f})")
    for k, y in zip(keys, ys):
        d = eq[k]
        ax.plot(d["bootstrap_ci_on_difference"], [y, y], c="0.35", lw=2.5)
        ax.plot(d["difference_of_means"], y, "o", ms=9,
                c="#1b7837" if d["equivalent_at_margin"] else "#b2182b")
    ax.axvline(0, c="k", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([eq[k]["label"] for k in keys], fontsize=9)
    ax.set_xlabel("difference of means, antibody minus protein-protein")
    ax.set_title("Equivalence testing of the reported null results\n"
                 "green = statistically equivalent (TOST p < 0.05)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "add7_equivalence.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, arm in zip(axes, C.ARMS):
        d = ch[arm]
        ax.axvspan(d["null_ci"][0], d["null_ci"][1], color="0.85",
                   label="null: moiety carries no grammar")
        ax.axvline(d["null_median"], c="0.4", ls=":", lw=1.5)
        ax.axvline(d["observed"], c="#b2182b", lw=2.5,
                   label=f"observed r = {d['observed']:.3f}")
        ax.set_xlabel("side-chain vs main-chain map correlation")
        ax.set_title(f"{C.ARM_LABEL[arm]}\np = {d['p_one_sided_observed_below_null']:.4f}",
                     fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.set_yticks([])
    fig.suptitle("Are the two readout channels genuinely encoding different "
                 "preferences?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(C.OUT_FIG / "add7_channel_null.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    strata = [k for k in cb["strata"]]
    if strata:
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
        for ax, f, lab in zip(axes, ("probe_rg", "probe_planarity"),
                              ("probe patch radius of gyration",
                               "probe patch planarity")):
            xs, es = [], []
            for k in strata:
                b = cb["strata"][k].get(f)
                if b:
                    xs.append(k.replace("length_", "").replace("_", "-")
                              .replace("resolution-better-than-2A", "res<2A"))
                    es.append(b["effect_rank_biserial"])
            ax.bar(range(len(xs)), es, color="#2166ac")
            ax.axhline(0, c="k", lw=0.8)
            ax.set_xticks(range(len(xs)))
            ax.set_xticklabels(xs, fontsize=8, rotation=20, ha="right")
            ax.set_ylabel("effect (antibody vs PPI)")
            ax.set_title(lab, fontsize=10)
        fig.suptitle("Does the geometry difference survive matching on the size of "
                     "the bound protein and on resolution?", fontsize=11.5)
        fig.tight_layout(rect=[0, 0, 1, 0.9])
        fig.savefig(C.OUT_FIG / "add7_crystallisability.png", dpi=200,
                    bbox_inches="tight")
        plt.close(fig)

    print("\n--- equivalence ---")
    for k in keys:
        d = eq[k]
        print(f"{d['label']:34s} diff={d['difference_of_means']:+.5f} "
              f"CI=[{d['bootstrap_ci_on_difference'][0]:+.5f},"
              f"{d['bootstrap_ci_on_difference'][1]:+.5f}] "
              f"TOST p={d['tost_p']:.4f} "
              f"{'EQUIVALENT' if d['equivalent_at_margin'] else 'not established'}")
    print("\n--- channel null ---")
    for arm in C.ARMS:
        d = ch[arm]
        print(f"{arm:9s} observed={d['observed']:.3f} null={d['null_median']:.3f} "
              f"[{d['null_ci'][0]:.3f},{d['null_ci'][1]:.3f}] "
              f"p={d['p_one_sided_observed_below_null']:.4f}")
    print("\n--- crystallisability ---")
    t = cb["target_length"]
    print(f"target protein length: antibody {t['median_antibody']:.0f} vs ppi "
          f"{t['median_ppi']:.0f} effect {t['effect_rank_biserial']:+.3f} "
          f"p={t['p']:.2e}")
    for k, b in cb["strata"].items():
        rg = b.get("probe_rg", {})
        pl = b.get("probe_planarity", {})
        print(f"  {k:28s} n={b['n_antibody']:4d}/{b['n_ppi']:4d}  "
              f"rg_eff={rg.get('effect_rank_biserial', float('nan')):+.3f}  "
              f"planarity_eff={pl.get('effect_rank_biserial', float('nan')):+.3f}")
