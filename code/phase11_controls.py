"""Phase 11 - controls: could any of this be an artefact?

Four ways the headline result (antibody and PPI grammars differ; permutation
p = 0.001) could be wrong, each tested rather than argued:

  1. RESOLUTION. Poorly resolved structures have noisier contacts, and noise
     attenuates a correlation. The two arms are well matched overall (median 1.83 vs
     1.80 A), but matching medians is not matching distributions. The comparison is
     re-run within resolution strata.
  2. REDUNDANCY. 45 lysozyme and 39 spike complexes sit in the antibody corpus.
     Unit weighting is supposed to absorb them; this drops the largest clusters
     outright and re-measures, plus leave-one-cluster-out over the top clusters.
  3. H-BOND CRITERION. The 3.5 A cutoff inherited from the extraction is a choice.
     The H-bond map is re-derived from the ATOM tables at 3.0, 3.2 and 3.5 A using
     only N/O donor-acceptor pairs, and the maps are compared.
  4. CONSTRUCT TYPE. 203 of 512 antibodies are nanobodies with no light chain and a
     longer CDR-H3. If Fab and nanobody grammars disagree with each other as much as
     antibodies disagree with PPI, then "antibody grammar" is not one thing.

Outputs
  figures/phase11_controls.png
  metadata/analysis/controls.json
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C


def resolutions():
    out = {}
    for r in csv.DictReader((C.MET / "corpus_manifest.csv").open()):
        try:
            out[(r["group"], r["pdb_id"])] = float(r["resolution"])
        except (ValueError, KeyError):
            pass
    return out


def ctypes():
    return {r["pdb_id"]: r["antibody_type"]
            for r in csv.DictReader((C.MET / "corpus_manifest.csv").open())
            if r["group"] == "antibody"}


def grammar_r(rows_ab, rows_ppi):
    a = C.log2_odds(C.contact_matrix(rows_ab, symmetrise=False))
    b = C.log2_odds(C.contact_matrix(rows_ppi, symmetrise=True))
    return C.matrix_similarity(a, b)


def size_matched_null(ab_all, ppi_all, n_a, n_p, n_draw=40, seed=C.SEED):
    """Restricting shrinks the data; small samples attenuate r on their own."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_draw):
        ia = rng.choice(len(ab_all), size=min(n_a, len(ab_all)), replace=False)
        ip = rng.choice(len(ppi_all), size=min(n_p, len(ppi_all)), replace=False)
        vals.append(grammar_r([ab_all[i] for i in ia], [ppi_all[i] for i in ip]))
    v = np.array(vals)
    return float(np.median(v)), [float(np.percentile(v, 2.5)),
                                 float(np.percentile(v, 97.5))]


def hbond_maps_from_atoms(arm, cutoffs=(3.0, 3.2, 3.5)):
    """Re-derive H-bond maps straight from the atom table at several cutoffs.

    A hydrogen bond is approximated as any N/O-to-N/O heavy-atom pair within the
    cutoff. This is deliberately cruder than the inherited flag; the question is
    whether the MAP is stable to the criterion, not whether the criterion is ideal.
    """
    units = C.units(arm)
    idx = {a: i for i, a in enumerate(C.AAS)}
    mats = {c: np.zeros((20, 20)) for c in cutoffs}
    structs = defaultdict(set)
    recs = []
    with (C.INT / f"{arm}_atom_contacts.csv").open() as fh:
        for row in csv.DictReader(fh):
            if arm == "ppi" and row.get("direction") != "a_to_b":
                continue
            if row["ab_atom"][:1] not in ("N", "O") or row["ag_atom"][:1] not in ("N", "O"):
                continue
            u = units.get(row["pdb_id"])
            if u is None:
                continue
            if arm == "ppi":
                u = u.get(row["direction"])
                if u is None:
                    continue
            unit = u.get("cluster30")
            if not unit or unit == "NA":
                continue
            a1 = C.AA3_TO_1.get(row["ab_aa3"]); a2 = C.AA3_TO_1.get(row["ag_aa3"])
            if a1 is None or a2 is None:
                continue
            structs[unit].add(row["pdb_id"])
            recs.append((unit, idx[a1], idx[a2], float(row["dist"])))
    for unit, i, j, d in recs:
        w = 1.0 / len(structs[unit])
        for c in cutoffs:
            if d <= c:
                mats[c][i, j] += w
    if arm == "ppi":
        mats = {c: (m + m.T) / 2 for c, m in mats.items()}
    return {c: C.log2_odds(m) for c, m in mats.items()}


if __name__ == "__main__":
    res = resolutions()
    ct = ctypes()
    ab = C.dedup_rows("antibody", C.load_contacts("antibody"))
    pp = C.dedup_rows("ppi", C.load_contacts("ppi"))
    out = {"baseline_r": grammar_r(ab, pp)}
    print(f"baseline grammar r = {out['baseline_r']:.3f}")

    # ---------------------------------------------------------- 1. resolution
    strata = {"high (<2.0 A)": lambda v: v is not None and v < 2.0,
              "low (>=2.0 A)": lambda v: v is not None and v >= 2.0}
    out["resolution"] = {}
    for name, pred in strata.items():
        a = [r for r in ab if pred(res.get(("antibody", r["pdb_id"])))]
        p = [r for r in pp if pred(res.get(("ppi", r["pdb_id"])))]
        if len(a) < 500 or len(p) < 500:
            continue
        r_obs = grammar_r(a, p)
        med, ci = size_matched_null(ab, pp, len(a), len(p))
        out["resolution"][name] = {
            "r": r_obs, "n_antibody": len(a), "n_ppi": len(p),
            "null_median": med, "null_ci": ci,
            "verdict": "consistent with baseline" if ci[0] <= r_obs <= ci[1]
                       else "differs from size-matched expectation"}
        print(f"  resolution {name:16s} r={r_obs:+.3f} null={med:+.3f} "
              f"[{ci[0]:+.3f},{ci[1]:+.3f}]")

    # ---------------------------------------------------------- 2. redundancy
    units_ab = C.units("antibody")
    cl = Counter(units_ab[p].get("cluster30") for p in {r["pdb_id"] for r in ab}
                 if p in units_ab)
    top = [c for c, _ in cl.most_common(5)]
    out["redundancy"] = {"top_clusters": cl.most_common(5), "drop_top_n": {}}
    for n in (1, 2, 3, 5):
        drop = set(top[:n])
        a = [r for r in ab
             if units_ab.get(r["pdb_id"], {}).get("cluster30") not in drop]
        out["redundancy"]["drop_top_n"][str(n)] = {
            "r": grammar_r(a, pp), "n_contacts": len(a),
            "dropped_clusters": list(drop)}
    print("  drop-top-N clusters:",
          {k: round(v["r"], 3) for k, v in out["redundancy"]["drop_top_n"].items()})

    # leave-one-cluster-out over the 20 largest
    loo = []
    for c, _ in cl.most_common(20):
        a = [r for r in ab if units_ab.get(r["pdb_id"], {}).get("cluster30") != c]
        loo.append(grammar_r(a, pp))
    out["redundancy"]["leave_one_cluster_out"] = {
        "n": len(loo), "min": float(np.min(loo)), "max": float(np.max(loo)),
        "median": float(np.median(loo))}
    print(f"  leave-one-cluster-out r range: {np.min(loo):.3f} - {np.max(loo):.3f}")

    # ------------------------------------------------------- 3. H-bond criterion
    hb = {arm: hbond_maps_from_atoms(arm) for arm in C.ARMS}
    out["hbond_criterion"] = {}
    cuts = sorted(hb["antibody"])
    for c in cuts:
        out["hbond_criterion"][str(c)] = {
            "cross_arm_r": C.matrix_similarity(hb["antibody"][c], hb["ppi"][c]),
            "antibody_self_r_vs_3.5": C.matrix_similarity(hb["antibody"][c],
                                                          hb["antibody"][3.5]),
            "ppi_self_r_vs_3.5": C.matrix_similarity(hb["ppi"][c], hb["ppi"][3.5])}
    print("  H-bond cutoff robustness:",
          {k: round(v["cross_arm_r"], 3)
           for k, v in out["hbond_criterion"].items()})

    # ----------------------------------------------------- 4. construct type
    fab = [r for r in ab if ct.get(r["pdb_id"]) == "FAB"]
    vhh = [r for r in ab if ct.get(r["pdb_id"]) == "SD-H"]
    m_fab = C.log2_odds(C.contact_matrix(fab, symmetrise=False))
    m_vhh = C.log2_odds(C.contact_matrix(vhh, symmetrise=False))
    out["construct_type"] = {
        "n_contacts_fab": len(fab), "n_contacts_nanobody": len(vhh),
        "fab_vs_nanobody_r": C.matrix_similarity(m_fab, m_vhh),
        "fab_vs_ppi_r": grammar_r(fab, pp),
        "nanobody_vs_ppi_r": grammar_r(vhh, pp),
        "interpretation": (
            "if fab_vs_nanobody_r is much HIGHER than either arm's correlation with "
            "PPI, the two antibody formats share one grammar and 'antibody grammar' "
            "is a coherent object")}
    print(f"  FAB vs nanobody r = {out['construct_type']['fab_vs_nanobody_r']:.3f}; "
          f"FAB vs PPI = {out['construct_type']['fab_vs_ppi_r']:.3f}; "
          f"nanobody vs PPI = {out['construct_type']['nanobody_vs_ppi_r']:.3f}")

    (C.OUT_MET / "controls.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))

    ax = axes[0]
    names = list(out["resolution"])
    if names:
        for i, n in enumerate(names):
            d = out["resolution"][n]
            ax.plot(d["null_ci"], [i, i], c="0.6", lw=6, solid_capstyle="butt")
            ax.plot(d["r"], i, "o", ms=9, c="#b2182b")
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
    ax.axvline(out["baseline_r"], ls="--", c="0.4", lw=1)
    ax.set_xlabel("grammar r"); ax.set_title("Resolution deconfound", fontsize=10)

    ax = axes[1]
    ks = list(out["redundancy"]["drop_top_n"])
    ax.plot([0] + [int(k) for k in ks],
            [out["baseline_r"]] + [out["redundancy"]["drop_top_n"][k]["r"] for k in ks],
            "o-", c="#b2182b")
    ax.set_xlabel("largest antigen clusters dropped")
    ax.set_ylabel("grammar r")
    ax.set_title("Redundancy control", fontsize=10)

    ax = axes[2]
    ax.plot(cuts, [out["hbond_criterion"][str(c)]["cross_arm_r"] for c in cuts],
            "o-", c="#b2182b", label="antibody vs PPI")
    ax.set_xlabel("H-bond distance cutoff ($\\AA$)")
    ax.set_ylabel("grammar r")
    ax.set_title("H-bond criterion robustness", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[3]
    d = out["construct_type"]
    vals = [d["fab_vs_nanobody_r"], d["fab_vs_ppi_r"], d["nanobody_vs_ppi_r"]]
    ax.bar(range(3), vals, color=["#1b7837", "#b2182b", "#2166ac"], width=0.55)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.012, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Fab vs\nnanobody", "Fab vs\nPPI", "nanobody\nvs PPI"],
                       fontsize=8)
    ax.set_ylabel("grammar r")
    ax.set_title("Is 'antibody grammar' one thing?", fontsize=10)

    fig.suptitle("Phase 11 - controls and robustness", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(C.OUT_FIG / "phase11_controls.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("\nwrote controls.json + phase11_controls.png")
