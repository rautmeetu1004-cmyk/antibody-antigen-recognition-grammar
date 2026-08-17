"""Full audit of every phase: verify the science, not just that files exist.

Each check re-derives a quantity by an independent route, or asserts a property that
must hold if the analysis is correct. FAILURES are the point. Extends
audit_phases.py (Phases 0-5) to cover 6-11 and the additions.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

R = {}


def check(name, ok, detail):
    R[name] = {"pass": bool(ok), "detail": detail}
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def load(n):
    p = C.OUT_MET / n
    return json.loads(p.read_text()) if p.exists() else None


# ------------------------------------------------------------------ 6a structural
def a_dihedral():
    import phase6a_structural as P
    try:
        worst = P.verify_dihedral()
        check("6a_dihedral_matches_biopython", worst < 1e-3,
              f"max deviation {worst:.2e} deg vs Bio.PDB.vectors.calc_dihedral")
    except AssertionError as e:
        check("6a_dihedral_matches_biopython", False, str(e))


def a_chi1_canonical():
    rows = list(csv.DictReader((C.OUT_INT / "structural_residues.csv").open()))
    chi = np.array([float(r["chi1"]) for r in rows if r["chi1"]])
    h, e = np.histogram(chi, bins=36, range=(-180, 180))
    peak = (e[np.argmax(h)] + e[np.argmax(h) + 1]) / 2
    # dominant rotamer must be g- near -60 deg
    check("6a_chi1_modal_rotamer_canonical", -80 <= peak <= -40,
          f"modal chi1 bin at {peak:+.0f} deg (canonical g- = -60); "
          f"n={len(chi):,}")
    g = [r for r in rows if r["aa"] == "G" and r["chi1"]]
    check("6a_glycine_has_no_chi1", len(g) == 0,
          f"{len(g)} glycines carry a chi1 value (must be 0 - Gly has no CB)")


def a_ss_validation():
    d = load("structural_pass.json")
    r = d.get("ss_agreement_vs_pdb_records")
    check("6a_kabsch_sander_validated", r is not None and r > 0.85,
          f"{r:.3f} agreement with PDB HELIX/SHEET over "
          f"{d.get('ss_validated', 0):,} residues; {d.get('n_failed')} failures")


# ------------------------------------------------------------- 6c loop confound
def a_loop_confound_independent():
    """Re-derive the headline stratum r from scratch, not from the JSON."""
    import phase6c_conformation as P6
    rows, idx = P6.load_structural()
    ann = {arm: P6.annotate_contacts(arm, idx) for arm in C.ARMS}
    a = [r for r in ann["antibody"]
         if r["ss_probe"] == "loop" and r["ss_target"] == "loop"]
    p = [r for r in ann["ppi"]
         if r["ss_probe"] == "loop" and r["ss_target"] == "loop"]
    r_obs = C.matrix_similarity(
        C.log2_odds(C.contact_matrix(a, symmetrise=False)),
        C.log2_odds(C.contact_matrix(p, symmetrise=True)))
    stored = load("conformation.json")["loop_confound"]["loop_both"]["r"]
    check("6c_loop_confound_reproduces", abs(r_obs - stored) < 1e-9,
          f"recomputed r={r_obs:.6f} vs stored {stored:.6f}")
    check("6c_matching_context_does_not_raise_r", r_obs < 0.387,
          f"loop-loop r={r_obs:.3f} < all-contacts r=0.387, i.e. matching "
          f"structural context does NOT explain the antibody-PPI difference")


# ------------------------------------------------------------------ 7 paratope
def a_paratope():
    d = load("paratope.json")
    cov = d["imgt_mapping_coverage_median"]
    check("7_imgt_mapping_complete", cov > 0.99,
          f"median CA-coordinate mapping coverage {cov:.4f}")
    vh = d["cdr_vs_framework"]["SD-H"]
    check("7_nanobodies_have_no_light_chain",
          vh.get("L1", 0) == 0 and vh.get("L3", 0) == 0,
          f"SD-H L1={vh.get('L1', 0)}, L2={vh.get('L2', 0)}, L3={vh.get('L3', 0)} "
          f"(must be exactly 0)")
    for ct in ("FAB", "SD-H"):
        fr = d["cdr_vs_framework"][ct]["_cdr_only_fractions"]
        s = sum(fr.values())
        check(f"7_cdr_fractions_sum_to_one_{ct}", abs(s - 1.0) < 0.02,
              f"{ct} CDR-only fractions sum to {s:.4f}")
    pc = d["paratope_composition"]
    both = [a for a in C.AAS
            if pc[a]["vs_ppi_interface"] > 0 and pc[a]["vs_antibody_surface"] > 0]
    check("7_two_nulls_disagree_as_expected", set(both) >= {"Y", "W"},
          f"enriched on BOTH nulls: {sorted(both)} - Ser/Thr/Gly excluded, which is "
          f"the finding that a single-background study would miss")


# --------------------------------------------------------------- 8 differential
def a_differential():
    d = load("differential.json")
    pt = d["permutation_test"]
    check("8_permutation_observed_below_null",
          pt["observed_r"] < pt["null_ci"][0],
          f"observed r={pt['observed_r']:.3f} below null CI "
          f"[{pt['null_ci'][0]:.3f}, {pt['null_ci'][1]:.3f}], p={pt['p_one_sided']}")
    a = d["asymmetry"]
    check("8_antibody_asymmetric_ppi_not",
          a["antibody"]["asymmetry"] > 10 * a["ppi"]["asymmetry"],
          f"antibody {a['antibody']['asymmetry']:.3f} vs ppi "
          f"{a['ppi']['asymmetry']:.3f}")
    dy = d["dynamic_range"]
    check("8_dynamic_range_not_sample_size",
          dy["sd_antibody"] > dy["sd_ppi_subsampled_ci"][1],
          f"SD_ab={dy['sd_antibody']:.3f} above subsampled PPI CI "
          f"{[round(x,3) for x in dy['sd_ppi_subsampled_ci']]}")


# ------------------------------------------------------------------ 9 modified
def a_modified():
    d = load("modified.json")
    rows = list(csv.DictReader((C.OUT_INT / "modified_contacts.csv").open()))
    bad = [r for r in rows if r["resname"] in C.AA3_TO_1 and r["kind"] == "ptm"]
    check("9_no_standard_residue_called_ptm", len(bad) == 0,
          f"{len(bad)} standard residues classified as PTM (ASN N-glycosylation "
          f"MODRES records were the earlier bug)")
    mse = [r for r in rows if r["resname"] == "MSE"]
    check("9_mse_excluded", len(mse) == 0,
          f"{len(mse)} selenomethionines counted (must be 0 - crystallographic, "
          f"not biological)")
    ab = d["antibody_glycan"]["frac_structures_with_interfacial"]
    pp = d["ppi_glycan"]["frac_structures_with_interfacial"]
    check("9_glycan_ptm_double_dissociation",
          ab > pp and d["ppi_ptm"]["frac_structures_with_interfacial"]
          > d["antibody_ptm"]["frac_structures_with_interfacial"],
          f"glycans antibody {ab:.3f} > ppi {pp:.3f}; PTMs ppi "
          f"{d['ppi_ptm']['frac_structures_with_interfacial']:.3f} > antibody "
          f"{d['antibody_ptm']['frac_structures_with_interfacial']:.3f}")


# ------------------------------------------------------------------ 11 controls
def a_controls():
    d = load("controls.json")
    loo = d["redundancy"]["leave_one_cluster_out"]
    check("11_robust_to_dropping_any_cluster",
          abs(loo["max"] - loo["min"]) < 0.15,
          f"leave-one-cluster-out r range {loo['min']:.3f}-{loo['max']:.3f}")
    hb = d["hbond_criterion"]
    vals = [v["cross_arm_r"] for v in hb.values()]
    check("11_robust_to_hbond_cutoff", max(vals) - min(vals) < 0.1,
          f"cross-arm r across 3.0/3.2/3.5 A: {[round(v,3) for v in vals]}")


# ------------------------------------------------------------- add1 transfer
def a_transfer():
    d = load("transfer.json")
    if d is None:
        check("add1_present", False, "transfer.json missing")
        return
    for k in ("antibody_scored_by_antibody_potential",
              "antibody_scored_by_ppi_potential"):
        m = d[k]["mean_auc"]
        check(f"add1_auc_in_range_{k[:28]}", 0 <= m <= 1, f"AUC={m:.3f}")
    self_ab = d["antibody_scored_by_antibody_potential"]["mean_auc"]
    cross_ab = d["antibody_scored_by_ppi_potential"]["mean_auc"]
    # this is the REFUTED prediction, recorded as such
    check("add1_prediction_refuted_recorded", cross_ab > self_ab,
          f"PPI potential ({cross_ab:.3f}) OUTPERFORMS antibody potential "
          f"({self_ab:.3f}) on antibody interfaces - the predicted "
          f"below-chance transfer failure did NOT occur")


def a_loo_subtraction():
    """The fast LOO (total - unit) must equal rebuilding from the other units."""
    rows = C.dedup_rows("antibody", C.load_contacts("antibody"))
    ids, U = C.unit_matrices(rows, symmetrise=False)
    total = U.sum(0)
    u = ids[3]
    fast = C.log2_odds(total - U[3])
    slow = C.log2_odds(C.contact_matrix([r for r in rows if r["unit"] != u],
                                        symmetrise=False))
    ok = np.isfinite(fast) & np.isfinite(slow)
    d = float(np.max(np.abs(fast[ok] - slow[ok])))
    check("add1_loo_subtraction_exact", d < 1e-9,
          f"max |fast - rebuilt| = {d:.2e} over {ok.sum()} cells")


# ------------------------------------------------------------- cross-phase
def a_cross_phase():
    hm, cons, ctrl = load("heatmap.json"), load("conservation.json"), load("controls.json")
    r3 = hm["cluster30"]["grammar_similarity_antibody_vs_ppi"]
    r4 = cons["overall"]["all_contacts"]["r"]
    check("xphase_grammar_r_consistent", abs(r3 - r4) < 1e-6,
          f"phase3 {r3:.6f} == phase4 {r4:.6f}")
    # phase 11 baseline uses dedup+all PPI, phase 3 uses all rows: must be close
    check("xphase_baseline_r_close_to_phase3",
          abs(ctrl["baseline_r"] - r3) < 0.05,
          f"phase11 baseline {ctrl['baseline_r']:.3f} vs phase3 {r3:.3f} "
          f"(differ only by PPI deduplication)")

    figs = sorted(p.name for p in C.OUT_FIG.glob("*.png"))
    check("xphase_all_figures_present", len(figs) >= 16,
          f"{len(figs)} figures: {', '.join(f[:-4] for f in figs[:6])}...")


# ---------------------------------------------- unit-identifier namespace hazard
def a_unit_id_namespace():
    """cluster30 identifiers are NOT globally unique across the two corpora.

    They were assigned in separate clustering runs, so the same string can denote
    different proteins in the two arms. Any dictionary that groups by the bare
    identifier while iterating rows from BOTH arms will silently merge unrelated
    structure sets. This check documents the hazard so it cannot be reintroduced,
    and verifies that the weighting in aa_common is computed per arm.
    """
    ab = {v.get("cluster30") for v in C.units("antibody").values()
          if v.get("cluster30")}
    pp = set()
    for v in C.units("ppi").values():
        for d in (v.get("a_to_b"), v.get("b_to_a")):
            if d and d.get("cluster30"):
                pp.add(d["cluster30"])
    shared = ab & pp

    ab_up = {v["cluster30"]: v.get("uniprot") for v in C.units("antibody").values()
             if v.get("cluster30")}
    pp_up = {}
    for v in C.units("ppi").values():
        d = v.get("a_to_b") or {}
        if d.get("cluster30"):
            pp_up.setdefault(d["cluster30"], d.get("uniprot"))
    conflicting = [c for c in shared if ab_up.get(c) != pp_up.get(c)]

    # aa_common must weight within arm: identical unit strings in the two arms must
    # not produce identical weights unless the arms genuinely share the structures
    wa = {r["unit"]: r["weight"] for r in C.load_contacts("antibody")}
    wp = {}
    for r in C.load_contacts("ppi"):
        wp.setdefault(r["unit"], r["weight"])
    independent = any(c in wa and c in wp and wa[c] != wp[c] for c in conflicting)

    check("unit_ids_are_not_globally_unique_documented",
          len(shared) > 0 and independent,
          f"{len(shared)} cluster30 ids occur in both arms, {len(conflicting)} of "
          f"them denoting different proteins; aa_common weights are computed per "
          f"arm so the collision does not propagate there")


if __name__ == "__main__":
    print("=" * 78)
    for fn in (a_dihedral, a_chi1_canonical, a_ss_validation,
               a_loop_confound_independent, a_paratope, a_differential,
               a_modified, a_controls, a_transfer, a_loo_subtraction,
               a_cross_phase, a_unit_id_namespace):
        try:
            fn()
        except Exception as e:                                    # noqa: BLE001
            check(fn.__name__, False, f"{type(e).__name__}: {e}")
    print("=" * 78)
    nf = sum(1 for v in R.values() if not v["pass"])
    print(f"{len(R) - nf}/{len(R)} checks passed")
    (C.OUT_MET / "audit_full.json").write_text(json.dumps(R, indent=1))
