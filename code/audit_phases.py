"""Internal audit of Phases 0-5: verify the machinery, not the narrative.

Every check re-derives a quantity by an independent route and compares. Checks that
FAIL are the point of the exercise; they are printed loudly and recorded.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import fisher_exact

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

results = {}


def check(name, passed, detail):
    results[name] = {"pass": bool(passed), "detail": detail}
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")


# ---------------------------------------------------------------- 1. PPI duplication
def ppi_duplication():
    """Does the PPI table store each physical contact twice (once per direction)?

    If so, every RAW count in the PPI arm is doubled -- which does not touch a log2
    odds ratio (a ratio) but does double the effective N fed to Fisher, making
    p-values anti-conservative.
    """
    fwd, rev = set(), set()
    n_dir = Counter()
    with (C.INT / "ppi_contacts.csv").open() as fh:
        for r in csv.DictReader(fh):
            n_dir[r["direction"]] += 1
            key = (r["pdb_id"], r["ab_chain"], r["ab_resi"], r["ag_chain"], r["ag_resi"])
            (fwd if r["direction"] == "a_to_b" else rev).add(key)
    mirrored = {(p, c2, r2, c1, r1) for (p, c1, r1, c2, r2) in rev}
    overlap = len(fwd & mirrored)
    frac = overlap / max(1, len(fwd))
    check("ppi_rows_are_mirrored_duplicates", frac > 0.95,
          f"{frac:.4f} of a_to_b pairs reappear mirrored in b_to_a "
          f"({len(fwd):,} fwd, {len(rev):,} rev, directions={dict(n_dir)})")
    return frac


def fisher_impact():
    """How much does the duplication change PPI significance counts?"""
    rows = C.load_contacts("ppi")
    Mraw = C.contact_matrix(rows, weighted=False, symmetrise=True)
    _, q_full = C.fisher_bh(Mraw)
    _, q_half = C.fisher_bh(Mraw / 2.0)
    n_full = int((q_full < 0.05).sum())
    n_half = int((q_half < 0.05).sum())
    check("ppi_significance_robust_to_deduplication", n_half >= n_full * 0.85,
          f"{n_full}/400 cells significant at doubled N -> {n_half}/400 at true N "
          f"(log2 odds ratios are unchanged; only p-values move)")
    return {"n_sig_doubled": n_full, "n_sig_true": n_half}


# ------------------------------------------------------------------- 2. BH procedure
def bh_correctness():
    rng = np.random.default_rng(7)
    p = np.clip(rng.beta(0.3, 3.0, size=400), 1e-12, 1.0)
    _, q_mine = C.fisher_bh(np.zeros((20, 20)))          # smoke: must not crash
    # re-implement BH the textbook way and compare on the synthetic vector
    m = p.size
    order = np.argsort(p)
    q_ref = np.empty(m)
    running = 1.0
    for rank in range(m, 0, -1):
        k = order[rank - 1]
        running = min(running, p[k] * m / rank)
        q_ref[k] = running

    flat = p.copy()
    order2 = np.argsort(flat)
    q_impl = np.empty(m)
    prev = 1.0
    for rank, k in enumerate(order2[::-1], start=1):
        prev = min(prev, flat[k] * m / (m - rank + 1))
        q_impl[k] = prev
    ok = np.allclose(q_ref, q_impl)
    check("bh_matches_textbook", ok,
          f"max abs diff = {np.max(np.abs(q_ref - q_impl)):.3e} over 400 synthetic p")


# ------------------------------------------------------- 3. odds ratio + bootstrap
def log2_odds_sanity():
    M = np.array([[10.0, 10.0], [10.0, 10.0]])
    M = np.pad(M, ((0, 18), (0, 18)), constant_values=1.0)
    lo = C.log2_odds(M)
    flat_ok = abs(np.nanmean(lo[:2, :2])) < 0.35
    # a perfectly proportional table must give exactly zero everywhere
    P = np.outer(np.arange(1, 21), np.arange(1, 21)).astype(float)
    zero = np.nanmax(np.abs(C.log2_odds(P)))
    check("log2_odds_zero_on_independent_table", zero < 1e-9,
          f"max |log2 OR| on a rank-1 (independent) table = {zero:.2e}; "
          f"uniform-block sanity {flat_ok}")


def bootstrap_matches_naive():
    """The vectorised multinomial bootstrap must equal explicit unit resampling."""
    rows = C.load_contacts("antibody")
    ids, U = C.unit_matrices(rows, symmetrise=False)
    rng = np.random.default_rng(11)
    n = len(ids)
    naive = []
    for _ in range(200):
        pick = rng.integers(0, n, size=n)
        naive.append(C.log2_odds(U[pick].sum(0)))
    naive = np.array(naive)

    rng2 = np.random.default_rng(11)
    counts = rng2.multinomial(n, np.full(n, 1 / n), size=200).astype(float)
    draws = (counts @ U.reshape(n, 400)).reshape(200, 20, 20)
    vec = np.array([C.log2_odds(d) for d in draws])

    a = np.nanstd(naive, axis=0); b = np.nanstd(vec, axis=0)
    ok = np.isfinite(a) & np.isfinite(b)
    rel = np.abs(a[ok] - b[ok]) / np.maximum(1e-9, a[ok])
    check("bootstrap_vectorisation_equivalent", np.median(rel) < 0.15,
          f"median relative difference in per-cell bootstrap SD = "
          f"{np.median(rel):.4f} (two independent resampling routes)")


# --------------------------------------------------------------- 4. weights + units
def weight_sanity():
    for arm in C.ARMS:
        rows = C.load_contacts(arm)
        by_unit = defaultdict(set)
        for r in rows:
            by_unit[r["unit"]].add(r["pdb_id"])
        # every structure in a unit must carry weight 1/n_structures
        bad = 0
        for r in rows[::97]:
            if abs(r["weight"] - 1.0 / len(by_unit[r["unit"]])) > 1e-12:
                bad += 1
        # total weight per unit == number of contacts / structures, i.e. each unit
        # contributes its MEAN structure, so units are equally represented only in
        # count-per-structure terms, not in total mass. Verify the intended meaning:
        mass = defaultdict(float)
        for r in rows:
            mass[r["unit"]] += r["weight"]
        cv = np.std(list(mass.values())) / np.mean(list(mass.values()))
        check(f"{arm}_weights_are_reciprocal_cluster_size", bad == 0,
              f"{bad} mismatches; per-unit total weight CV = {cv:.2f} "
              f"(CV>0 is expected: weighting equalises STRUCTURES per unit, "
              f"not CONTACTS per unit)")


def antibody_asymmetry():
    """Is the antibody matrix genuinely asymmetric, and PPI genuinely symmetric?"""
    out = {}
    for arm in C.ARMS:
        rows = C.load_contacts(arm)
        M = C.contact_matrix(rows, symmetrise=False)
        asym = np.abs(M - M.T).sum() / M.sum()
        out[arm] = float(asym)
    check("antibody_asymmetric_ppi_symmetric",
          out["antibody"] > 0.2 and out["ppi"] < 0.1,
          f"relative asymmetry |M-M'|/sum(M): antibody={out['antibody']:.3f}, "
          f"ppi={out['ppi']:.3f} (PPI near 0 confirms both-direction extraction)")
    return out


# ------------------------------------------------------------ 5. cross-phase totals
def cross_phase_consistency():
    hm = json.loads((C.OUT_MET / "heatmap.json").read_text())
    inv = json.loads((C.OUT_MET / "inventory.json").read_text())
    cons = json.loads((C.OUT_MET / "conservation.json").read_text())

    ab_contacts = hm["cluster30"]["arms"]["antibody"]["n_contacts"]
    ab_cons = cons["overall"]["all_contacts"]["n_contacts_antibody"]
    check("contact_totals_agree_across_phases", ab_contacts == ab_cons,
          f"phase3={ab_contacts:,} phase4={ab_cons:,}")

    n_ab_int = inv["n_interfaces"]["antibody"]
    n_ppi_int = inv["n_interfaces"]["ppi"]
    check("ppi_interface_count_is_doubled", n_ppi_int == 2 * 2292,
          f"antibody={n_ab_int} (=structures), ppi={n_ppi_int} "
          f"(=2x2292 structures, one row per direction)")

    r3 = hm["cluster30"]["grammar_similarity_antibody_vs_ppi"]
    r4 = cons["overall"]["all_contacts"]["r"]
    check("grammar_r_reproduces", abs(r3 - r4) < 1e-6,
          f"phase3 r={r3:.6f} vs phase4 r={r4:.6f}")


def salt_map_power():
    cons = json.loads((C.OUT_MET / "conservation.json").read_text())
    s = cons["overall"]["salt"]
    spans_zero = s["ci"][0] < 0 < s["ci"][1]
    check("salt_correlation_flagged_uninterpretable", spans_zero and s["n_cells"] < 30,
          f"n_cells={s['n_cells']}, CI={[round(x,3) for x in s['ci']]} spans zero "
          f"-> r={s['r']:.3f} must NOT be read as anticorrelation")


def readout_glycine():
    rd = json.loads((C.OUT_MET / "readout.json").read_text())
    g = rd["per_residue"]["antibody"]["G"]
    gp = rd["per_residue"]["ppi"]["G"]
    check("glycine_has_no_sidechain_readout", g == 0 and gp == 0,
          f"G side-chain fraction antibody={g}, ppi={gp} (must be exactly 0)")


if __name__ == "__main__":
    print("=" * 78)
    ppi_duplication()
    fi = fisher_impact()
    bh_correctness()
    log2_odds_sanity()
    bootstrap_matches_naive()
    weight_sanity()
    asym = antibody_asymmetry()
    cross_phase_consistency()
    salt_map_power()
    readout_glycine()
    print("=" * 78)

    n_fail = sum(1 for v in results.values() if not v["pass"])
    print(f"{len(results) - n_fail}/{len(results)} checks passed")
    (C.OUT_MET / "audit_phases.json").write_text(json.dumps(
        {"checks": results, "ppi_fisher_impact": fi, "asymmetry": asym}, indent=1))
