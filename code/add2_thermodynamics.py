"""Addition 2 - thermodynamics from buried surface: is antibody binding enthalpic?

Phase 2 found antibody interfaces are aromatic-rich and apolar-poor, and I inferred
from that composition that antibody binding leans enthalpic while general PPI leans
entropic (classical hydrophobic effect). That was an INFERENCE FROM COMPOSITION, not
a thermodynamic measurement. This turns it into one.

The heat-capacity change on binding is the classical discriminator between
hydrophobically-driven and polar-driven association, and it is calculable from
buried surface alone:

    dCp = 0.45 * dA_apolar - 0.26 * dA_polar        (cal/mol/K)

(Murphy & Freire 1992; the same parameterisation family as Spolar & Record 1994.)
A large negative dCp is the signature of burying apolar surface -- the hydrophobic
effect, entropy-driven. A dCp near zero, or dominated by the polar term, indicates
polar/H-bond-driven, enthalpy-leaning association.

Rigour note. The lazy way to split buried surface into apolar and polar is to apply
a per-residue-type apolar fraction to the per-residue dSASA already computed in
Phase 6a. That is an approximation, and this study has already been bitten twice by
approximations that looked fine. So this does the real thing: a fresh Shrake-Rupley
pass at ATOM level, classifying each atom by element (C and S apolar; N and O
polar), in the complex and in isolation.

Also reported: apolar/polar burial ratio, and dCp normalised per unit buried area,
because antibody interfaces are smaller and a raw dCp comparison would partly be a
size comparison.

Outputs
  interactions/analysis/thermodynamics.csv
  metadata/analysis/thermodynamics.json
  figures/add2_thermodynamics.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import random
import sys
import time
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB.SASA import ShrakeRupley

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import phase6a_structural as P

SEED = 20260816
PPI_SAMPLE = 1000
SASA_POINTS = 64
_sr = ShrakeRupley(n_points=SASA_POINTS)

APOLAR = {"C", "S", "SE"}
POLAR = {"N", "O"}

# Murphy & Freire 1992 coefficients, cal/(mol*K*A^2)
C_APOLAR = 0.45
C_POLAR = -0.26


def atom_sasa_by_polarity(model, chain_ids):
    """(apolar_A2, polar_A2) for the requested chains, atom-level Shrake-Rupley."""
    sub = P.clean_copy(model, set(chain_ids))
    if not len(list(sub)):
        return 0.0, 0.0
    _sr.compute(sub, level="A")
    ap = po = 0.0
    for ch in sub:
        for res in ch:
            for at in res:
                el = (at.element or "").strip().upper()
                if el in APOLAR:
                    ap += at.sasa
                elif el in POLAR:
                    po += at.sasa
    return ap, po


def process(arm, pdb_id, probe_chains, target_chains):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return None
    model = P.load_model(path)

    ap_c, po_c = atom_sasa_by_polarity(model, set(probe_chains) | set(target_chains))
    ap_a, po_a = atom_sasa_by_polarity(model, probe_chains)
    ap_b, po_b = atom_sasa_by_polarity(model, target_chains)
    if ap_c == 0 and po_c == 0:
        return None

    # buried on complexation = (free A + free B) - complex
    d_ap = (ap_a + ap_b) - ap_c
    d_po = (po_a + po_b) - po_c
    d_tot = d_ap + d_po
    if d_tot <= 0:
        return None

    # Murphy & Freire is written in terms of the CHANGE in accessible area, which is
    # NEGATIVE on binding. d_ap / d_po above are buried (positive) areas, so they
    # enter with a minus sign. Without it dCp comes out positive, which would invert
    # the entire thermodynamic interpretation.
    dcp = C_APOLAR * (-d_ap) + C_POLAR * (-d_po)
    return {
        "arm": arm, "pdb_id": pdb_id,
        "dA_apolar": round(d_ap, 1), "dA_polar": round(d_po, 1),
        "dA_total": round(d_tot, 1),
        "frac_apolar_buried": round(d_ap / d_tot, 4),
        "dCp_cal_mol_K": round(dcp, 1),
        "dCp_per_1000A2": round(dcp / (d_tot / 1000.0), 1),
    }


if __name__ == "__main__":
    rng = random.Random(SEED)
    tg = P.targets("antibody")
    plan = [("antibody", p, v[0], v[1]) for p, v in tg.items()]
    tp = P.targets("ppi")
    ids = sorted(tp)
    if len(ids) > PPI_SAMPLE:
        ids = sorted(rng.sample(ids, PPI_SAMPLE))
    plan += [("ppi", p, tp[p][0], tp[p][1]) for p in ids]
    print(f"{len(plan)} structures", flush=True)

    rows, n_fail = [], 0
    t0 = time.time()
    for k, (arm, pdb, pc, tc) in enumerate(plan, 1):
        try:
            r = process(arm, pdb, pc, tc)
            if r:
                rows.append(r)
            else:
                n_fail += 1
        except Exception:                                        # noqa: BLE001
            n_fail += 1
        if k % 100 == 0:
            el = time.time() - t0
            print(f"  {k}/{len(plan)}  {el:5.0f}s  "
                  f"eta {el/k*(len(plan)-k)/60:4.1f} min", flush=True)

    with (C.OUT_INT / "thermodynamics.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    # unit-weighted comparison: average within cluster30, then across units
    units = {arm: {} for arm in C.ARMS}
    for arm in C.ARMS:
        u = C.units(arm)
        for p, v in u.items():
            units[arm][p] = (v.get("cluster30") if arm == "antibody"
                             else (v.get("a_to_b") or {}).get("cluster30"))

    from scipy.stats import mannwhitneyu
    fields = ["dA_total", "dA_apolar", "dA_polar", "frac_apolar_buried",
              "dCp_cal_mol_K", "dCp_per_1000A2"]
    stats = {}
    per_unit = {arm: defaultdict(lambda: defaultdict(list)) for arm in C.ARMS}
    for r in rows:
        u = units[r["arm"]].get(r["pdb_id"])
        if not u:
            continue
        for f in fields:
            per_unit[r["arm"]][u][f].append(r[f])

    vals = {arm: {f: np.array([np.mean(d[f]) for d in per_unit[arm].values() if d[f]])
                  for f in fields} for arm in C.ARMS}
    for f in fields:
        a, b = vals["antibody"][f], vals["ppi"][f]
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        stats[f] = {"median_antibody": float(np.median(a)),
                    "median_ppi": float(np.median(b)),
                    "effect_rank_biserial": float(2 * u / (len(a) * len(b)) - 1),
                    "p": float(p),
                    "n_units_antibody": int(len(a)), "n_units_ppi": int(len(b))}

    out = {"n_structures": len(rows), "n_failed": n_fail, "stats": stats,
           "model": "dCp = 0.45*dA_apolar - 0.26*dA_polar (cal/mol/K), "
                    "Murphy & Freire 1992",
           "method_note": (
               "atom-level Shrake-Rupley (64 points), elements C/S/Se counted apolar "
               "and N/O polar; buried area = (free A + free B) - complex. No "
               "residue-composition approximation is used. dCp_per_1000A2 is the "
               "size-normalised comparison, because antibody interfaces are smaller "
               "and a raw dCp difference would partly restate that."),
           "interpretation": (
               "A more negative dCp indicates greater apolar burial and a more "
               "hydrophobically (entropically) driven association. If antibody "
               "interfaces are enthalpy-leaning as Phase 2 inferred, their "
               "size-normalised dCp should be LESS negative than general PPI.")}
    (C.OUT_MET / "thermodynamics.json").write_text(json.dumps(out, indent=1))

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    for ax, f, lab in zip(axes,
                          ("frac_apolar_buried", "dCp_per_1000A2", "dA_total"),
                          ("fraction of buried surface that is apolar",
                           "$\\Delta C_p$ per 1000 $\\AA^2$ buried (cal/mol/K)",
                           "total buried surface ($\\AA^2$)")):
        for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
            v = vals[arm][f]
            ax.hist(v, bins=35, density=True, alpha=0.55, color=col,
                    label=C.ARM_LABEL[arm])
            ax.axvline(np.median(v), color=col, ls="--", lw=1.3)
        ax.set_xlabel(lab, fontsize=9)
        s = stats[f]
        ax.set_title(f"effect {s['effect_rank_biserial']:+.2f}, p={s['p']:.1e}",
                     fontsize=9)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Addition 2 - burial thermodynamics: is antibody binding "
                 "enthalpy-leaning?", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(C.OUT_FIG / "add2_thermodynamics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\n{len(rows)} structures ({n_fail} failed)")
    for f in fields:
        s = stats[f]
        print(f"{f:22s} antibody {s['median_antibody']:10.3f}  ppi "
              f"{s['median_ppi']:10.3f}  effect {s['effect_rank_biserial']:+.3f}  "
              f"p={s['p']:.2e}")
