"""Addition 3 - shape complementarity (Sc), after Lawrence & Colman 1993.

Required rather than optional: Lawrence & Colman DEFINED Sc while comparing
antibody-antigen interfaces against protease-inhibitor and oligomeric interfaces,
which is exactly the comparison this study makes. A 2026 paper on the same question
that does not report Sc will be sent back.

It is also a genuine test here, not a formality. Phase 2 found antibody interfaces
are smaller (59 vs 87 residue pairs) yet more tightly packed (closest approach 3.587
vs 3.641 A). Those pull Sc in opposite directions, so the result is not predetermined.

Method and its honest limitation. Lawrence & Colman compute Sc on the Connolly
MOLECULAR surface (probe-rolled, re-entrant surfaces included). No Connolly surface
generator is available here, so this uses the solvent-ACCESSIBLE surface: points on
spheres of radius (vdW + probe), keeping those not occluded by the same side's other
atoms, and taking the outward normal from the atom centre. For each interface point
on A, the nearest interface point on B gives

    S(x) = -(n_A . n_B) * exp(-w * d^2),      w = 0.5 A^-2

and Sc is the median of S over both directions. Because the underlying surface
differs from Connolly's, ABSOLUTE VALUES ARE NOT DIRECTLY COMPARABLE to the numbers
in Lawrence & Colman; the antibody-vs-PPI comparison, computed identically for both
arms, is what this supports. That limitation is stated rather than glossed.

Outputs
  interactions/analysis/shape_complementarity.csv
  metadata/analysis/shape_complementarity.json
  figures/add3_shape_complementarity.png
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
from scipy.spatial import cKDTree

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import phase6a_structural as P

SEED = 20260816
PPI_SAMPLE = 600            # Sc is the most expensive analysis; sample is reported
PROBE = 1.4
W = 0.5
N_SPHERE = 120
INTERFACE_CUT = 8.0         # atoms within this of the partner define the interface

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90, "P": 1.80}


def sphere_points(n=N_SPHERE):
    """Golden-spiral points on the unit sphere."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


UNIT = sphere_points()


def side_atoms(model, chains):
    xyz, rad = [], []
    for ch in model:
        if ch.id not in chains:
            continue
        for res in ch:
            if res.get_resname().strip() not in C.AA3_TO_1:
                continue
            for at in res:
                el = (at.element or "C").strip().upper()
                if el == "H":
                    continue
                xyz.append(at.coord)
                rad.append(VDW.get(el, 1.7) + PROBE)
    if not xyz:
        return np.zeros((0, 3)), np.zeros(0)
    return np.array(xyz, dtype=float), np.array(rad, dtype=float)


def surface_points(xyz, rad, near_mask):
    """Accessible surface points (and normals) for atoms flagged by near_mask."""
    if xyz.shape[0] == 0 or not near_mask.any():
        return np.zeros((0, 3)), np.zeros((0, 3))
    tree = cKDTree(xyz)
    pts, nrm = [], []
    idxs = np.where(near_mask)[0]
    for i in idxs:
        # Accessibility is decided on the EXPANDED (vdW+probe) sphere -- that is
        # where a probe centre may sit -- but the surface point itself is reported
        # on the vdW sphere. That is the contact portion of the Connolly molecular
        # surface. Reporting the expanded point instead (plain SAS) leaves the two
        # partners' surfaces several angstroms apart, and exp(-w d^2) then drives
        # every Sc toward zero, which is exactly what happened first time.
        p = xyz[i] + rad[i] * UNIT
        p_vdw = xyz[i] + (rad[i] - PROBE) * UNIT
        # occluded if inside any OTHER atom's expanded sphere on the same side
        cand = tree.query_ball_point(xyz[i], rad[i] + rad.max())
        cand = [j for j in cand if j != i]
        if cand:
            d = np.linalg.norm(p[:, None, :] - xyz[cand][None, :, :], axis=2)
            keep = np.all(d >= rad[cand][None, :] - 1e-6, axis=1)
        else:
            keep = np.ones(len(p), dtype=bool)
        if keep.any():
            pts.append(p_vdw[keep])
            nrm.append(UNIT[keep])
    if not pts:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.vstack(pts), np.vstack(nrm)


def sc_for(model, ca, cb):
    xa, ra = side_atoms(model, set(ca))
    xb, rb = side_atoms(model, set(cb))
    if xa.shape[0] == 0 or xb.shape[0] == 0:
        return None
    ta, tb = cKDTree(xa), cKDTree(xb)
    near_a = np.array([len(tb.query_ball_point(p, INTERFACE_CUT)) > 0 for p in xa])
    near_b = np.array([len(ta.query_ball_point(p, INTERFACE_CUT)) > 0 for p in xb])
    if not near_a.any() or not near_b.any():
        return None

    pa, na = surface_points(xa, ra, near_a)
    pb, nb = surface_points(xb, rb, near_b)
    if pa.shape[0] < 50 or pb.shape[0] < 50:
        return None

    # THE BURIED PATCH: a point counts only if the PARTNER occludes it, i.e. it lies
    # inside some partner atom's expanded sphere. Taking every accessible point of an
    # interface atom instead (including those facing solvent, away from the partner)
    # mixes in near-perpendicular normals and collapses the median toward zero.
    def buried(pts, xyz_other, rad_other):
        t = cKDTree(xyz_other)
        keep = np.zeros(len(pts), dtype=bool)
        rmax = rad_other.max()
        for k, q in enumerate(t.query_ball_point(pts, rmax)):
            if q:
                d = np.linalg.norm(xyz_other[q] - pts[k], axis=1)
                keep[k] = np.any(d < rad_other[q])
        return keep

    ka_sel = buried(pa, xb, rb)
    kb_sel = buried(pb, xa, ra)
    if ka_sel.sum() < 30 or kb_sel.sum() < 30:
        return None
    pa, na = pa[ka_sel], na[ka_sel]
    pb, nb = pb[kb_sel], nb[kb_sel]

    da, ia = cKDTree(pb).query(pa)
    db, ib = cKDTree(pa).query(pb)

    s_a = -np.einsum("ij,ij->i", na, nb[ia]) * np.exp(-W * da ** 2)
    s_b = -np.einsum("ij,ij->i", nb, na[ib]) * np.exp(-W * db ** 2)
    allv = np.concatenate([s_a, s_b])
    return {"sc": float(np.median(allv)),
            "sc_a_to_b": float(np.median(s_a)), "sc_b_to_a": float(np.median(s_b)),
            "n_points": int(len(allv))}


def process(arm, pdb_id, ca, cb):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return None
    model = P.load_model(path)
    r = sc_for(model, ca, cb)
    if r is None:
        return None
    return {"arm": arm, "pdb_id": pdb_id, **r}


if __name__ == "__main__":
    rng = random.Random(SEED)
    ta = P.targets("antibody")
    plan = [("antibody", p, v[0], v[1]) for p, v in ta.items()]
    tp = P.targets("ppi")
    ids = sorted(tp)
    if len(ids) > PPI_SAMPLE:
        ids = sorted(rng.sample(ids, PPI_SAMPLE))
    plan += [("ppi", p, tp[p][0], tp[p][1]) for p in ids]
    print(f"{len(plan)} structures", flush=True)

    rows, n_fail = [], 0
    t0 = time.time()
    for k, (arm, pdb, ca, cb) in enumerate(plan, 1):
        try:
            r = process(arm, pdb, ca, cb)
            rows.append(r) if r else None
            if r is None:
                n_fail += 1
        except Exception:                                        # noqa: BLE001
            n_fail += 1
        if k % 50 == 0:
            el = time.time() - t0
            print(f"  {k}/{len(plan)}  {el:5.0f}s  "
                  f"eta {el/k*(len(plan)-k)/60:5.1f} min", flush=True)

    with (C.OUT_INT / "shape_complementarity.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    units = {}
    for arm in C.ARMS:
        u = C.units(arm)
        for p, v in u.items():
            units[(arm, p)] = (v.get("cluster30") if arm == "antibody"
                               else (v.get("a_to_b") or {}).get("cluster30"))
    per = {arm: defaultdict(list) for arm in C.ARMS}
    for r in rows:
        cl = units.get((r["arm"], r["pdb_id"]))
        if cl:
            per[r["arm"]][cl].append(r["sc"])
    vals = {arm: np.array([np.mean(v) for v in per[arm].values() if v])
            for arm in C.ARMS}

    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(vals["antibody"], vals["ppi"], alternative="two-sided")
    out = {
        "n_structures": len(rows), "n_failed": n_fail,
        "ppi_sampled": PPI_SAMPLE,
        "median_sc_antibody": float(np.median(vals["antibody"])),
        "median_sc_ppi": float(np.median(vals["ppi"])),
        "n_units_antibody": int(len(vals["antibody"])),
        "n_units_ppi": int(len(vals["ppi"])),
        "effect_rank_biserial": float(
            2 * u / (len(vals["antibody"]) * len(vals["ppi"])) - 1),
        "p": float(p),
        "limitation": (
            "Solvent-accessible-surface approximation to Lawrence & Colman's "
            "Connolly-molecular-surface Sc. Absolute values are NOT comparable to "
            "their published figures; the antibody-vs-PPI contrast, computed "
            "identically for both arms, is what this supports."),
    }
    (C.OUT_MET / "shape_complementarity.json").write_text(json.dumps(out, indent=1))

    fig, ax = plt.subplots(figsize=(7.2, 5))
    for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
        ax.hist(vals[arm], bins=35, density=True, alpha=0.55, color=col,
                label=f"{C.ARM_LABEL[arm]} (median {np.median(vals[arm]):.3f})")
        ax.axvline(np.median(vals[arm]), color=col, ls="--", lw=1.3)
    ax.set_xlabel("shape complementarity $S_c$ (SAS approximation)")
    ax.set_ylabel("density")
    ax.set_title(f"Addition 3 - shape complementarity\n"
                 f"effect {out['effect_rank_biserial']:+.2f}, p={out['p']:.2e}",
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "add3_shape_complementarity.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    print(f"\nSc antibody {out['median_sc_antibody']:.3f} "
          f"(n={out['n_units_antibody']} units) vs PPI "
          f"{out['median_sc_ppi']:.3f} (n={out['n_units_ppi']}), "
          f"effect {out['effect_rank_biserial']:+.3f}, p={out['p']:.2e}")
