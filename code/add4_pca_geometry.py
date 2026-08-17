"""Addition 4 - interface geometry and the PCA view: are antibody interfaces separable?

Two related questions the study has not yet answered visually.

  GEOMETRIC SPREAD. Phase 7 found epitopes average four sequence segments, but
  sequence discontinuity says nothing about SPATIAL arrangement: four segments can
  fold into one compact patch or scatter across a face. Measured here as the radius
  of gyration of the interface patch, its planarity and elongation (from the
  eigenvalues of the patch coordinate covariance), and the number of spatially
  separate sub-patches. Also relevant to the nanobody-binds-clefts claim: a concave
  cleft binder should present a less planar, more compact paratope than a Fab.

  PCA / SEPARABILITY. Every result so far is a single statistic comparing two
  distributions. This asks the aggregate question directly: given ~18 measured
  interface properties, do antibody interfaces occupy their own region of descriptor
  space, or do they sit inside the general-PPI cloud? A classifier AUC quantifies
  what the scatter shows. Cross-validation is GROUPED BY cluster30 unit, so no unit
  appears in both train and test -- without that, redundancy inflates the AUC.

Outputs
  interactions/analysis/interface_geometry.csv
  metadata/analysis/pca_geometry.json
  figures/add4_pca.png
  figures/add4_geometry.png
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
PPI_SAMPLE = 1000
PATCH_LINK = 8.0        # CA-CA distance defining one contiguous spatial patch


def interface_coords(arm, pdb_id, iface_keys):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return None
    model = P.load_model(path)
    out = {"probe": [], "target": []}
    for ch in model:
        for r in ch:
            if "CA" not in r:
                continue
            key = (ch.id, P.residue_key(r))
            side = iface_keys.get(key)
            if side:
                out[side].append(r["CA"].coord)
    return {k: np.array(v, dtype=float) for k, v in out.items() if len(v) >= 3}


def patch_stats(xyz):
    """Rg, planarity, elongation and spatial sub-patch count for a set of CA points."""
    if xyz is None or len(xyz) < 3:
        return None
    cen = xyz.mean(0)
    d = xyz - cen
    rg = float(np.sqrt((d ** 2).sum(1).mean()))
    ev = np.linalg.eigvalsh(np.cov(d.T))
    ev = np.sort(ev)[::-1]
    tot = ev.sum() or 1.0
    # planarity: how little variance lies along the smallest axis
    planarity = float(1 - ev[2] / tot * 3)
    elongation = float(ev[0] / max(ev[1], 1e-9))
    # spatial sub-patches: single-linkage components at PATCH_LINK
    t = cKDTree(xyz)
    pairs = t.query_pairs(PATCH_LINK)
    parent = list(range(len(xyz)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    n_patch = len({find(i) for i in range(len(xyz))})
    return {"rg": round(rg, 2), "planarity": round(planarity, 4),
            "elongation": round(elongation, 3), "n_spatial_patches": n_patch,
            "n_residues": len(xyz)}


def build_geometry():
    rng = random.Random(SEED)
    plan = []
    for arm in C.ARMS:
        tg = P.targets(arm)
        ids = sorted(tg)
        if arm == "ppi" and len(ids) > PPI_SAMPLE:
            ids = sorted(rng.sample(ids, PPI_SAMPLE))
        for p in ids:
            pc, tc, iface = tg[p]
            keys = {}
            for (ch, resi) in iface:
                keys[(ch, resi)] = "probe" if ch in pc else "target"
            plan.append((arm, p, keys))

    rows, t0 = [], time.time()
    for k, (arm, pdb, keys) in enumerate(plan, 1):
        try:
            co = interface_coords(arm, pdb, keys)
        except Exception:                                        # noqa: BLE001
            co = None
        if co:
            rec = {"arm": arm, "pdb_id": pdb}
            for side in ("probe", "target"):
                st = patch_stats(co.get(side))
                if st:
                    for kk, vv in st.items():
                        rec[f"{side}_{kk}"] = vv
            if len(rec) > 2:
                rows.append(rec)
        if k % 200 == 0:
            el = time.time() - t0
            print(f"  geometry {k}/{len(plan)}  {el:4.0f}s", flush=True)
    return rows


def load_descriptors(geom_rows):
    """Join interface inventory + geometry into one descriptor table."""
    inv = {}
    for r in csv.DictReader((C.OUT_INT / "interface_inventory.csv").open()):
        if r["arm"] == "ppi" and r["direction"] != "a_to_b":
            continue
        inv[(r["arm"], r["pdb_id"])] = r
    units = {}
    for arm in C.ARMS:
        u = C.units(arm)
        for p, v in u.items():
            units[(arm, p)] = (v.get("cluster30") if arm == "antibody"
                               else (v.get("a_to_b") or {}).get("cluster30"))

    inv_fields = ["f_hbond", "f_salt", "f_vdw", "f_aromatic", "f_apolar_apolar",
                  "f_cation_pi", "f_opp_charge", "f_like_charge",
                  "f_target_sidechain_read", "contacts_per_probe_residue",
                  "mean_min_dist", "n_pairs", "n_probe_residues"]
    geo_fields = ["probe_rg", "probe_planarity", "probe_elongation",
                  "probe_n_spatial_patches", "target_rg", "target_planarity",
                  "target_n_spatial_patches"]

    X, y, groups, names = [], [], [], inv_fields + geo_fields
    for g in geom_rows:
        key = (g["arm"], g["pdb_id"])
        i = inv.get(key)
        u = units.get(key)
        if i is None or u is None:
            continue
        try:
            row = [float(i[f]) for f in inv_fields] + [float(g[f]) for f in geo_fields]
        except (KeyError, ValueError):
            continue
        if not all(np.isfinite(row)):
            continue
        X.append(row)
        y.append(0 if g["arm"] == "antibody" else 1)
        groups.append(u)
    return np.array(X), np.array(y), np.array(groups), names


if __name__ == "__main__":
    geom = build_geometry()
    with (C.OUT_INT / "interface_geometry.csv").open("w", newline="") as fh:
        keys = sorted({k for r in geom for k in r})
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(geom)
    print(f"geometry computed for {len(geom)} interfaces")

    X, y, groups, names = load_descriptors(geom)
    print(f"descriptor matrix: {X.shape}, antibody={int((y==0).sum())}, "
          f"ppi={int((y==1).sum())}")

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score

    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=5).fit(Xs)
    Z = pca.transform(Xs)

    clf = LogisticRegression(max_iter=3000, C=1.0)
    cv = GroupKFold(n_splits=5)
    pred = cross_val_predict(clf, Xs, y, cv=cv, groups=groups,
                             method="predict_proba")[:, 1]
    auc = float(roc_auc_score(y, pred))

    clf.fit(Xs, y)
    coefs = sorted(zip(names, clf.coef_[0]), key=lambda kv: -abs(kv[1]))

    # overlap: what fraction of antibody points fall inside the PPI 95% ellipse of PC1/2
    ppi_mu = Z[y == 1, :2].mean(0)
    ppi_cov = np.cov(Z[y == 1, :2].T)
    inv = np.linalg.inv(ppi_cov)
    md = np.einsum("ij,jk,ik->i", Z[y == 0, :2] - ppi_mu, inv, Z[y == 0, :2] - ppi_mu)
    frac_inside = float(np.mean(md < 5.991))          # chi2 2df 95%

    out = {
        "n_interfaces": int(len(y)),
        "explained_variance_ratio": [round(float(v), 4)
                                     for v in pca.explained_variance_ratio_],
        "classifier_auc_group_cv": round(auc, 4),
        "frac_antibody_inside_ppi_95pct_ellipse": round(frac_inside, 4),
        "top_discriminative_descriptors": [{"descriptor": n, "coef": round(float(c), 3)}
                                           for n, c in coefs[:10]],
        "geometry_medians": {},
        "note": ("Cross-validation is grouped by cluster30 so no unit appears in "
                 "both train and test; an ungrouped CV would report an inflated AUC "
                 "driven by redundancy. AUC near 1 would mean antibody interfaces "
                 "are a separable class; AUC near 0.5 would mean they are "
                 "indistinguishable in these descriptors."),
    }
    for f in ("probe_rg", "probe_planarity", "probe_n_spatial_patches",
              "target_rg", "target_planarity", "target_n_spatial_patches"):
        vals = {arm: np.array([g[f] for g in geom
                               if g["arm"] == arm and f in g], dtype=float)
                for arm in C.ARMS}
        from scipy.stats import mannwhitneyu
        try:
            u_, p_ = mannwhitneyu(vals["antibody"], vals["ppi"],
                                  alternative="two-sided")
            eff = 2 * u_ / (len(vals["antibody"]) * len(vals["ppi"])) - 1
        except ValueError:
            p_, eff = float("nan"), float("nan")
        out["geometry_medians"][f] = {
            "antibody": float(np.median(vals["antibody"])),
            "ppi": float(np.median(vals["ppi"])),
            "effect_rank_biserial": float(eff), "p": float(p_)}
    (C.OUT_MET / "pca_geometry.json").write_text(json.dumps(out, indent=1))

    # ---------------------------------------------------------------- PCA figure
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
    ax = axes[0]
    ax.scatter(Z[y == 1, 0], Z[y == 1, 1], s=7, alpha=0.35, color="#2166ac",
               label=f"Protein-protein (n={(y==1).sum()})", edgecolors="none")
    ax.scatter(Z[y == 0, 0], Z[y == 0, 1], s=9, alpha=0.55, color="#b2182b",
               label=f"Antibody-antigen (n={(y==0).sum()})", edgecolors="none")
    for lab, m, col in (("PPI", Z[y == 1, :2].mean(0), "#2166ac"),
                        ("Ab", Z[y == 0, :2].mean(0), "#b2182b")):
        ax.plot(m[0], m[1], marker="X", ms=15, color=col, mec="k", mew=1.2)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.0%} variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.0%} variance)")
    ax.set_title(f"Interface descriptor space\n"
                 f"{frac_inside:.0%} of antibody interfaces fall inside the PPI 95% "
                 f"ellipse\ngrouped-CV classifier AUC = {auc:.3f}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    top = coefs[:10][::-1]
    ax.barh(range(len(top)), [c for _, c in top],
            color=["#b2182b" if c < 0 else "#2166ac" for _, c in top])
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([n for n, _ in top], fontsize=8)
    ax.axvline(0, c="0.4", lw=0.8)
    ax.set_xlabel("logistic coefficient\n(negative = antibody-like, "
                  "positive = PPI-like)")
    ax.set_title("What separates them", fontsize=10)
    fig.suptitle("Addition 4 - are antibody interfaces a separable class?",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(C.OUT_FIG / "add4_pca.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ----------------------------------------------------------- geometry figure
    gf = ["probe_rg", "probe_planarity", "probe_n_spatial_patches",
          "target_rg", "target_planarity", "target_n_spatial_patches"]
    labs = ["paratope radius of gyration ($\\AA$)", "paratope planarity",
            "paratope spatial patches", "epitope radius of gyration ($\\AA$)",
            "epitope planarity", "epitope spatial patches"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.6))
    for ax, f, lab in zip(axes.ravel(), gf, labs):
        for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
            v = np.array([g[f] for g in geom if g["arm"] == arm and f in g],
                         dtype=float)
            if f.endswith("patches"):
                ax.hist(v, bins=np.arange(0.5, 12.5), density=True, alpha=0.55,
                        color=col, label=C.ARM_LABEL[arm])
            else:
                ax.hist(v, bins=35, density=True, alpha=0.55, color=col,
                        label=C.ARM_LABEL[arm])
            ax.axvline(np.median(v), color=col, ls="--", lw=1.2)
        s = out["geometry_medians"][f]
        ax.set_xlabel(lab, fontsize=9)
        ax.set_title(f"effect {s['effect_rank_biserial']:+.2f}, p={s['p']:.1e}",
                     fontsize=9)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("Addition 4 - geometric spread of paratope and epitope",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(C.OUT_FIG / "add4_geometry.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\nPCA: PC1 {pca.explained_variance_ratio_[0]:.1%}, "
          f"PC2 {pca.explained_variance_ratio_[1]:.1%}")
    print(f"grouped-CV classifier AUC = {auc:.3f}")
    print(f"{frac_inside:.1%} of antibody interfaces inside PPI 95% ellipse")
    print("top discriminative descriptors:")
    for n, c in coefs[:8]:
        print(f"   {n:32s} {c:+.3f}")
    print("\ngeometry:")
    for f in gf:
        s = out["geometry_medians"][f]
        print(f"   {f:26s} ab={s['antibody']:7.2f} ppi={s['ppi']:7.2f} "
              f"eff={s['effect_rank_biserial']:+.3f} p={s['p']:.2e}")
