"""Addition 5 - epitope logic, made symmetric with the paratope analysis.

Phase 7 scored PARATOPE composition against two nulls -- general-PPI interfaces and
the antibody's own accessible surface -- and the second null overturned the standard
Ser/Thr/Gly enrichment story. The equivalent EPITOPE question was never asked, and it
is the more interesting one:

  Is an epitope chemically distinguishable from the rest of the antigen's own
  accessible surface, or is it simply whichever patch was available?

If epitopes are NOT distinguishable from their own antigen's surface while PPI
interfaces ARE distinguishable from theirs, that is a clean structural statement of
the difference between a co-evolved interface and an imposed one: co-evolved partners
bind a surface that was shaped for binding; antibodies impose their chemistry on an
arbitrary patch.

Three measurements, each computed identically for all three interface roles
(antibody paratope, antibody epitope, PPI partner) so the comparison is symmetric:

  COMPOSITION   log2 enrichment of each residue at the interface relative to the
                SAME protein's own non-interface accessible surface.
  ACCESSIBILITY relative solvent accessibility of interface vs non-interface surface
                residues -- are epitopes on unusually exposed positions?
  PROTRUSION    distance of each residue from its chain's centroid, normalised by
                the chain's radius -- the classic "epitopes protrude" claim, which
                needs coordinates rather than SASA alone.

Outputs
  metadata/analysis/epitope_logic.json
  figures/add5_epitope_logic.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
import time
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import phase6a_structural as P

SURFACE_SASA = 20.0        # A^2 in the isolated chain to count as "surface"

# Tien et al. 2013 theoretical maxima, for relative accessibility
MAXA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
        "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
        "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}

ROLES = [("antibody", "probe", "Antibody paratope"),
         ("antibody", "target", "Antibody epitope"),
         ("ppi", "probe", "PPI partner interface")]


def load_rows():
    rows = list(csv.DictReader((C.OUT_INT / "structural_residues.csv").open()))
    for r in rows:
        r["is_interface"] = int(r["is_interface"])
        for k in ("sasa_free", "sasa_complex"):
            r[k] = float(r[k]) if r[k] else np.nan
    return rows


def protrusion():
    """{(arm,pdb,chain,resi) -> normalised distance from that chain's centroid}."""
    out = {}
    t0 = time.time()
    plan = []
    for arm in C.ARMS:
        for p in sorted(P.targets(arm)):
            plan.append((arm, p))
    done = 0
    for k, (arm, pdb) in enumerate(plan, 1):
        d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
        path = next((d / nm for nm in (f"{pdb}.pdb", f"{pdb}.cif")
                     if (d / nm).exists()), None)
        if path is None:
            continue
        try:
            model = P.load_model(path)
        except Exception:                                        # noqa: BLE001
            continue
        for ch in model:
            coords, keys = [], []
            for r in ch:
                if "CA" in r and r.get_resname().strip() in C.AA3_TO_1:
                    coords.append(r["CA"].coord)
                    keys.append((arm, pdb, ch.id, P.residue_key(r)))
            if len(coords) < 10:
                continue
            xyz = np.array(coords, dtype=float)
            cen = xyz.mean(0)
            dist = np.linalg.norm(xyz - cen, axis=1)
            rad = dist.max() or 1.0
            for kk, dd in zip(keys, dist):
                out[kk] = float(dd / rad)
        done += 1
        if k % 300 == 0:
            print(f"  protrusion {k}/{len(plan)}  {time.time()-t0:4.0f}s", flush=True)
    print(f"  protrusion computed for {done} structures")
    return out


def unit_of_structure():
    """{(arm, pdb) -> cluster30 unit}, so residue counts can be unit-weighted."""
    out = {}
    for arm in C.ARMS:
        for p, v in C.units(arm).items():
            out[(arm, p)] = (v.get("cluster30") if arm == "antibody"
                             else (v.get("a_to_b") or {}).get("cluster30"))
    return out


UNIT_OF = unit_of_structure()


def composition(rows, arm, role, interface):
    """Unit-weighted residue composition.

    Weighting matters here for the same reason it does everywhere else in this
    study: 45 lysozyme complexes would otherwise contribute 45 times the evidence
    of a singleton. Each residue is weighted by the reciprocal of the number of
    structures in its cluster30 unit.
    """
    per_unit_structs = defaultdict(set)
    for r in rows:
        u = UNIT_OF.get((r["arm"], r["pdb_id"]))
        if u:
            per_unit_structs[u].add(r["pdb_id"])

    w = defaultdict(float)
    n_raw = 0
    for r in rows:
        if r["arm"] != arm or r["role"] != role:
            continue
        if r["is_interface"] != (1 if interface else 0):
            continue
        if not interface and not (r["sasa_free"] > SURFACE_SASA):
            continue                     # non-interface must still be SURFACE
        if r["aa"] not in C.AAS:
            continue
        u = UNIT_OF.get((r["arm"], r["pdb_id"]))
        if not u:
            continue
        w[r["aa"]] += 1.0 / len(per_unit_structs[u])
        n_raw += 1
    tot = sum(w.values()) or 1.0
    return {a: w.get(a, 0.0) / tot for a in C.AAS}, n_raw


if __name__ == "__main__":
    rows = load_rows()
    prot = protrusion()

    out = {"surface_definition_A2": SURFACE_SASA, "roles": {}}
    enrich = {}
    for arm, role, label in ROLES:
        iface, n_i = composition(rows, arm, role, True)
        surf, n_s = composition(rows, arm, role, False)
        e = {a: (float(np.log2(iface[a] / surf[a]))
                 if iface[a] > 0 and surf[a] > 0 else np.nan) for a in C.AAS}
        enrich[label] = e

        # Accessibility and protrusion are compared at UNIT level: each structure's
        # interface and surface medians are averaged within cluster30 unit, and the
        # test runs over units. A residue-level test here would report n in the tens
        # of thousands, give p values below 1e-60 for trivial effects, and ignore
        # the correlation between residues of the same structure.
        sub = [r for r in rows if r["arm"] == arm and r["role"] == role]
        by_struct = defaultdict(lambda: {"rsa_i": [], "rsa_s": [],
                                         "pr_i": [], "pr_s": []})
        for r in sub:
            if r["aa"] not in MAXA or not np.isfinite(r["sasa_free"]):
                continue
            rsa = r["sasa_free"] / MAXA[r["aa"]]
            pv = prot.get((r["arm"], r["pdb_id"], r["chain"], r["resi"]))
            acc = by_struct[(r["arm"], r["pdb_id"])]
            if r["is_interface"]:
                acc["rsa_i"].append(rsa)
                if pv is not None:
                    acc["pr_i"].append(pv)
            elif r["sasa_free"] > SURFACE_SASA:
                acc["rsa_s"].append(rsa)
                if pv is not None:
                    acc["pr_s"].append(pv)

        per_unit = defaultdict(lambda: {"rsa_i": [], "rsa_s": [],
                                        "pr_i": [], "pr_s": []})
        for key, rec in by_struct.items():
            u = UNIT_OF.get(key)
            if not u:
                continue
            for f in ("rsa_i", "rsa_s", "pr_i", "pr_s"):
                if rec[f]:
                    per_unit[u][f].append(float(np.median(rec[f])))

        def paired(f_i, f_s):
            a, b = [], []
            for u, e in per_unit.items():
                if e[f_i] and e[f_s]:
                    a.append(np.mean(e[f_i]))
                    b.append(np.mean(e[f_s]))
            return np.array(a), np.array(b)

        rsa_i, rsa_s = paired("rsa_i", "rsa_s")
        pr_i, pr_s = paired("pr_i", "pr_s")

        def cmp(a, b):
            """Wilcoxon signed-rank: interface and surface are PAIRED within unit,
            since both come from the same proteins."""
            from scipy.stats import wilcoxon
            if len(a) < 10 or len(b) < 10:
                return {}
            try:
                stat, p = wilcoxon(a, b)
            except ValueError:
                return {}
            d = a - b
            # matched-pairs rank-biserial: proportion of units favouring interface
            eff = float(np.mean(np.sign(d)))
            return {"median_interface": float(np.median(a)),
                    "median_surface": float(np.median(b)),
                    "effect_rank_biserial": eff,
                    "p": float(p), "n_units": int(len(a))}

        vals = np.array([v for v in e.values() if np.isfinite(v)])
        out["roles"][label] = {
            "n_interface_residues": n_i, "n_surface_residues": n_s,
            "composition_log2_vs_own_surface": {a: round(e[a], 3) for a in C.AAS},
            "mean_abs_enrichment": float(np.mean(np.abs(vals))),
            "max_abs_enrichment": float(np.max(np.abs(vals))),
            "relative_accessibility": cmp(rsa_i, rsa_s),
            "protrusion": cmp(pr_i, pr_s),
        }

    out["interpretation"] = (
        "mean_abs_enrichment is the headline: how chemically distinguishable an "
        "interface is from its OWN protein's accessible surface. A co-evolved PPI "
        "interface should be distinguishable; if an antibody epitope is much less "
        "so, antibodies are imposing chemistry on an arbitrary patch rather than "
        "selecting a chemically special one.")
    (C.OUT_MET / "epitope_logic.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.2))
    ax = axes[0]
    x = np.arange(20)
    cols = {"Antibody paratope": "#b2182b", "Antibody epitope": "#1b7837",
            "PPI partner interface": "#2166ac"}
    order = sorted(C.AAS, key=lambda a: -enrich["Antibody paratope"].get(a, 0))
    for i, (lab, col) in enumerate(cols.items()):
        ax.bar(x + (i - 1) * 0.27, [enrich[lab].get(a, np.nan) for a in order],
               width=0.26, color=col, label=lab)
    ax.axhline(0, c="0.4", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(order, fontsize=8)
    for t, a in zip(ax.get_xticklabels(), order):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    ax.set_ylabel("log$_2$ enrichment vs the SAME protein's own surface")
    ax.set_title("Is the interface chemically special\nrelative to its own surface?",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    labs = list(cols)
    vals = [out["roles"][l]["mean_abs_enrichment"] for l in labs]
    ax.bar(range(3), vals, color=[cols[l] for l in labs], width=0.55)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=9)
    ax.set_xticks(range(3))
    ax.set_xticklabels([l.replace(" ", "\n") for l in labs], fontsize=8)
    ax.set_ylabel("mean |log$_2$ enrichment| vs own surface")
    ax.set_title("How distinguishable is each interface\nfrom its own surface?",
                 fontsize=10)

    ax = axes[2]
    w = 0.35
    for i, (key, lab) in enumerate((("relative_accessibility", "rel. accessibility"),
                                    ("protrusion", "protrusion"))):
        vs = [out["roles"][l][key].get("effect_rank_biserial", np.nan) for l in labs]
        ax.bar(np.arange(3) + (i - 0.5) * w, vs, width=w, label=lab)
    ax.axhline(0, c="0.4", lw=0.8)
    ax.set_xticks(range(3))
    ax.set_xticklabels([l.replace(" ", "\n") for l in labs], fontsize=8)
    ax.set_ylabel("effect: interface vs non-interface surface")
    ax.set_title("Are interface residues more exposed\nand more protruding?",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Addition 5 - epitope and paratope logic against each protein's "
                 "own surface", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.savefig(C.OUT_FIG / "add5_epitope_logic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for lab in cols:
        d = out["roles"][lab]
        print(f"\n{lab}: n_iface={d['n_interface_residues']:,} "
              f"n_surface={d['n_surface_residues']:,}")
        print(f"   mean |log2 enrichment| vs own surface = "
              f"{d['mean_abs_enrichment']:.3f} (max {d['max_abs_enrichment']:.3f})")
        ra, pr = d["relative_accessibility"], d["protrusion"]
        if ra:
            print(f"   accessibility: iface {ra['median_interface']:.3f} vs surface "
                  f"{ra['median_surface']:.3f}  effect "
                  f"{ra['effect_rank_biserial']:+.3f} p={ra['p']:.1e}")
        if pr:
            print(f"   protrusion:    iface {pr['median_interface']:.3f} vs surface "
                  f"{pr['median_surface']:.3f}  effect "
                  f"{pr['effect_rank_biserial']:+.3f} p={pr['p']:.1e}")
