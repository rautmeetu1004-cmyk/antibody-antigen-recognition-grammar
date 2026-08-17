"""Phase 2 - interface inventory: what chemistry and physics govern these interfaces?

Answers the question directly and quantitatively, for both groups, so every number
about antibodies has a protein-protein baseline sitting next to it.

The forces available to a protein-protein interface are few, and each leaves a
measurable signature in the contact tables:

  electrostatic     salt bridges (charged N to carboxylate O, <=4.0 A)
  hydrogen bonding  directional donor-acceptor pairs (<=3.5 A)
  van der Waals     close-packed apolar surface (<=4.5 A)
  aromatic/pi       ring stacking (<=5.5 A, <=35 deg) and cation-pi
  hydrophobic       not a force but an entropic effect; read out as apolar-apolar
                    contact enrichment and apolar burial
  shape/entropy     contact density, side-chain vs main-chain mediation

This script quantifies all of them per interface, unit-weighted, and reports the
distributions rather than only the means -- the mean interface does not exist, and
antibody interfaces in particular are bimodal (nanobody vs Fab).

Outputs
  figures/phase2_interface_inventory.png
  figures/phase2_physicochemistry.png
  interactions/analysis/interface_inventory.csv     one row per interface
  metadata/analysis/inventory.json
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

# physicochemical classes used for the force decomposition
CHARGED_POS = set("RK")          # H is left out: mostly neutral at pH 7
CHARGED_NEG = set("DE")
AROMATIC = set("FWYH")
APOLAR = set("AVLIMFWCPG")
POLAR = set("STNQYCHW")


def per_interface(arm: str, level: str = "cluster30",
                  dedup: bool = True) -> list[dict]:
    """Collapse the contact table to one row per physical interface.

    dedup=True keeps only the a_to_b direction for PPI. Without it every PPI
    interface is counted twice, and for the many metrics that are symmetric
    (interface size, bond fractions, closest approach) the two copies are
    identical -- pseudo-replication that halves p-values without adding evidence.
    The audit showed this flipped one metric ('contacts per probe residue') from
    p=5e-4 to p=0.50, so it is corrected here rather than annotated.
    """
    rows = C.load_contacts(arm, level)
    if dedup:
        rows = C.dedup_rows(arm, rows)
    agg: dict = defaultdict(lambda: {
        "n_pairs": 0, "salt": 0, "hbond": 0, "vdw": 0, "stack": 0,
        "n_atom": 0, "probe_res": set(), "target_res": set(),
        "side_hb": 0, "main_hb": 0, "side_vdw": 0, "main_vdw": 0,
        "side_contacts": 0, "main_contacts": 0,
        "opp_charge": 0, "like_charge": 0, "apolar_apolar": 0,
        "aromatic_any": 0, "cation_pi": 0, "dists": [],
    })
    for r in rows:
        key = (r["pdb_id"], r.get("direction", "")) if arm == "ppi" else (r["pdb_id"], "")
        a = agg[key]
        a["unit"] = r["unit"]
        a["weight"] = r["weight"]
        a["n_pairs"] += 1
        for b in C.BOND_COLS:
            a[b] += int(bool(r[b]))
        a["n_atom"] += int(r["n_atom"] or 0)
        a["probe_res"].add((r["ab_chain"], r["ab_resi"]))
        a["target_res"].add((r["ag_chain"], r["ag_resi"]))
        if not np.isnan(r["min_dist"]):
            a["dists"].append(r["min_dist"])

        # moiety mediation: which part of the TARGET residue is being read
        for col, k in (("hbond_ag_side", "side_hb"), ("hbond_ag_main", "main_hb"),
                       ("vdw_ag_side", "side_vdw"), ("vdw_ag_main", "main_vdw"),
                       ("ag_side_contacts", "side_contacts"),
                       ("ag_main_contacts", "main_contacts")):
            v = r.get(col)
            a[k] += int(v) if v not in (None, "") else 0

        p, t = r["aa_probe"], r["aa_target"]
        pos = (p in CHARGED_POS and t in CHARGED_NEG) or (p in CHARGED_NEG and t in CHARGED_POS)
        like = ((p in CHARGED_POS and t in CHARGED_POS)
                or (p in CHARGED_NEG and t in CHARGED_NEG))
        a["opp_charge"] += int(pos)
        a["like_charge"] += int(like)
        a["apolar_apolar"] += int(p in APOLAR and t in APOLAR)
        a["aromatic_any"] += int(p in AROMATIC or t in AROMATIC)
        a["cation_pi"] += int((p in CHARGED_POS and t in AROMATIC)
                              or (t in CHARGED_POS and p in AROMATIC))

    out = []
    for (pdb, direction), a in agg.items():
        n = max(1, a["n_pairs"])
        side_tot = a["side_contacts"] + a["main_contacts"]
        out.append({
            "pdb_id": pdb, "direction": direction, "arm": arm,
            "unit": a["unit"], "weight": a["weight"],
            "n_pairs": a["n_pairs"],
            "n_probe_residues": len(a["probe_res"]),
            "n_target_residues": len(a["target_res"]),
            "n_atom_contacts": a["n_atom"],
            "contacts_per_probe_residue": a["n_pairs"] / max(1, len(a["probe_res"])),
            "mean_min_dist": float(np.mean(a["dists"])) if a["dists"] else np.nan,
            "f_salt": a["salt"] / n,
            "f_hbond": a["hbond"] / n,
            "f_vdw": a["vdw"] / n,
            "f_stack": a["stack"] / n,
            "f_opp_charge": a["opp_charge"] / n,
            "f_like_charge": a["like_charge"] / n,
            "f_apolar_apolar": a["apolar_apolar"] / n,
            "f_aromatic": a["aromatic_any"] / n,
            "f_cation_pi": a["cation_pi"] / n,
            # what fraction of the target is read through its SIDE CHAIN (identity)
            # rather than its BACKBONE (shape) -- the sequence-vs-shape readout axis
            "f_target_sidechain_read": a["side_contacts"] / max(1, side_tot),
            "f_hbond_to_sidechain": a["side_hb"] / max(1, a["side_hb"] + a["main_hb"]),
        })
    return out


def wmean(vals, wts):
    v = np.asarray(vals, float); w = np.asarray(wts, float)
    ok = np.isfinite(v)
    return float(np.sum(v[ok] * w[ok]) / np.sum(w[ok])) if ok.any() else float("nan")


def unit_level_values(recs, field):
    """Average within unit first, then across units -- one vote per unit."""
    by_unit = defaultdict(list)
    for r in recs:
        if np.isfinite(r[field]):
            by_unit[r["unit"]].append(r[field])
    return np.array([np.mean(v) for v in by_unit.values() if v])


def compare(a_vals, b_vals):
    """Mann-Whitney U on unit-level values + rank-biserial effect size."""
    from scipy.stats import mannwhitneyu
    if len(a_vals) < 3 or len(b_vals) < 3:
        return {"p": None, "effect": None}
    u, p = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
    r = 2 * u / (len(a_vals) * len(b_vals)) - 1      # -1..1, 0 = no separation
    return {"p": float(p), "effect_rank_biserial": float(r),
            "median_antibody": float(np.median(a_vals)),
            "median_ppi": float(np.median(b_vals)),
            "n_units_antibody": len(a_vals), "n_units_ppi": len(b_vals)}


FORCE_FIELDS = [
    ("f_hbond", "H-bond"), ("f_salt", "Salt bridge"), ("f_vdw", "van der Waals"),
    ("f_stack", "Aromatic stacking"), ("f_opp_charge", "Opposite-charge pair"),
    ("f_like_charge", "Like-charge pair"), ("f_apolar_apolar", "Apolar-apolar"),
    ("f_aromatic", "Aromatic involved"), ("f_cation_pi", "Cation-pi"),
    ("f_target_sidechain_read", "Read via target side chain"),
]
SIZE_FIELDS = [
    ("n_pairs", "Residue pairs per interface"),
    ("n_probe_residues", "Probe-side interface residues"),
    ("n_target_residues", "Target-side interface residues"),
    ("contacts_per_probe_residue", "Contacts per probe residue"),
    ("n_atom_contacts", "Atom contacts per interface"),
    ("mean_min_dist", "Mean closest-approach distance (A)"),
]


def figure(recs, fields, fname, suptitle, logx=()):
    n = len(fields)
    ncol = 5 if n > 6 else 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.9 * nrow))
    axes = np.atleast_1d(axes).ravel()
    stats = {}
    for ax, (f, label) in zip(axes, fields):
        vals = {arm: unit_level_values([r for r in recs if r["arm"] == arm], f)
                for arm in C.ARMS}
        bins = 30
        if f in logx:
            allv = np.concatenate(list(vals.values()))
            allv = allv[allv > 0]
            bins = np.logspace(np.log10(allv.min()), np.log10(allv.max()), 30)
            ax.set_xscale("log")
        for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
            ax.hist(vals[arm], bins=bins, density=True, alpha=0.55, color=col,
                    label=C.ARM_LABEL[arm] if f == fields[0][0] else None)
            ax.axvline(np.median(vals[arm]), color=col, ls="--", lw=1.2)
        ax.set_title(label, fontsize=9)
        ax.tick_params(labelsize=7)
        stats[f] = compare(vals["antibody"], vals["ppi"])
    for ax in axes[len(fields):]:
        ax.axis("off")
    fig.legend(loc="lower right", fontsize=9, frameon=False)
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(C.OUT_FIG / fname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return stats


if __name__ == "__main__":
    recs = per_interface("antibody") + per_interface("ppi")

    path = C.OUT_INT / "interface_inventory.csv"
    keys = list(recs[0])
    with path.open("w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in recs:
            fh.write(",".join(str(r[k]) for k in keys) + "\n")

    size_stats = figure(recs, SIZE_FIELDS, "phase2_interface_inventory.png",
                        "Interface size and packing (unit-level distributions)",
                        logx=("n_pairs", "n_atom_contacts", "n_probe_residues",
                              "n_target_residues"))
    force_stats = figure(recs, FORCE_FIELDS, "phase2_physicochemistry.png",
                         "The forces at work: bond and chemistry composition "
                         "per interface")

    summary = {
        "n_interfaces": {arm: sum(1 for r in recs if r["arm"] == arm) for arm in C.ARMS},
        "size": size_stats,
        "physicochemistry": force_stats,
        "note": ("all tests are Mann-Whitney U on UNIT-level means (one vote per "
                 "cluster30 unit), not per-interface values, so hyper-solved "
                 "antigens cannot inflate significance. effect_rank_biserial is "
                 "the interpretable quantity: 0 = distributions overlap "
                 "completely, +-1 = complete separation."),
        "salt_bridge_caveat": ("salt-bridge counts are not comparable to the "
                               "nucleic-acid arm: carboxylate vs phosphate acceptor"),
    }
    (C.OUT_MET / "inventory.json").write_text(json.dumps(summary, indent=1))

    print(f"interfaces: {summary['n_interfaces']}")
    print(f"\n{'field':32s} {'Ab median':>10s} {'PPI median':>11s} {'effect':>8s} {'p':>10s}")
    for f, label in SIZE_FIELDS + FORCE_FIELDS:
        s = (size_stats | force_stats)[f]
        if s["p"] is None:
            continue
        print(f"{label:32s} {s['median_antibody']:10.3f} {s['median_ppi']:11.3f} "
              f"{s['effect_rank_biserial']:+8.3f} {s['p']:10.2e}")
