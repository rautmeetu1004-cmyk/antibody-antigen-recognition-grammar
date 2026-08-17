"""Addition 1 - the predictive payoff: does a general-PPI potential work on antibodies?

Everything so far is descriptive: two grammars differ (permutation p = 0.001) and the
difference survives secondary-structure, fold, resolution and redundancy controls.
This asks the question that makes the difference matter:

  IF YOU SCORE AN ANTIBODY INTERFACE WITH A POTENTIAL LEARNED FROM GENERAL
  PROTEIN-PROTEIN INTERFACES, DOES IT STILL WORK?

If it does, the grammar difference is a curiosity. If it does not, then antibody
design tools built on general interface statistics are mis-specified, and that is a
concrete, testable consequence.

Method. Each 20x20 log-odds map is used as a threading potential: an interface's
score is the mean log-odds over its residue-residue contacts. Decoys keep the
interface GEOMETRY exactly and permute the residue IDENTITIES on the target side
among its own interface positions -- so a decoy has the same number of contacts, the
same contact graph, and the same residue composition, differing only in which
residue sits at which position. That isolates recognition chemistry from interface
size, which a rigid-body docking decoy would confound (and no docking software is
available here anyway).

Discrimination is measured as AUC over decoys per interface, then averaged over
cluster30 UNITS so hyper-solved antigens cannot dominate.

Two potentials are applied to each arm:
  SELF   built leave-one-unit-out from the same arm (no interface scores itself)
  CROSS  built from the other arm entirely

Outputs
  figures/add1_transfer.png
  metadata/analysis/transfer.json
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

N_DECOY = 50
SEED = C.SEED
IDX = {a: i for i, a in enumerate(C.AAS)}


def potential(rows, symmetrise):
    return C.log2_odds(C.contact_matrix(rows, symmetrise=symmetrise))


def interfaces(arm):
    """{(pdb, direction) -> {unit, contacts:[(probe_aa, target_res_key)]}}"""
    out = defaultdict(lambda: {"unit": None, "contacts": [], "target_res": {}})
    for r in C.dedup_rows(arm, C.load_contacts(arm)):
        key = (r["pdb_id"], r.get("direction", ""))
        e = out[key]
        e["unit"] = r["unit"]
        tkey = (r["ag_chain"], r["ag_resi"])
        e["target_res"][tkey] = r["aa_target"]
        e["contacts"].append((r["aa_probe"], tkey))
    return {k: v for k, v in out.items() if len(v["contacts"]) >= 10}


def score(contacts, target_ids, pot):
    """Mean log-odds over contacts, given an assignment of identities to positions."""
    tot, n = 0.0, 0
    for probe_aa, tkey in contacts:
        t = target_ids.get(tkey)
        if t is None:
            continue
        v = pot[IDX[probe_aa], IDX[t]]
        if np.isfinite(v):
            tot += v
            n += 1
    return tot / n if n else np.nan


def auc_for_interface(entry, pot, rng, n_decoy=N_DECOY):
    contacts = entry["contacts"]
    native = entry["target_res"]
    s_nat = score(contacts, native, pot)
    if not np.isfinite(s_nat):
        return np.nan
    keys = list(native)
    vals = [native[k] for k in keys]
    wins = ties = 0
    for _ in range(n_decoy):
        perm = list(vals)
        rng.shuffle(perm)
        s_dec = score(contacts, dict(zip(keys, perm)), pot)
        if not np.isfinite(s_dec):
            continue
        if s_nat > s_dec:
            wins += 1
        elif s_nat == s_dec:
            ties += 1
    return (wins + 0.5 * ties) / n_decoy


def evaluate(arm, pot_self_rows, pot_cross, rng):
    """AUC per unit, for a leave-one-unit-out self potential and a cross potential.

    The LOO potential is the TOTAL unit-weighted matrix minus the held-out unit's
    own matrix, which is exactly equivalent to rebuilding from the other units but
    costs one subtraction instead of a full pass over ~400k rows per unit.
    """
    ifaces = interfaces(arm)
    sym = (arm == "ppi")
    by_unit = defaultdict(list)
    for k, e in ifaces.items():
        by_unit[e["unit"]].append(k)

    ids, U = C.unit_matrices(C.dedup_rows(arm, C.load_contacts(arm)), symmetrise=sym)
    pos = {u: i for i, u in enumerate(ids)}
    total = U.sum(0)

    self_auc, cross_auc = defaultdict(list), defaultdict(list)
    for unit, keys in by_unit.items():
        i = pos.get(unit)
        pot_loo = C.log2_odds(total - U[i]) if i is not None else C.log2_odds(total)
        for k in keys:
            e = ifaces[k]
            a = auc_for_interface(e, pot_loo, rng)
            b = auc_for_interface(e, pot_cross, rng)
            if np.isfinite(a):
                self_auc[unit].append(a)
            if np.isfinite(b):
                cross_auc[unit].append(b)
    su = np.array([np.mean(v) for v in self_auc.values() if v])
    cu = np.array([np.mean(v) for v in cross_auc.values() if v])
    return su, cu, len(ifaces)


def summarise(vals, label):
    from scipy.stats import wilcoxon
    m = float(np.mean(vals))
    ci = [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]
    try:
        p = float(wilcoxon(vals - 0.5).pvalue)
    except ValueError:
        p = float("nan")
    return {"label": label, "mean_auc": m, "unit_ci": ci,
            "n_units": int(len(vals)), "p_vs_chance": p,
            "frac_units_below_chance": float(np.mean(vals < 0.5))}


if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    ab_rows = C.dedup_rows("antibody", C.load_contacts("antibody"))
    pp_rows = C.dedup_rows("ppi", C.load_contacts("ppi"))
    pot_ab = potential(ab_rows, False)
    pot_pp = potential(pp_rows, True)

    print("scoring antibody interfaces...", flush=True)
    ab_self, ab_cross, n_ab = evaluate("antibody", ab_rows, pot_pp, rng)
    print("scoring PPI interfaces...", flush=True)
    pp_self, pp_cross, n_pp = evaluate("ppi", pp_rows, pot_ab, rng)

    from scipy.stats import mannwhitneyu
    out = {
        "n_interfaces": {"antibody": n_ab, "ppi": n_pp},
        "n_decoys_per_interface": N_DECOY,
        "antibody_scored_by_antibody_potential": summarise(ab_self, "self (LOO)"),
        "antibody_scored_by_ppi_potential": summarise(ab_cross, "cross (PPI)"),
        "ppi_scored_by_ppi_potential": summarise(pp_self, "self (LOO)"),
        "ppi_scored_by_antibody_potential": summarise(pp_cross, "cross (antibody)"),
        "transfer_loss_antibody": float(np.mean(ab_self) - np.mean(ab_cross)),
        "transfer_loss_ppi": float(np.mean(pp_self) - np.mean(pp_cross)),
        "antibody_self_vs_cross_p": float(
            mannwhitneyu(ab_self, ab_cross, alternative="two-sided").pvalue),
        "ppi_self_vs_cross_p": float(
            mannwhitneyu(pp_self, pp_cross, alternative="two-sided").pvalue),
        "method_note": (
            "Decoys preserve interface geometry, contact count and residue "
            "composition, permuting only WHICH residue occupies WHICH position. AUC "
            "is averaged within cluster30 unit before any statistic. The self "
            "potential is leave-one-unit-out, so no interface contributes to the "
            "potential that scores it."),
    }
    (C.OUT_MET / "transfer.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (arm, s, c) in zip(axes, (("antibody", ab_self, ab_cross),
                                      ("ppi", pp_self, pp_cross))):
        bins = np.linspace(0, 1, 41)
        ax.hist(s, bins=bins, alpha=0.6, color="#1b7837",
                label=f"self potential (LOO), mean {np.mean(s):.3f}")
        ax.hist(c, bins=bins, alpha=0.6, color="#b2182b",
                label=f"cross potential, mean {np.mean(c):.3f}")
        ax.axvline(0.5, ls="--", c="0.35", lw=1.2, label="chance")
        ax.set_xlabel("native-vs-decoy AUC (per cluster30 unit)")
        ax.set_ylabel("units")
        ax.set_title(f"{C.ARM_LABEL[arm]} interfaces", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Addition 1 - does a general-PPI potential transfer to antibody "
                 "interfaces?", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(C.OUT_FIG / "add1_transfer.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for k in ("antibody_scored_by_antibody_potential",
              "antibody_scored_by_ppi_potential",
              "ppi_scored_by_ppi_potential",
              "ppi_scored_by_antibody_potential"):
        d = out[k]
        print(f"{k:42s} AUC={d['mean_auc']:.3f} "
              f"[{d['unit_ci'][0]:.3f},{d['unit_ci'][1]:.3f}] "
              f"units={d['n_units']:4d} p={d['p_vs_chance']:.2e} "
              f"below_chance={d['frac_units_below_chance']:.2f}")
    print(f"\ntransfer loss antibody = {out['transfer_loss_antibody']:+.3f} "
          f"(p={out['antibody_self_vs_cross_p']:.2e})")
    print(f"transfer loss ppi      = {out['transfer_loss_ppi']:+.3f} "
          f"(p={out['ppi_self_vs_cross_p']:.2e})")
