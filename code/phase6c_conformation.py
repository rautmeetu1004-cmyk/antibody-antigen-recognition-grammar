"""Phase 6c - secondary structure, dihedrals and rotamers at the interface.

Consumes the Phase 6a structural pass and answers three things:

  1. THE LOOP CONFOUND (the reason this phase exists). Antibody paratopes are CDR
     loops. If the "antibody grammar" of Phases 3-6 is really LOOP grammar, then
     restricting BOTH arms to loop-loop contacts should collapse the difference.
     This is the same class of test as the nucleic-acid study's target-composition
     control, which reframed that paper's headline. It is run before any
     interpretation of the antibody signature is allowed to stand.

  2. TORSION ANGLES (your question). phi/psi of the residues being read -- is the
     antibody reading backbone in conformations that general interfaces do not? --
     and chi1 rotamer distributions.

  3. ENTROPIC COST. Rotamer restriction at the interface versus away from it. If
     antibodies bind enthalpically (Phase 2's hypothesis), their interface residues
     should not be paying an unusually large side-chain entropy penalty.

Outputs
  figures/phase6c_conformation.png
  figures/phase6c_loop_confound.png
  metadata/analysis/conformation.json
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
from scipy.stats import chi2_contingency, mannwhitneyu

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

SS3 = ("helix", "strand", "loop")


def load_structural():
    """{(arm, pdb, chain, resi) -> row}, plus the full row list."""
    path = C.OUT_INT / "structural_residues.csv"
    rows = list(csv.DictReader(path.open()))
    for r in rows:
        for k in ("phi", "psi", "chi1", "chi2", "sasa_complex", "sasa_free", "d_sasa"):
            r[k] = float(r[k]) if r[k] not in ("", None) else np.nan
        r["is_interface"] = int(r["is_interface"])
    idx = {(r["arm"], r["pdb_id"], r["chain"], r["resi"]): r for r in rows}
    return rows, idx


def annotate_contacts(arm, idx, level="cluster30"):
    """Attach probe-side and target-side secondary structure to each contact."""
    out = []
    for r in C.dedup_rows(arm, C.load_contacts(arm, level)):
        a = idx.get((arm, r["pdb_id"], r["ab_chain"], r["ab_resi"].strip()))
        b = idx.get((arm, r["pdb_id"], r["ag_chain"], r["ag_resi"].strip()))
        if a is None or b is None:
            continue                     # structure not in the Phase 6a sample
        r = dict(r)
        r["ss_probe"], r["ss_target"] = a["ss3"], b["ss3"]
        out.append(r)
    return out


def ss_composition(rows_by_arm):
    """Interface secondary-structure composition, per arm and per side."""
    out = {}
    for arm, rows in rows_by_arm.items():
        for side, key in (("probe", "ss_probe"), ("target", "ss_target")):
            # weight by unit so hyper-solved antigens cannot dominate
            w = defaultdict(float)
            for r in rows:
                w[r[key]] += r["weight"]
            tot = sum(w.values()) or 1.0
            out[f"{arm}_{side}"] = {s: round(w.get(s, 0.0) / tot, 4) for s in SS3}
    return out


def loop_confound(rows_by_arm):
    """Recompute the cross-arm grammar correlation within each SS stratum.

    If the antibody-vs-PPI divergence is really a loop-vs-not-loop artefact, the
    correlation should rise sharply once both arms are restricted to the same
    structural context.
    """
    res = {}
    for stratum, pred in (
        ("all", lambda r: True),
        ("loop_probe", lambda r: r["ss_probe"] == "loop"),
        ("loop_both", lambda r: r["ss_probe"] == "loop" and r["ss_target"] == "loop"),
        ("target_loop", lambda r: r["ss_target"] == "loop"),
        ("target_strand", lambda r: r["ss_target"] == "strand"),
        ("target_helix", lambda r: r["ss_target"] == "helix"),
    ):
        mats, ns = {}, {}
        for arm, rows in rows_by_arm.items():
            sub = [r for r in rows if pred(r)]
            ns[arm] = len(sub)
            mats[arm] = C.log2_odds(
                C.contact_matrix(sub, symmetrise=(arm == "ppi")))
        res[stratum] = {
            "r": C.matrix_similarity(mats["antibody"], mats["ppi"]),
            "n_antibody": ns["antibody"], "n_ppi": ns["ppi"],
        }
    return res


def torsion_stats(rows):
    """phi/psi and chi1 summaries for interface vs non-interface residues."""
    out = {}
    for arm in C.ARMS:
        a = [r for r in rows if r["arm"] == arm]
        for label, sel in (("interface", [r for r in a if r["is_interface"]]),
                           ("background", [r for r in a if not r["is_interface"]])):
            chi = np.array([r["chi1"] for r in sel if np.isfinite(r["chi1"])])
            rot = Counter(r["rotamer"] for r in sel if r["rotamer"])
            tot = sum(rot.values()) or 1
            ds = np.array([r["d_sasa"] for r in sel if np.isfinite(r["d_sasa"])])
            out[f"{arm}_{label}"] = {
                "n": len(sel),
                "rotamer_fractions": {k: round(v / tot, 4)
                                      for k, v in sorted(rot.items())},
                "chi1_circular_spread": float(round(
                    1 - abs(np.mean(np.exp(1j * np.radians(chi)))), 4)) if len(chi) else None,
                "d_sasa_median": float(np.median(ds)) if len(ds) else None,
            }
    return out


def figure_confound(conf, comp):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))

    ax = axes[0]
    labels = ["all contacts", "probe loop", "loop-loop",
              "target loop", "target strand", "target helix"]
    keys = ["all", "loop_probe", "loop_both", "target_loop", "target_strand",
            "target_helix"]
    vals = [conf[k]["r"] for k in keys]
    cols = ["0.35"] + ["#1b7837"] * 2 + ["#2166ac"] * 3
    ax.bar(range(len(keys)), vals, color=cols)
    for i, (k, v) in enumerate(zip(keys, vals)):
        ax.text(i, v + 0.012, f"{v:.2f}", ha="center", fontsize=8)
        ax.text(i, 0.01, f"n={conf[k]['n_antibody']:,}", ha="center", fontsize=6,
                rotation=90, color="w")
    ax.axhline(conf["all"]["r"], ls="--", c="0.5", lw=1)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("antibody vs PPI grammar correlation (r)")
    ax.set_title("The loop confound: does matching structural context\n"
                 "erase the antibody-PPI difference?", fontsize=10)

    ax = axes[1]
    groups = [f"{a}_{s}" for a in C.ARMS for s in ("probe", "target")]
    names = ["Antibody\nparatope", "Antibody\nepitope", "PPI\npartner A",
             "PPI\npartner B"]
    bottom = np.zeros(len(groups))
    for s, col in zip(SS3, ("#b2182b", "#2166ac", "#1b7837")):
        v = np.array([comp[g][s] for g in groups])
        ax.bar(range(len(groups)), v, bottom=bottom, color=col, label=s)
        for i, (b, x) in enumerate(zip(bottom, v)):
            if x > 0.06:
                ax.text(i, b + x / 2, f"{x:.0%}", ha="center", va="center",
                        fontsize=8, color="w")
        bottom += v
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("fraction of interface contacts (unit-weighted)")
    ax.set_title("Secondary structure at the interface", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Phase 6c - is the antibody signature really a loop signature?",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(C.OUT_FIG / "phase6c_loop_confound.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


def figure_torsion(rows):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.4))
    for r, arm in enumerate(C.ARMS):
        a = [x for x in rows if x["arm"] == arm]
        iface = [x for x in a if x["is_interface"]]
        bg = [x for x in a if not x["is_interface"]]

        ax = axes[r, 0]
        phi = np.array([x["phi"] for x in iface if np.isfinite(x["phi"])])
        psi = np.array([x["psi"] for x in iface if np.isfinite(x["psi"])])
        n = min(len(phi), len(psi))
        ax.hexbin(phi[:n], psi[:n], gridsize=45, cmap="viridis", bins="log",
                  extent=(-180, 180, -180, 180))
        ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
        ax.set_xlabel("$\\phi$ (deg)"); ax.set_ylabel("$\\psi$ (deg)")
        ax.set_title(f"{C.ARM_LABEL[arm]}\nRamachandran, interface residues",
                     fontsize=9)

        ax = axes[r, 1]
        for sel, lab, col in ((iface, "interface", "#b2182b"),
                              (bg, "background", "#2166ac")):
            ch = np.array([x["chi1"] for x in sel if np.isfinite(x["chi1"])])
            if len(ch):
                ax.hist(ch, bins=72, range=(-180, 180), density=True, alpha=0.55,
                        color=col, label=lab)
        ax.set_xlabel("$\\chi_1$ (deg)"); ax.set_ylabel("density")
        ax.set_title(f"{C.ARM_LABEL[arm]}\nside-chain $\\chi_1$", fontsize=9)
        ax.legend(fontsize=7, frameon=False)

        ax = axes[r, 2]
        for sel, lab, col in ((iface, "interface", "#b2182b"),
                              (bg, "background", "#2166ac")):
            ds = np.array([x["d_sasa"] for x in sel if np.isfinite(x["d_sasa"])])
            ds = ds[ds >= 0]
            if len(ds):
                ax.hist(ds, bins=50, range=(0, 150), density=True, alpha=0.55,
                        color=col, label=lab)
        ax.set_xlabel("$\\Delta$SASA buried by partner ($\\AA^2$)")
        ax.set_title(f"{C.ARM_LABEL[arm]}\nburial", fontsize=9)
        ax.legend(fontsize=7, frameon=False)

    fig.suptitle("Phase 6c - backbone dihedrals, side-chain rotamers, and burial",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(C.OUT_FIG / "phase6c_conformation.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    rows, idx = load_structural()
    print(f"structural rows: {len(rows):,}  "
          f"structures: {len({(r['arm'], r['pdb_id']) for r in rows}):,}")

    ann = {arm: annotate_contacts(arm, idx) for arm in C.ARMS}
    for arm in C.ARMS:
        print(f"  {arm}: {len(ann[arm]):,} contacts annotated with SS")

    comp = ss_composition(ann)
    conf = loop_confound(ann)
    tors = torsion_stats(rows)

    # is paratope SS composition different from PPI partner SS composition?
    tab = np.array([[sum(1 for r in ann[arm] if r["ss_probe"] == s) for s in SS3]
                    for arm in C.ARMS])
    chi2, p, _, _ = chi2_contingency(tab)
    ss_test = {"chi2": float(chi2), "p": float(p),
               "cramers_v": float(np.sqrt(chi2 / tab.sum()))}

    out = {"ss_composition": comp, "loop_confound": conf, "torsion": tors,
           "ss_composition_test": ss_test,
           "note": ("loop_confound compares the SAME statistic (cross-arm grammar r) "
                    "within matched structural context. A large rise from 'all' to "
                    "'loop_both' would mean the antibody signature is substantially a "
                    "loop-context artefact; a flat profile means it is not.")}
    (C.OUT_MET / "conformation.json").write_text(json.dumps(out, indent=1))

    figure_confound(conf, comp)
    figure_torsion(rows)

    print("\nsecondary structure at the interface (unit-weighted):")
    for k, v in comp.items():
        print(f"  {k:18s} " + "  ".join(f"{s}={v[s]:.3f}" for s in SS3))
    print(f"\nparatope vs PPI-probe SS composition: Cramer's V="
          f"{ss_test['cramers_v']:.3f}, p={ss_test['p']:.2e}")
    print("\nloop confound - cross-arm grammar r within matched context:")
    for k, v in conf.items():
        print(f"  {k:14s} r={v['r']:+.3f}  (n_ab={v['n_antibody']:,}, "
              f"n_ppi={v['n_ppi']:,})")
    print("\nrotamer / burial:")
    for k, v in tors.items():
        print(f"  {k:22s} n={v['n']:6d}  rot={v['rotamer_fractions']}  "
              f"dSASA_med={v['d_sasa_median']}")
