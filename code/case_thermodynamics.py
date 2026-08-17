"""Case study - thermodynamic signature of a convergent immune response.

The corpus-wide result is that antibody association buries proportionally less
apolar surface than co-evolved association, giving a less negative heat-capacity
change. That is a statement about two distributions. It leaves open a question that
only a single-antigen series can answer:

  When many independent antibodies converge on one face of one antigen, do they also
  converge on one thermodynamic strategy, or do chemically different solutions to the
  same surface coexist?

The 42 lysozyme complexes provide the series. For each we already have buried apolar
and polar surface at atom resolution, and the set of antigen positions contacted. The
analysis asks three things: how tightly the thermodynamic signature clusters within
one antigen relative to the corpus at large, whether antibodies engaging the
immunodominant face differ thermodynamically from those engaging elsewhere, and
whether buried surface tracks the number of dominant positions contacted.

Outputs
  metadata/analysis/case_thermodynamics.json
  figures/case_thermodynamics.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import immunodominance as IM

DOMINANT_TOP_N = 12          # positions defining the immunodominant face


def thermo_by_pdb():
    out = {}
    path = C.OUT_INT / "thermodynamics.csv"
    for r in csv.DictReader(path.open()):
        if r["arm"] != "antibody":
            continue
        ap, po = float(r["dA_apolar"]), float(r["dA_polar"])
        tot = ap + po
        if tot <= 0:
            continue
        dcp = 0.45 * (-ap) + (-0.26) * (-po)
        out[r["pdb_id"]] = {"dA_apolar": ap, "dA_polar": po, "dA_total": tot,
                            "frac_apolar": ap / tot,
                            "dCp": dcp, "dCp_per_1000": dcp / (tot / 1000.0)}
    return out


if __name__ == "__main__":
    units = C.units("antibody")
    rows = C.dedup_rows("antibody", C.load_contacts("antibody"))
    pdbs = {r["pdb_id"] for r in rows}
    from collections import Counter
    cluster = Counter(units[p]["cluster30"] for p in pdbs
                      if p in units).most_common(1)[0][0]

    # epitope positions in canonical numbering, per structure
    by_pdb = defaultdict(lambda: defaultdict(set))
    for r in rows:
        if units.get(r["pdb_id"], {}).get("cluster30") == cluster:
            by_pdb[r["pdb_id"]][r["ag_chain"]].add(r["ag_resi"].strip())

    pos_by_pdb, hits = {}, np.zeros(len(IM.HEWL))
    for pdb, chains in sorted(by_pdb.items()):
        got = set()
        for ch, resis in chains.items():
            rm = IM.chain_residues(pdb, ch)
            off, score = IM.best_offset(rm)
            if off is None or score < 0.8:
                continue
            for x in resis:
                v = rm.get(x)
                if v:
                    pp = v[0] - off + 1
                    if 1 <= pp <= len(IM.HEWL):
                        got.add(pp)
        if got:
            pos_by_pdb[pdb] = got
            for pp in got:
                hits[pp - 1] += 1

    dominant = set(int(i + 1) for i in np.argsort(-hits)[:DOMINANT_TOP_N])
    th = thermo_by_pdb()

    recs = []
    for pdb, pos in pos_by_pdb.items():
        t = th.get(pdb)
        if not t:
            continue
        recs.append({"pdb_id": pdb, **t,
                     "n_dominant_contacted": len(pos & dominant),
                     "frac_dominant": len(pos & dominant) / len(dominant),
                     "n_epitope_positions": len(pos)})
    if not recs:
        raise SystemExit("no lysozyme structures with thermodynamic data")

    dcp = np.array([r["dCp_per_1000"] for r in recs])
    fap = np.array([r["frac_apolar"] for r in recs])
    ndom = np.array([r["n_dominant_contacted"] for r in recs])
    tot = np.array([r["dA_total"] for r in recs])

    # corpus-wide antibody distribution, for a dispersion comparison
    corpus = np.array([v["dCp_per_1000"] for k, v in th.items()])
    cv_case = float(np.std(dcp) / abs(np.mean(dcp)))
    cv_corpus = float(np.std(corpus) / abs(np.mean(corpus)))

    # engaging the dominant face vs not
    hi = dcp[ndom >= np.median(ndom)]
    lo = dcp[ndom < np.median(ndom)]
    if len(hi) > 3 and len(lo) > 3:
        u, p_face = mannwhitneyu(hi, lo, alternative="two-sided")
        eff_face = float(2 * u / (len(hi) * len(lo)) - 1)
    else:
        p_face, eff_face = float("nan"), float("nan")

    rho_dom, p_dom = spearmanr(ndom, dcp)
    rho_size, p_size = spearmanr(tot, dcp)

    out = {
        "antigen": "Lysozyme C",
        "n_structures": len(recs),
        "dominant_positions": sorted(dominant),
        "dCp_per_1000A2": {
            "median": float(np.median(dcp)), "iqr": [float(np.percentile(dcp, 25)),
                                                     float(np.percentile(dcp, 75))],
            "range": [float(dcp.min()), float(dcp.max())],
            "coefficient_of_variation": cv_case},
        "corpus_comparison": {
            "corpus_median": float(np.median(corpus)),
            "corpus_coefficient_of_variation": cv_corpus,
            "note": ("A case-series CV substantially below the corpus CV would mean "
                     "that antibodies against one antigen converge thermodynamically "
                     "as well as spatially. A comparable CV would mean the same face "
                     "is engaged by thermodynamically diverse solutions.")},
        "frac_apolar_buried": {"median": float(np.median(fap)),
                               "range": [float(fap.min()), float(fap.max())]},
        "dominant_face_engagement": {
            "median_n_dominant_contacted": float(np.median(ndom)),
            "effect_high_vs_low": eff_face, "p": float(p_face),
            "spearman_rho_dCp_vs_n_dominant": float(rho_dom), "p_spearman": float(p_dom)},
        "dCp_vs_interface_size": {"spearman_rho": float(rho_size),
                                  "p": float(p_size)},
        "per_structure": sorted(recs, key=lambda r: r["dCp_per_1000"]),
    }
    (C.OUT_MET / "case_thermodynamics.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    ax = axes[0]
    bins = np.linspace(min(corpus.min(), dcp.min()), max(corpus.max(), dcp.max()), 40)
    ax.hist(corpus, bins=bins, density=True, alpha=0.45, color="#2166ac",
            label=f"all antibody interfaces (n={len(corpus)})")
    ax.hist(dcp, bins=bins, density=True, alpha=0.7, color="#b2182b",
            label=f"lysozyme series (n={len(dcp)})")
    ax.axvline(np.median(corpus), color="#2166ac", ls="--", lw=1.2)
    ax.axvline(np.median(dcp), color="#b2182b", ls="--", lw=1.2)
    ax.set_xlabel("$\\Delta C_p$ per 1000 $\\AA^2$ buried (cal mol$^{-1}$ K$^{-1}$)")
    ax.set_ylabel("density")
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    ax.scatter(ndom, dcp, s=34, color="#762a83", alpha=0.8, edgecolors="none")
    if len(ndom) > 3:
        z = np.polyfit(ndom, dcp, 1)
        xs = np.linspace(ndom.min(), ndom.max(), 20)
        ax.plot(xs, np.polyval(z, xs), c="0.35", lw=1.3, ls="--")
    ax.set_xlabel(f"immunodominant positions contacted (of {DOMINANT_TOP_N})")
    ax.set_ylabel("$\\Delta C_p$ per 1000 $\\AA^2$")

    ax = axes[2]
    ax.scatter(fap, tot, s=34, color="#1b7837", alpha=0.8, edgecolors="none")
    ax.set_xlabel("apolar fraction of buried surface")
    ax.set_ylabel("total buried surface ($\\AA^2$)")

    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "case_thermodynamics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"lysozyme series: n={len(recs)}")
    print(f"dCp per 1000 A2: median {np.median(dcp):.1f}, IQR "
          f"[{np.percentile(dcp,25):.1f}, {np.percentile(dcp,75):.1f}], range "
          f"[{dcp.min():.1f}, {dcp.max():.1f}]")
    print(f"CV case series {cv_case:.3f} vs corpus {cv_corpus:.3f}")
    print(f"apolar fraction: median {np.median(fap):.3f}, range "
          f"[{fap.min():.3f}, {fap.max():.3f}]")
    print(f"dCp vs n dominant contacted: rho={rho_dom:+.3f} p={p_dom:.3f}")
    print(f"dominant-face high vs low: effect {eff_face:+.3f} p={p_face:.3f}")
    print(f"dCp vs interface size: rho={rho_size:+.3f} p={p_size:.3g}")
