"""Addition 6 - do interfaces overlap functional sites?

The lysozyme case study found antibodies converging on the RIM of the active-site
cleft (Asp101, 69% of antibodies) while largely avoiding the catalytic pair (Glu35
12%, Asp52 19%). This asks the same question corpus-wide, for both arms.

Functional site definition. UniProt residue-level features would need SIFTS
numbering reconciliation for every chain, which is a large surface for silent error.
Instead the functional site is defined structurally and self-containedly: residues
within 4.5 A of a bound BIOLOGICAL LIGAND -- a heteroatom group that is not water,
not a buffer/cryoprotectant, and not a lone ion. This is a proxy, and is labelled as
one, but it is computed identically for both arms so the comparison is fair.

The question it answers: are antibody epitopes drawn to functional sites (relevant
to neutralisation) or merely to accessible surface? And do natural partners -- which
often ARE the functional interaction -- overlap them more?

Outputs
  metadata/analysis/functional_sites.json
  figures/add6_functional_sites.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import random
import sys
import time
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import NeighborSearch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import phase6a_structural as P

CUTOFF = 4.5
PPI_SAMPLE = 1000
SEED = 20260816

# crystallisation additives, cryoprotectants and lone ions: present because of how
# the crystal was grown, not because they are the protein's biological ligand
IGNORE = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "PGE", "P6G", "1PE", "MPD", "DMS",
    "TRS", "MES", "EPE", "ACT", "ACY", "FMT", "CIT", "TLA", "SO4", "PO4", "NO3",
    "CL", "BR", "IOD", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "NI", "CU", "CD",
    "HG", "CO", "CS", "RB", "SR", "BA", "LI", "F", "OXY", "UNX", "UNL", "DTT",
    "BME", "IMD", "AZI", "SCN", "NH4", "PGO", "BU3", "MRD", "XYP",
}
MIN_LIGAND_ATOMS = 6          # exclude tiny fragments


def scan(arm, pdb_id, probe_chains, target_chains, iface):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return None
    model = P.load_model(path)

    lig_atoms, lig_names = [], []
    prot_atoms = []
    for ch in model:
        for res in ch:
            name = res.get_resname().strip()
            het = res.get_id()[0]
            if name in C.AA3_TO_1:
                if ch.id in target_chains:
                    for at in res:
                        if at.element != "H":
                            prot_atoms.append((at, ch.id, res))
                continue
            if het == " " or name in IGNORE:
                continue
            atoms = [a for a in res if a.element != "H"]
            if len(atoms) < MIN_LIGAND_ATOMS:
                continue
            lig_atoms += atoms
            lig_names.append(name)

    if not lig_atoms or not prot_atoms:
        return None

    ns = NeighborSearch([a[0] for a in prot_atoms])
    owner = {id(a[0]): (a[1], P.residue_key(a[2])) for a in prot_atoms}
    site = set()
    for at in lig_atoms:
        for hit in ns.search(at.coord, CUTOFF):
            site.add(owner[id(hit)])
    if len(site) < 3:
        return None

    epitope = {(ch, resi) for (ch, resi) in iface if ch in target_chains}
    if len(epitope) < 5:
        return None
    inter = site & epitope
    return {
        "arm": arm, "pdb_id": pdb_id,
        "n_site": len(site), "n_interface": len(epitope), "n_shared": len(inter),
        "jaccard": round(len(inter) / len(site | epitope), 4),
        "frac_interface_in_site": round(len(inter) / len(epitope), 4),
        "frac_site_in_interface": round(len(inter) / len(site), 4),
        "ligands": "|".join(sorted(set(lig_names))[:5]),
    }


if __name__ == "__main__":
    rng = random.Random(SEED)
    plan = []
    for arm in C.ARMS:
        tg = P.targets(arm)
        ids = sorted(tg)
        if arm == "ppi" and len(ids) > PPI_SAMPLE:
            ids = sorted(rng.sample(ids, PPI_SAMPLE))
        for p in ids:
            pc, tc, iface = tg[p]
            plan.append((arm, p, pc, tc, iface))
    print(f"scanning {len(plan)} structures for bound biological ligands", flush=True)

    rows, t0 = [], time.time()
    for k, (arm, pdb, pc, tc, iface) in enumerate(plan, 1):
        try:
            r = scan(arm, pdb, pc, tc, iface)
            if r:
                rows.append(r)
        except Exception:                                        # noqa: BLE001
            pass
        if k % 300 == 0:
            print(f"  {k}/{len(plan)}  {time.time()-t0:4.0f}s", flush=True)

    units = {}
    for arm in C.ARMS:
        for p, v in C.units(arm).items():
            units[(arm, p)] = (v.get("cluster30") if arm == "antibody"
                               else (v.get("a_to_b") or {}).get("cluster30"))
    per = {arm: defaultdict(list) for arm in C.ARMS}
    for r in rows:
        u = units.get((r["arm"], r["pdb_id"]))
        if u:
            per[r["arm"]][u].append(r["frac_interface_in_site"])
    vals = {arm: np.array([np.mean(v) for v in per[arm].values() if v])
            for arm in C.ARMS}

    from scipy.stats import mannwhitneyu
    stat = {}
    if len(vals["antibody"]) > 5 and len(vals["ppi"]) > 5:
        u_, p_ = mannwhitneyu(vals["antibody"], vals["ppi"], alternative="two-sided")
        stat = {"median_antibody": float(np.median(vals["antibody"])),
                "median_ppi": float(np.median(vals["ppi"])),
                "effect_rank_biserial": float(
                    2 * u_ / (len(vals["antibody"]) * len(vals["ppi"])) - 1),
                "p": float(p_),
                "n_units_antibody": int(len(vals["antibody"])),
                "n_units_ppi": int(len(vals["ppi"]))}

    out = {
        "n_structures_with_ligand": len(rows),
        "n_scanned": len(plan),
        "by_arm": {arm: int(sum(1 for r in rows if r["arm"] == arm))
                   for arm in C.ARMS},
        "frac_interface_in_functional_site": stat,
        "commonest_ligands": Counter(
            l for r in rows for l in r["ligands"].split("|") if l).most_common(12),
        "definition": (f"functional site = target-chain residues within {CUTOFF} A "
                       f"of a bound heteroatom group with >={MIN_LIGAND_ATOMS} heavy "
                       f"atoms, excluding water, buffers, cryoprotectants and lone "
                       f"ions. A structural PROXY for a functional site, computed "
                       f"identically for both arms."),
        "caveat": ("Only structures that happen to be solved WITH a ligand "
                   "contribute, which is not a random subset -- enzymes and "
                   "receptors are over-represented."),
    }
    (C.OUT_MET / "functional_sites.json").write_text(json.dumps(out, indent=1))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    for arm, col in (("antibody", "#b2182b"), ("ppi", "#2166ac")):
        v = vals[arm]
        if len(v):
            ax.hist(v, bins=np.linspace(0, 1, 26), density=True, alpha=0.55,
                    color=col, label=f"{C.ARM_LABEL[arm]} (median "
                                     f"{np.median(v):.3f}, n={len(v)} units)")
            ax.axvline(np.median(v), color=col, ls="--", lw=1.3)
    ax.set_xlabel("fraction of interface residues inside a ligand-binding site")
    ax.set_ylabel("density")
    if stat:
        ax.set_title(f"effect {stat['effect_rank_biserial']:+.2f}, "
                     f"p={stat['p']:.2e}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    labs, cnts = zip(*out["commonest_ligands"][:10]) if out["commonest_ligands"] \
        else ((), ())
    if labs:
        ax.barh(range(len(labs)), cnts, color="#1b7837")
        ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=8)
        ax.invert_yaxis()
    ax.set_xlabel("structures")
    ax.set_title("Commonest bound ligands defining the sites", fontsize=10)

    fig.suptitle("Addition 6 - do interfaces overlap ligand-binding (functional) "
                 "sites?", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(C.OUT_FIG / "add6_functional_sites.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    print(f"\n{len(rows)} structures had a usable ligand "
          f"({out['by_arm']})")
    if stat:
        print(f"fraction of interface inside a functional site: antibody "
              f"{stat['median_antibody']:.3f} vs PPI {stat['median_ppi']:.3f}  "
              f"effect {stat['effect_rank_biserial']:+.3f}  p={stat['p']:.2e}")
    print("commonest ligands:", out["commonest_ligands"][:8])
