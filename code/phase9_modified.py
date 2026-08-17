"""Phase 9 - modified amino acids and glycans at the interface.

Established during the design phase: the inherited contact tables contain NO
modified residues. Only the 20 standard amino acids survived extraction, plus MSE
(selenomethionine) in the PPI arm -- and MSE is a crystallographic substitution for
Met, not a biological modification. So this question cannot be answered from those
tables at all; it requires going back to the coordinates.

Two categories are recovered here, and they are different things:

  PTMs        post-translationally modified residues that remain part of the
              polymer chain -- phospho (SEP/TPO/PTR), methylated and acetylated
              lysine (MLY/M3L/MLZ/ALY), oxidised cysteine (CSO/CSD/OCS/CME),
              hydroxyproline (HYP), carboxylated residues (KCX/CGU), sulfotyrosine
              (TYS), citrulline (CIR), pyroglutamate (PCA), and anything the
              structure's own MODRES records declare.
  GLYCANS     N-/O-linked sugars (NAG, MAN, BMA, FUC, GAL, SIA, ...). These matter
              disproportionately here: the glycan shield is a defining feature of
              viral antigens, and antibodies must either avoid it or bind it.

MODRES records are used as the authoritative per-structure declaration, unioned
with a curated code list so that structures lacking MODRES are still covered.

The honest expected outcome is a LOW frequency. If so, the deliverable is a
well-powered null with a stated detection limit, not a grammar -- and that is
reported as such rather than dressed up.

Outputs
  interactions/analysis/modified_contacts.csv
  metadata/analysis/modified.json
  figures/phase9_modified.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
import time
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

CUTOFF = 4.5

PTM_CODES = {
    "SEP": ("phospho-Ser", "S"), "TPO": ("phospho-Thr", "T"),
    "PTR": ("phospho-Tyr", "Y"), "MLY": ("methyl-Lys", "K"),
    "M3L": ("trimethyl-Lys", "K"), "MLZ": ("methyl-Lys", "K"),
    "ALY": ("acetyl-Lys", "K"), "KCX": ("carboxy-Lys", "K"),
    "CSO": ("oxidised-Cys", "C"), "CSD": ("oxidised-Cys", "C"),
    "OCS": ("oxidised-Cys", "C"), "CME": ("modified-Cys", "C"),
    "CSS": ("modified-Cys", "C"), "SNC": ("S-nitroso-Cys", "C"),
    "HYP": ("hydroxy-Pro", "P"), "HY3": ("hydroxy-Pro", "P"),
    "CGU": ("carboxy-Glu", "E"), "TYS": ("sulfo-Tyr", "Y"),
    "CIR": ("citrulline", "R"), "PCA": ("pyroglutamate", "Q"),
    "NEP": ("phospho-His", "H"), "HIC": ("methyl-His", "H"),
    "SEC": ("selenocysteine", "C"), "MHO": ("oxidised-Met", "M"),
    "FME": ("formyl-Met", "M"), "DDE": ("modified-His", "H"),
    "AAR": ("modified-Arg", "R"), "MEA": ("methyl-Phe", "F"),
    "SAC": ("acetyl-Ser", "S"), "TPQ": ("modified-Tyr", "Y"),
}
# MSE is a crystallographic Met substitute, deliberately NOT counted as a PTM
EXCLUDE = {"MSE"}

GLYCANS = {"NAG", "NDG", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA", "FUC", "FUL",
           "SIA", "NGA", "XYS", "RAM", "A2G", "GCU", "IDS", "SGN"}


def modres_from_file(path: pathlib.Path) -> dict:
    """{resname -> parent} declared by the structure itself (PDB MODRES records)."""
    out = {}
    if path.suffix != ".pdb":
        return out
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("MODRES"):
            name = line[12:15].strip()
            parent = line[24:27].strip()
            if name and name not in EXCLUDE:
                out[name] = parent
        elif line.startswith("ATOM"):
            break                       # MODRES precedes coordinates
    return out


def targets():
    """{arm -> {pdb -> (probe_chains, target_chains)}}"""
    out = {arm: defaultdict(lambda: [set(), set()]) for arm in C.ARMS}
    for arm in C.ARMS:
        for r in C.dedup_rows(arm, C.load_contacts(arm)):
            e = out[arm][r["pdb_id"]]
            e[0].add(r["ab_chain"]); e[1].add(r["ag_chain"])
    return out


_pdb = PDBParser(QUIET=True)
_cif = MMCIFParser(QUIET=True)


def scan(arm, pdb_id, probe_chains, target_chains):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return [], {}
    model = (_cif if path.suffix == ".cif" else _pdb).get_structure(
        pdb_id, str(path))[0]
    declared = modres_from_file(path)

    both = set(probe_chains) | set(target_chains)
    # every heavy atom of the two interacting sides, standard residues included
    side_atoms = {"probe": [], "target": []}
    special = []          # (side, chain, resi, resname, kind, parent)
    for ch in model:
        side = ("probe" if ch.id in probe_chains
                else "target" if ch.id in target_chains else None)
        for res in ch:
            name = res.get_resname().strip()
            het = res.get_id()[0]
            is_std = name in C.AA3_TO_1 and name not in EXCLUDE
            # A standard residue is never a PTM, even when MODRES names it: PDB
            # declares N-glycosylation SITES as "MODRES ... ASN ... ASN", which
            # would otherwise make asparagine the commonest modification in the
            # corpus. The attached sugar is counted separately, as a glycan.
            is_ptm = (not is_std and name not in EXCLUDE
                      and (name in PTM_CODES or name in declared))
            is_gly = name in GLYCANS
            if side and (is_std or is_ptm):
                for at in res:
                    if at.element != "H":
                        side_atoms[side].append((at, ch.id, res))
            if is_ptm or is_gly:
                kind = "glycan" if is_gly else "ptm"
                parent = (PTM_CODES.get(name, (None, declared.get(name)))[1]
                          if kind == "ptm" else None)
                label = (PTM_CODES.get(name, (name, None))[0] if kind == "ptm"
                         else name)
                special.append((side, ch.id, res, name, kind, label, parent, het))

    if not special or not side_atoms["probe"] or not side_atoms["target"]:
        return [], {"n_special": len(special)}

    ns = {s: NeighborSearch([a[0] for a in v]) for s, v in side_atoms.items() if v}
    rows = []
    for side, ch_id, res, name, kind, label, parent, het in special:
        # which side is it near? a glycan may hang off either partner
        best = {}
        for s, tree in ns.items():
            n_contacts = 0
            for at in res:
                if at.element == "H":
                    continue
                hits = tree.search(at.coord, CUTOFF)
                n_contacts += sum(1 for h in hits
                                  if h.get_parent() is not res)
            best[s] = n_contacts
        # interfacial = contacts the OPPOSITE side from where it sits
        opp = ("target" if side == "probe" else "probe") if side else None
        n_opp = best.get(opp, 0) if opp else max(best.values())
        rows.append({
            "arm": arm, "pdb_id": pdb_id, "chain": ch_id,
            "resi": f"{res.get_id()[1]}{res.get_id()[2].strip()}",
            "resname": name, "kind": kind, "label": label,
            "parent_aa": parent or "", "side": side or "other",
            "n_cross_contacts": n_opp,
            "at_interface": int(n_opp > 0),
        })
    return rows, {"n_special": len(special)}


if __name__ == "__main__":
    tg = targets()
    plan = [(arm, p, v[0], v[1]) for arm in C.ARMS for p, v in tg[arm].items()]
    print(f"scanning {len(plan):,} structures for PTMs and glycans", flush=True)

    out_rows, n_fail = [], 0
    t0 = time.time()
    for k, (arm, pdb, pc, tc) in enumerate(plan, 1):
        try:
            rows, _ = scan(arm, pdb, pc, tc)
            out_rows += rows
        except Exception:                                        # noqa: BLE001
            n_fail += 1
        if k % 200 == 0:
            el = time.time() - t0
            print(f"  {k}/{len(plan)}  {el:5.0f}s  "
                  f"eta {el/k*(len(plan)-k)/60:4.1f} min", flush=True)

    path = C.OUT_INT / "modified_contacts.csv"
    if out_rows:
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
            w.writeheader(); w.writerows(out_rows)

    summary = {"n_structures_scanned": len(plan), "n_failed": n_fail}
    for arm in C.ARMS:
        a = [r for r in out_rows if r["arm"] == arm]
        n_struct = len({r["pdb_id"] for r in a})
        for kind in ("ptm", "glycan"):
            k_rows = [r for r in a if r["kind"] == kind]
            iface = [r for r in k_rows if r["at_interface"]]
            n_tot = len({(arm, p) for arm, p in
                         [(r["arm"], r["pdb_id"]) for r in a]}) or 1
            summary[f"{arm}_{kind}"] = {
                "n_residues_found": len(k_rows),
                "n_at_interface": len(iface),
                "n_structures_with_any": len({r["pdb_id"] for r in k_rows}),
                "n_structures_with_interfacial": len({r["pdb_id"] for r in iface}),
                "frac_structures_with_interfacial": round(
                    len({r["pdb_id"] for r in iface})
                    / len(tg[arm]), 4),
                "top_types": Counter(r["label"] for r in iface).most_common(10),
            }
    # detection limit: with N structures and k observations, the upper 95% bound
    for arm in C.ARMS:
        n = len(tg[arm])
        k = summary[f"{arm}_ptm"]["n_structures_with_interfacial"]
        summary[f"{arm}_ptm"]["detection_note"] = (
            f"{k} of {n} structures carry an interfacial PTM. With this corpus a "
            f"per-structure rate below ~{3/n:.4f} (rule of three) could not have "
            f"been distinguished from zero.")
    summary["method_note"] = (
        "MODRES records are unioned with a curated PTM code list so structures "
        "without MODRES are still covered. MSE (selenomethionine) is EXCLUDED: it "
        "is a crystallographic phasing substitute for Met, not a biological "
        "modification, and counting it would inflate every PTM statistic.")
    (C.OUT_MET / "modified.json").write_text(json.dumps(summary, indent=1))

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, kind, title in zip(axes, ("ptm", "glycan"),
                               ("Post-translational modifications",
                                "Glycans")):
        labels, vals = [], []
        for arm in C.ARMS:
            d = summary[f"{arm}_{kind}"]
            labels.append(C.ARM_LABEL[arm])
            vals.append(d["frac_structures_with_interfacial"])
        ax.bar(range(2), vals, color=["#b2182b", "#2166ac"], width=0.5)
        for i, v in enumerate(vals):
            ax.text(i, v + 0.005, f"{v:.1%}", ha="center", fontsize=10)
        ax.set_xticks(range(2)); ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("fraction of structures with one at the interface")
        top = summary[f"antibody_{kind}"]["top_types"][:5]
        ax.set_title(f"{title}\ncommonest (antibody): " +
                     (", ".join(f"{t}" for t, _ in top) if top else "none"),
                     fontsize=9)
    fig.suptitle("Phase 9 - modified residues and glycans at the interface",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(C.OUT_FIG / "phase9_modified.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps({k: v for k, v in summary.items()
                      if not k.endswith("_note")}, indent=1)[:2600])
