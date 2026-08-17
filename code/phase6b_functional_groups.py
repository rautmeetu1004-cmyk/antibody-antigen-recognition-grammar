"""Phase 6b - functional-group chemistry: what actually touches what.

The 20x20 residue matrix blurs chemistry that the atoms keep separate. Asp and Glu
are both "acidic" but present their carboxylate at different distances from the
backbone; Tyr is scored as "aromatic" though its hydroxyl is a strong H-bond donor
and its ring is a pi surface, and those two moieties do different jobs. This phase
drops to the atom tables (211k antibody / 4.49M PPI atom-atom contacts) and asks
the chemically honest question: which FUNCTIONAL GROUP contacts which.

Groups are defined by (residue, atom name) so every heavy atom lands in exactly one
group. Backbone is split into its amide N-H and carbonyl C=O rather than lumped,
because those are a donor and an acceptor and behave differently.

Outputs
  figures/phase6b_functional_groups.png
  interactions/analysis/fgroup_matrix_{arm}.csv
  metadata/analysis/functional_groups.json
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
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

# ---------------------------------------------------------------- group definitions
BACKBONE = {"N": "bb amide N-H", "O": "bb carbonyl C=O", "C": "bb carbonyl C=O",
            "CA": "aliphatic C", "OXT": "bb carbonyl C=O"}

SIDE: dict[str, dict[str, str]] = {
    "ARG": {"NE": "guanidinium", "NH1": "guanidinium", "NH2": "guanidinium",
            "CZ": "guanidinium", "CB": "aliphatic C", "CG": "aliphatic C",
            "CD": "aliphatic C"},
    "LYS": {"NZ": "ammonium", "CB": "aliphatic C", "CG": "aliphatic C",
            "CD": "aliphatic C", "CE": "aliphatic C"},
    "ASP": {"OD1": "carboxylate", "OD2": "carboxylate", "CG": "carboxylate",
            "CB": "aliphatic C"},
    "GLU": {"OE1": "carboxylate", "OE2": "carboxylate", "CD": "carboxylate",
            "CB": "aliphatic C", "CG": "aliphatic C"},
    "ASN": {"OD1": "amide (Asn/Gln)", "ND2": "amide (Asn/Gln)",
            "CG": "amide (Asn/Gln)", "CB": "aliphatic C"},
    "GLN": {"OE1": "amide (Asn/Gln)", "NE2": "amide (Asn/Gln)",
            "CD": "amide (Asn/Gln)", "CB": "aliphatic C", "CG": "aliphatic C"},
    "SER": {"OG": "hydroxyl", "CB": "aliphatic C"},
    "THR": {"OG1": "hydroxyl", "CB": "aliphatic C", "CG2": "aliphatic C"},
    "TYR": {"OH": "phenol -OH", "CZ": "aromatic ring (Phe/Tyr)",
            "CG": "aromatic ring (Phe/Tyr)", "CD1": "aromatic ring (Phe/Tyr)",
            "CD2": "aromatic ring (Phe/Tyr)", "CE1": "aromatic ring (Phe/Tyr)",
            "CE2": "aromatic ring (Phe/Tyr)", "CB": "aliphatic C"},
    "PHE": {"CZ": "aromatic ring (Phe/Tyr)", "CG": "aromatic ring (Phe/Tyr)",
            "CD1": "aromatic ring (Phe/Tyr)", "CD2": "aromatic ring (Phe/Tyr)",
            "CE1": "aromatic ring (Phe/Tyr)", "CE2": "aromatic ring (Phe/Tyr)",
            "CB": "aliphatic C"},
    "TRP": {"NE1": "indole N-H", "CG": "indole ring", "CD1": "indole ring",
            "CD2": "indole ring", "CE2": "indole ring", "CE3": "indole ring",
            "CZ2": "indole ring", "CZ3": "indole ring", "CH2": "indole ring",
            "CB": "aliphatic C"},
    "HIS": {"ND1": "imidazole", "NE2": "imidazole", "CG": "imidazole",
            "CD2": "imidazole", "CE1": "imidazole", "CB": "aliphatic C"},
    "CYS": {"SG": "thiol/thioether", "CB": "aliphatic C"},
    "MET": {"SD": "thiol/thioether", "CB": "aliphatic C", "CG": "aliphatic C",
            "CE": "aliphatic C"},
    "PRO": {"CB": "aliphatic C", "CG": "aliphatic C", "CD": "aliphatic C"},
}
for _r in ("ALA", "VAL", "LEU", "ILE", "GLY"):
    SIDE[_r] = {}
SIDE["MSE"] = SIDE["MET"] | {"SE": "thiol/thioether"}

GROUPS = ["guanidinium", "ammonium", "carboxylate", "amide (Asn/Gln)", "hydroxyl",
          "phenol -OH", "aromatic ring (Phe/Tyr)", "indole N-H", "indole ring",
          "imidazole", "thiol/thioether", "aliphatic C",
          "bb amide N-H", "bb carbonyl C=O"]
GIDX = {g: i for i, g in enumerate(GROUPS)}

GROUP_KIND = {
    "guanidinium": "cationic", "ammonium": "cationic",
    "carboxylate": "anionic",
    "amide (Asn/Gln)": "polar", "hydroxyl": "polar", "phenol -OH": "polar",
    "indole N-H": "polar", "imidazole": "polar",
    "aromatic ring (Phe/Tyr)": "aromatic", "indole ring": "aromatic",
    "thiol/thioether": "sulfur", "aliphatic C": "apolar",
    "bb amide N-H": "backbone", "bb carbonyl C=O": "backbone",
}
KIND_COLOR = {"cationic": "#2166ac", "anionic": "#b2182b", "polar": "#1b7837",
              "aromatic": "#762a83", "apolar": "#b8860b", "sulfur": "#e08214",
              "backbone": "#4d4d4d"}


def group_of(res3: str, atom: str) -> str | None:
    if atom in BACKBONE:
        return BACKBONE[atom]
    return SIDE.get(res3, {}).get(atom)


def build(arm: str, level: str = "cluster30"):
    """Weighted group-x-group matrices, overall and per bond type."""
    units = C.units(arm)
    unit_of: dict[tuple, str] = {}
    structs = defaultdict(set)
    path = C.INT / f"{arm}_atom_contacts.csv"
    with path.open() as fh:
        for row in csv.DictReader(fh):
            if arm == "ppi" and row.get("direction") != "a_to_b":
                continue                      # dedup: one row per physical contact
            u = units.get(row["pdb_id"])
            if u is None:
                continue
            if arm == "ppi":
                u = u.get(row["direction"])
                if u is None:
                    continue
            unit = u.get(level)
            if not unit or unit == "NA":
                continue
            key = (row["pdb_id"], row.get("direction", ""))
            unit_of[key] = unit
            structs[unit].add(row["pdb_id"])

    n = len(GROUPS)
    M = np.zeros((n, n))
    Mb = {b: np.zeros((n, n)) for b in ("salt", "hbond", "vdw")}
    raw = np.zeros((n, n))
    n_rows = n_mapped = 0

    with path.open() as fh:
        for row in csv.DictReader(fh):
            if arm == "ppi" and row.get("direction") != "a_to_b":
                continue
            n_rows += 1
            key = (row["pdb_id"], row.get("direction", ""))
            unit = unit_of.get(key)
            if unit is None:
                continue
            ga = group_of(row["ab_aa3"], row["ab_atom"])
            gb = group_of(row["ag_aa3"], row["ag_atom"])
            if ga is None or gb is None:
                continue
            n_mapped += 1
            w = 1.0 / len(structs[unit])
            i, j = GIDX[ga], GIDX[gb]
            M[i, j] += w
            raw[i, j] += 1
            for b in Mb:
                if row.get(b) == "1":
                    Mb[b][i, j] += w

    if arm == "ppi":                          # co-evolved partners are exchangeable
        M = (M + M.T) / 2.0
        raw = (raw + raw.T) / 2.0
        Mb = {b: (v + v.T) / 2.0 for b, v in Mb.items()}
    return {"M": M, "raw": raw, "bonds": Mb,
            "n_rows": n_rows, "n_mapped": n_mapped,
            "coverage": n_mapped / max(1, n_rows)}


def draw(ax, lo2, title, show_y=True):
    n = len(GROUPS)
    im = ax.imshow(np.ma.masked_invalid(lo2), cmap="RdBu_r",
                   norm=TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5))
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(GROUPS, fontsize=6.5, rotation=45, ha="right")
    ax.set_yticklabels(GROUPS if show_y else [""] * n, fontsize=6.5)
    for t, g in zip(ax.get_xticklabels(), GROUPS):
        t.set_color(KIND_COLOR[GROUP_KIND[g]])
    if show_y:
        for t, g in zip(ax.get_yticklabels(), GROUPS):
            t.set_color(KIND_COLOR[GROUP_KIND[g]])
    for i in range(n):
        for j in range(n):
            if np.isfinite(lo2[i, j]) and abs(lo2[i, j]) >= 0.6:
                ax.text(j, i, f"{lo2[i,j]:+.1f}", ha="center", va="center",
                        fontsize=5.2,
                        color="w" if abs(lo2[i, j]) > 1.0 else "k")
    ax.set_title(title, fontsize=10)
    return im


def top_cells(lo2, M, n=12, enriched=True, min_count=200):
    cells = []
    for i in range(len(GROUPS)):
        for j in range(len(GROUPS)):
            if np.isfinite(lo2[i, j]) and M[i, j] >= min_count:
                cells.append((lo2[i, j], GROUPS[i], GROUPS[j]))
    cells.sort(reverse=enriched)
    return [{"pair": f"{a} <-> {b}", "log2_or": round(float(v), 3)}
            for v, a, b in cells[:n]]


if __name__ == "__main__":
    data = {arm: build(arm) for arm in C.ARMS}

    fig, axes = plt.subplots(1, 3, figsize=(21, 6.6))
    out = {}
    lo2s = {}
    for k, arm in enumerate(C.ARMS):
        d = data[arm]
        lo2 = C.log2_odds(d["M"])
        lo2s[arm] = lo2
        im = draw(axes[k], lo2, f"{C.ARM_LABEL[arm]}\n"
                                f"{int(d['raw'].sum()):,} atom-atom contacts mapped",
                  show_y=(k == 0))
        out[arm] = {
            "coverage_of_atom_rows": round(d["coverage"], 4),
            "n_atom_contacts_mapped": int(d["raw"].sum()),
            "top_enriched": top_cells(lo2, d["raw"], 12, True),
            "top_depleted": top_cells(lo2, d["raw"], 12, False),
        }
        with (C.OUT_INT / f"fgroup_matrix_{arm}.csv").open("w") as fh:
            fh.write("group_a,group_b,log2_or,weighted,raw\n")
            for i, a in enumerate(GROUPS):
                for j, b in enumerate(GROUPS):
                    fh.write(f'"{a}","{b}",{lo2[i,j]:.4f},{d["M"][i,j]:.4f},'
                             f'{int(d["raw"][i,j])}\n')

    diff = lo2s["antibody"] - lo2s["ppi"]
    imd = axes[2].imshow(np.ma.masked_invalid(diff), cmap="PuOr_r",
                         norm=TwoSlopeNorm(vmin=-1.2, vcenter=0, vmax=1.2))
    axes[2].set_xticks(range(len(GROUPS))); axes[2].set_yticks(range(len(GROUPS)))
    axes[2].set_xticklabels(GROUPS, fontsize=6.5, rotation=45, ha="right")
    axes[2].set_yticklabels([""] * len(GROUPS))
    for t, g in zip(axes[2].get_xticklabels(), GROUPS):
        t.set_color(KIND_COLOR[GROUP_KIND[g]])
    for i in range(len(GROUPS)):
        for j in range(len(GROUPS)):
            if np.isfinite(diff[i, j]) and abs(diff[i, j]) >= 0.4:
                axes[2].text(j, i, f"{diff[i,j]:+.1f}", ha="center", va="center",
                             fontsize=5.2,
                             color="w" if abs(diff[i, j]) > 0.9 else "k")
    axes[2].set_title("Difference (antibody - PPI)\n"
                      "orange = antibody-enriched, purple = PPI-enriched", fontsize=10)

    fig.colorbar(im, ax=axes[:2], fraction=0.02, pad=0.02).set_label(
        "log$_2$(obs/exp)", fontsize=8)
    fig.colorbar(imd, ax=axes[2], fraction=0.04, pad=0.02).set_label(
        "$\\Delta$ log$_2$(obs/exp)", fontsize=8)
    fig.suptitle("Phase 6b - functional-group recognition chemistry", fontsize=13)
    fig.savefig(C.OUT_FIG / "phase6b_functional_groups.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    out["similarity_antibody_vs_ppi"] = C.matrix_similarity(lo2s["antibody"],
                                                            lo2s["ppi"])
    flat = [(diff[i, j], GROUPS[i], GROUPS[j])
            for i in range(len(GROUPS)) for j in range(len(GROUPS))
            if np.isfinite(diff[i, j]) and data["antibody"]["raw"][i, j] >= 200]
    flat.sort()
    out["most_ppi_enriched"] = [{"pair": f"{a} <-> {b}", "delta": round(float(v), 3)}
                                for v, a, b in flat[:10]]
    out["most_antibody_enriched"] = [{"pair": f"{a} <-> {b}",
                                      "delta": round(float(v), 3)}
                                     for v, a, b in flat[::-1][:10]]
    (C.OUT_MET / "functional_groups.json").write_text(json.dumps(out, indent=1))

    for arm in C.ARMS:
        print(f"{arm}: {out[arm]['n_atom_contacts_mapped']:,} atom contacts, "
              f"coverage {out[arm]['coverage_of_atom_rows']:.3f}")
    print(f"\nfunctional-group map similarity r = "
          f"{out['similarity_antibody_vs_ppi']:.3f}")
    print("\nMost antibody-enriched group pairs:")
    for r in out["most_antibody_enriched"][:8]:
        print(f"   {r['delta']:+.2f}  {r['pair']}")
    print("\nMost PPI-enriched group pairs:")
    for r in out["most_ppi_enriched"][:8]:
        print(f"   {r['delta']:+.2f}  {r['pair']}")
