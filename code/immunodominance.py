"""Immunodominance: where do 45 independent antibodies actually hit lysozyme?

Every other analysis in this study spends effort SUPPRESSING corpus redundancy --
45 lysozyme complexes are down-weighted to one cluster30 unit so they cannot
dominate the grammar. This analysis inverts that: the redundancy is the data. Forty-
five independently raised antibodies against one antigen is an experiment in which
surface positions the immune system finds, and how unevenly.

The numbering problem, and why it cannot be ignored. Antigen residue numbers in the
corpus run from 10 to 1584 for what is a 129-residue protein: different depositions
use different offsets (chain-concatenated constructs, 600/1000-series numbering).
Pooling those naively would smear every epitope. Each antigen chain is therefore
re-aligned to the canonical hen egg-white lysozyme sequence by finding the integer
offset that maximises sequence identity, and chains that do not align are dropped
and counted rather than silently included.

Outputs
  metadata/analysis/immunodominance.json
  figures/immunodominance_profile.png
  figures/pymol_immunodominance.png   (rendered separately by PyMOL)
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
from Bio.PDB import MMCIFParser, PDBParser

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

# canonical hen egg-white lysozyme (P00698 mature chain, 129 aa)
HEWL = ("KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRW"
        "WCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")
_pdb = PDBParser(QUIET=True)
_cif = MMCIFParser(QUIET=True)


def chain_residues(pdb_id, chain):
    d = C.ROOT / "pdb_antibody"
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return {}
    model = (_cif if path.suffix == ".cif" else _pdb).get_structure(
        pdb_id, str(path))[0]
    if chain not in [c.id for c in model]:
        return {}
    out = {}
    for r in model[chain]:
        aa = C.AA3_TO_1.get(r.get_resname().strip())
        if aa:
            out[f"{r.get_id()[1]}{r.get_id()[2].strip()}"] = (r.get_id()[1], aa)
    return out


def best_offset(res_map):
    """Integer offset o such that (resnum - o) indexes into HEWL, maximising identity."""
    if not res_map:
        return None, 0.0
    nums = [v[0] for v in res_map.values()]
    best, best_score = None, 0.0
    lo, hi = min(nums), max(nums)
    for o in range(lo - 2, lo + 3):                       # numbering starts near 1
        hits = tot = 0
        for num, aa in res_map.values():
            i = num - o
            if 0 <= i < len(HEWL):
                tot += 1
                hits += int(HEWL[i] == aa)
        if tot >= 50 and hits / tot > best_score:
            best, best_score = o, hits / tot
    # also try offsets implied by 100/600/1000-series numbering
    for base in (0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1400, 1500):
        o = lo - 1 - base + base
        for cand in (lo - 1, lo - 1 + base):
            hits = tot = 0
            for num, aa in res_map.values():
                i = num - cand
                if 0 <= i < len(HEWL):
                    tot += 1
                    hits += int(HEWL[i] == aa)
            if tot >= 50 and hits / tot > best_score:
                best, best_score = cand, hits / tot
    return best, best_score


if __name__ == "__main__":
    units = C.units("antibody")
    rows = C.dedup_rows("antibody", C.load_contacts("antibody"))
    # the most-SOLVED antigen, i.e. most distinct structures -- counting contacts
    # instead would pick whichever antigen is largest, not the one with the most
    # independent antibodies raised against it
    target_cluster = Counter(
        units[p]["cluster30"] for p in {r["pdb_id"] for r in rows}
        if p in units).most_common(1)[0][0]

    by_pdb = defaultdict(lambda: defaultdict(set))
    for r in rows:
        if units.get(r["pdb_id"], {}).get("cluster30") != target_cluster:
            continue
        by_pdb[r["pdb_id"]][r["ag_chain"]].add(r["ag_resi"].strip())

    hits = Counter()
    used, dropped = [], []
    for pdb, chains in sorted(by_pdb.items()):
        mapped_any = False
        for ch, resis in chains.items():
            res_map = chain_residues(pdb, ch)
            off, score = best_offset(res_map)
            if off is None or score < 0.8:
                continue
            for x in resis:
                v = res_map.get(x)
                if not v:
                    continue
                pos = v[0] - off + 1                       # 1-based canonical
                if 1 <= pos <= len(HEWL):
                    hits[pos] += 1
                    mapped_any = True
        (used if mapped_any else dropped).append(pdb)

    n = len(used)
    prof = np.array([hits.get(i, 0) for i in range(1, len(HEWL) + 1)])
    frac = prof / max(1, n)
    order = np.argsort(-prof)
    top = [{"position": int(i + 1), "residue": HEWL[i],
            "n_antibodies": int(prof[i]), "fraction": round(float(frac[i]), 3)}
           for i in order[:20]]

    # how concentrated is the response? Gini over surface positions touched at all
    touched = prof[prof > 0]
    s = np.sort(touched)
    gini = float((2 * np.arange(1, len(s) + 1) - len(s) - 1).dot(s)
                 / (len(s) * s.sum())) if len(s) else float("nan")

    out = {
        "antigen": units[used[0]]["antigen_description"] if used else None,
        "cluster30": target_cluster,
        "n_structures_used": n,
        "n_structures_dropped_unmappable": len(dropped),
        "dropped": dropped,
        "n_positions_contacted": int((prof > 0).sum()),
        "protein_length": len(HEWL),
        "frac_surface_ever_contacted": round(float((prof > 0).mean()), 3),
        "gini_of_contact_frequency": round(gini, 3),
        "top_positions": top,
        "note": ("Residue numbering was normalised to canonical HEWL by offset "
                 "alignment (>=80% sequence identity required); unmappable chains "
                 "are dropped and counted. Counts are numbers of INDEPENDENT "
                 "antibody complexes contacting that position."),
    }
    (C.OUT_MET / "immunodominance.json").write_text(json.dumps(out, indent=1))

    fig, ax = plt.subplots(figsize=(14, 4.2))
    ax.bar(range(1, len(HEWL) + 1), frac, color="#2166ac", width=0.9)
    for i in order[:8]:
        ax.text(i + 1, frac[i] + 0.015, f"{HEWL[i]}{i+1}", ha="center", fontsize=8)
    ax.set_xlabel("lysozyme residue (canonical numbering)")
    ax.set_ylabel(f"fraction of {n} antibodies contacting")
    ax.set_xlim(0, len(HEWL) + 1)
    ax.set_title(f"Immunodominance: where {n} independently raised antibodies bind "
                 f"one antigen\n"
                 f"{out['frac_surface_ever_contacted']:.0%} of the sequence is "
                 f"contacted by at least one antibody; "
                 f"Gini {out['gini_of_contact_frequency']:.2f}", fontsize=11)
    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "immunodominance_profile.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    print(f"antigen: {out['antigen']}  ({n} structures used, "
          f"{len(dropped)} dropped)")
    print(f"positions contacted: {out['n_positions_contacted']}/{len(HEWL)} "
          f"({out['frac_surface_ever_contacted']:.0%}), Gini "
          f"{out['gini_of_contact_frequency']:.2f}")
    print("most immunodominant positions:")
    for t in top[:10]:
        print(f"   {t['residue']}{t['position']:<4d} {t['n_antibodies']:3d} "
              f"antibodies ({t['fraction']:.0%})")
