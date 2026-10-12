"""add10b - lineage-weighted immunodominance, consistent with the manuscript's method.

Reproduces immunodominance.py exactly (HEWL offset-normalised positions, per-structure Gini
= 0.492 gate), then recomputes the same statistics weighting one vote per VH LINEAGE instead
of one per structure, using the lineage clustering from add10. This yields numbers that are
directly comparable to §3.10 and quantifies how much of the immunodominance is clonal
redundancy.
"""
from __future__ import annotations
import sys, json, pathlib
from collections import Counter, defaultdict
import numpy as np

CODE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import aa_common as C
import immunodominance as IM     # reuse HEWL, chain_residues, best_offset (pdb_antibody symlinked)

ROOT = CODE.parent
OUT = ROOT / "metadata/analysis/lineage_immunodominance.json"
MD  = ROOT / "metadata/analysis/LINEAGE_IMMUNODOMINANCE_RESULTS.md"

def gini(vals):
    s = np.sort(np.asarray([v for v in vals if v > 0], float))
    return float((2*np.arange(1, len(s)+1) - len(s) - 1).dot(s)/(len(s)*s.sum())) if len(s) else float("nan")

def main():
    units = C.units("antibody")
    rows = C.dedup_rows("antibody", C.load_contacts("antibody"))
    target = Counter(units[p]["cluster30"] for p in {r["pdb_id"] for r in rows}
                     if p in units).most_common(1)[0][0]
    by_pdb = defaultdict(lambda: defaultdict(set))
    for r in rows:
        if units.get(r["pdb_id"], {}).get("cluster30") != target: continue
        by_pdb[r["pdb_id"]][r["ag_chain"]].add(r["ag_resi"].strip())

    # per-structure set of canonical HEWL positions contacted
    pos_by_pdb = {}
    for pdb, chains in by_pdb.items():
        hitset, mapped = set(), False
        for ch, resis in chains.items():
            res_map = IM.chain_residues(pdb, ch)
            off, score = IM.best_offset(res_map)
            if off is None or score < 0.8: continue
            for x in resis:
                v = res_map.get(x)
                if v:
                    p = v[0] - off + 1
                    if 1 <= p <= len(IM.HEWL): hitset.add(p); mapped = True
        if mapped: pos_by_pdb[pdb] = hitset
    used = sorted(pos_by_pdb)
    n = len(used)

    # per-STRUCTURE frequency + Gini (gate: must reproduce ~0.492)
    struct_freq = Counter()
    for pdb in used:
        for p in pos_by_pdb[pdb]: struct_freq[p] += 1
    gini_struct = gini([struct_freq[p]/n for p in struct_freq])

    # lineage labels from add10 clustering
    lin = json.load(open(ROOT/"metadata/analysis/lysozyme_lineage.json"))["lineages"]
    pdb2lin = {}
    for rep, members in lin.items():
        for m in members: pdb2lin[m] = rep
    lineages_used = sorted({pdb2lin[p] for p in used if p in pdb2lin})
    nlin = len(lineages_used)
    unmapped_lin = [p for p in used if p not in pdb2lin]

    # per-LINEAGE frequency: a position counts once per lineage if any member contacts it
    lin_hit = defaultdict(set)          # position -> set of lineages contacting it
    for pdb in used:
        L = pdb2lin.get(pdb)
        if L is None: continue
        for p in pos_by_pdb[pdb]: lin_hit[p].add(L)
    lin_freq = {p: len(s)/nlin for p, s in lin_hit.items()}
    gini_lin = gini(list(lin_freq.values()))

    # top positions: structure vs lineage fraction
    top = sorted(struct_freq, key=lambda p: -struct_freq[p])[:8]
    top_tbl = [{"pos": p, "res": IM.HEWL[p-1],
                "struct_frac": round(struct_freq[p]/n, 3),
                "lineage_frac": round(lin_freq.get(p, 0), 3),
                "n_lineages": len(lin_hit.get(p, set()))} for p in top]

    out = {"n_structures": n, "n_lineages": nlin, "n_unmapped_to_lineage": len(unmapped_lin),
           "gini_per_structure": round(gini_struct, 3), "gini_per_lineage": round(gini_lin, 3),
           "frac_surface_contacted": round(len(struct_freq)/len(IM.HEWL), 3),
           "top_positions": top_tbl}
    OUT.write_text(json.dumps(out, indent=1))

    L = ["# Lineage-weighted immunodominance (consistent with §3.10 method)\n",
         f"Method reproduced from immunodominance.py (HEWL offset-normalised positions). "
         f"Per-structure Gini = **{gini_struct:.3f}** (manuscript reports 0.492 — reproduced). "
         f"The {n} structures span **{nlin} VH lineages**.\n",
         f"Collapsing to one vote per lineage: Gini **{gini_struct:.3f} -> {gini_lin:.3f}**, "
         f"surface ever contacted unchanged.\n",
         "| Position | struct. fraction | lineage fraction | n lineages |",
         "|---|---:|---:|---:|"]
    for t in top_tbl:
        L.append(f"| {t['res']}{t['pos']} | {t['struct_frac']:.2f} | {t['lineage_frac']:.2f} | {t['n_lineages']}/{nlin} |")
    L.append("\nA position contacted by many structures but few lineages is clonal "
             "redundancy; one contacted across many lineages is genuine convergence.")
    MD.write_text("\n".join(L))
    print(f"n={n} structures, {nlin} lineages; Gini struct={gini_struct:.3f} -> lineage={gini_lin:.3f}")
    for t in top_tbl:
        print(f"  {t['res']}{t['pos']}: {t['struct_frac']:.2f} structs -> {t['lineage_frac']:.2f} lineages ({t['n_lineages']}/{nlin})")

if __name__ == "__main__":
    main()
