"""Define what counts as an INDEPENDENT observation, from sequence rather than name strings.

Until now a "target" was the COMPND MOLECULE string. That over-counts independence: the
same protein deposited under two names is two targets, and close homologs are treated as
unrelated. Every statistic that treats targets as independent (bootstrap CI, per-target
replication, jackknife, saturation curve) is optimistic as a result.

This builds, per PDB entry:
  uniprot        accessions of all protein entities (identity of the protein)
  cluster95/50/30  RCSB precomputed sequence-identity cluster of the LARGEST protein
                 entity — cluster30 merges homologs and is the conservative unit
  na_cluster     aptamer/nucleic-acid sequence cluster (single-linkage at >=80% identity)
  resolution, method, r_free  for quality filtering

Run:  python3 code/build_independence_units.py [aptamer|natural|both]
Out:  metadata/independence_units_{corpus}.json
      interactions/residue_contacts_{corpus}_{level}.csv   (target column regrouped)
"""
from __future__ import annotations

# --- layout bootstrap -------------------------------------------------------------
# Scripts live in labelled subdirectories of code/ but keep a flat import namespace,
# and every path below is anchored at the heatmap_validation root (_HV), never at the
# script's own depth.
import sys as _sys, pathlib as _pl                                        # noqa: E402
_CODE = _pl.Path(__file__).resolve().parents[1]
_HV = _CODE.parent
_sys.path[:0] = [str(_CODE)] + [str(_p) for _p in sorted(_CODE.iterdir())
                                if _p.is_dir() and not _p.name.startswith(("_", "."))]
# ----------------------------------------------------------------------------------

import csv
import json
import pathlib
import sys
from collections import defaultdict
from difflib import SequenceMatcher

HV = _HV
MET = HV / "metadata"
INT = HV / "interactions"
LEVELS = ["uniprot", "cluster95", "cluster50", "cluster30"]


def largest_protein(entities):
    best, blen = None, -1
    for e in entities:
        if "polypeptide" in e["type"] and (e["length"] or 0) > blen:
            best, blen = e, e["length"] or 0
    return best


def na_clusters(meta, threshold=0.80):
    """Single-linkage clustering of nucleic-acid entity sequences by identity ratio."""
    seqs = {}
    for pid, m in meta.items():
        na = [e for e in m["entities"] if "ribonucleotide" in e["type"]]
        if na:
            seqs[pid] = max(na, key=lambda e: e["length"] or 0)["sequence"]
    ids = sorted(seqs)
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            s1, s2 = seqs[ids[a]], seqs[ids[b]]
            if not s1 or not s2:
                continue
            if min(len(s1), len(s2)) / max(len(s1), len(s2)) < 0.6:
                continue
            if SequenceMatcher(None, s1, s2).ratio() >= threshold:
                ra, rb = find(ids[a]), find(ids[b])
                if ra != rb:
                    parent[ra] = rb
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)
    lookup = {}
    for k, (root, members) in enumerate(sorted(groups.items()), 1):
        for m in members:
            lookup[m] = f"na{k:03d}"
    return lookup


def build(corpus):
    meta = json.loads((MET / f"entity_metadata_{corpus}.json").read_text())
    nac = na_clusters(meta)
    units = {}
    for pid, m in meta.items():
        prots = [e for e in m["entities"] if "polypeptide" in e["type"]]
        big = largest_protein(m["entities"])
        uni = sorted({u for e in prots for u in e["uniprot"]})
        cl = (big or {}).get("clusters", {}) or {}
        units[pid] = dict(
            uniprot="|".join(uni) if uni else f"noUP:{(big or {}).get('desc','?')[:28].lower()}",
            cluster95=f"c95_{cl.get('95')}" if cl.get("95") else None,
            cluster50=f"c50_{cl.get('50')}" if cl.get("50") else None,
            cluster30=f"c30_{cl.get('30')}" if cl.get("30") else None,
            na_cluster=nac.get(pid),
            n_protein_entities=len(prots),
            largest_protein=(big or {}).get("desc", ""),
            protein_length=(big or {}).get("length"),
            organism=(big or {}).get("organism", []),
            method=m["method"], resolution=m["resolution"], r_free=m["r_free"])
    # entries with no cluster id (short peptides are not clustered by RCSB): fall back
    for pid, u in units.items():
        for lvl in ("cluster95", "cluster50", "cluster30"):
            if not u[lvl]:
                u[lvl] = f"nc:{u['uniprot']}"
    (MET / f"independence_units_{corpus}.json").write_text(json.dumps(units, indent=2))

    src = INT / ("residue_contacts.csv" if corpus == "aptamer"
                 else "residue_contacts_natural.csv")
    rows = list(csv.DictReader(open(src)))
    counts = {}
    for lvl in LEVELS:
        out = INT / f"residue_contacts_{corpus}_{lvl}.csv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            for r in rows:
                u = units.get(r["pdb_id"])
                if not u:
                    continue
                r2 = dict(r); r2["target"] = u[lvl]
                w.writerow(r2)
        counts[lvl] = len({units[r["pdb_id"]][lvl] for r in rows if r["pdb_id"] in units})
    counts["compnd_string(old)"] = len({r["target"] for r in rows})
    counts["na_cluster"] = len(set(nac.values()))
    print(f"[{corpus}] {len(units)} entries — independent units by level:")
    for k, v in counts.items():
        print(f"    {k:22s} {v}")
    with_res = [u["resolution"] for u in units.values() if u["resolution"]]
    print(f"    resolution available for {len(with_res)}/{len(units)}, "
          f"median {sorted(with_res)[len(with_res)//2]:.2f} A" if with_res else "")
    return counts


def main(which="both"):
    for c in (["aptamer", "natural"] if which == "both" else [which]):
        build(c)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "both")
