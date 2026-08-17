"""Matched targets - the same protein read by an antibody and by a natural partner.

Every corpus-wide comparison in this study carries one irreducible objection: the two
groups bind DIFFERENT proteins, so a chemical difference could be a target
difference. Phases 6c/6d controlled that statistically (secondary structure, CATH
fold). This controls it by construction.

Proteins present as an antibody ANTIGEN in one arm and as a natural PPI PARTNER in
the other are matched on UniProt accession -- not on description string, which would
merge "Green fluorescent protein" with "GFP-like chromoprotein". For each matched
protein we ask two questions the corpus-wide analysis cannot:

  SAME SURFACE?   Do the antibody and the natural partner bind the same patch?
                  Measured as Jaccard overlap of contacted residue positions, after
                  normalising numbering by sequence alignment between the two
                  structures' copies of the protein.
  SAME CHEMISTRY? For the SAME protein, does the antibody read it with different
                  residue-pair preferences than its natural partner does?

This is the direct analogue of the nucleic-acid study's DUX4 case, but with 38
candidate proteins rather than one.

Outputs
  metadata/analysis/matched_targets.json
  figures/case_matched_targets.png
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

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import phase6a_structural as P


def antigen_uniprot():
    """{uniprot -> {'antibody': [(pdb, chain)], 'ppi': [(pdb, chain)]}}"""
    out = defaultdict(lambda: defaultdict(list))
    ab_units = C.units("antibody")
    for r in C.dedup_rows("antibody", C.load_contacts("antibody")):
        u = ab_units.get(r["pdb_id"], {}).get("uniprot")
        if u and u != "NA":
            out[u]["antibody"].append((r["pdb_id"], r["ag_chain"]))
    pp_units = C.units("ppi")
    for r in C.dedup_rows("ppi", C.load_contacts("ppi")):
        d = pp_units.get(r["pdb_id"], {})
        # ppi_units[pdb]["a_to_b"] describes the chain that is the TARGET in the
        # a_to_b direction (its "chains" list matches ag_chain). Reading b_to_a here
        # instead silently matched each protein against its PARTNER, which showed up
        # as ~20% sequence identity between two copies of the same protein.
        u = (d.get("a_to_b") or {}).get("uniprot")
        if u and u != "NA":
            out[u]["ppi"].append((r["pdb_id"], r["ag_chain"]))
    for u in out:
        for k in out[u]:
            out[u][k] = sorted(set(out[u][k]))
    return out


def chain_seq_and_positions(pdb_id, chain, arm):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif")
                 if (d / nm).exists()), None)
    if path is None:
        return "", {}
    try:
        model = P.load_model(path)
    except Exception:                                            # noqa: BLE001
        return "", {}
    if chain not in [c.id for c in model]:
        return "", {}
    seq, keys = [], []
    for r in model[chain]:
        aa = C.AA3_TO_1.get(r.get_resname().strip())
        if aa:
            seq.append(aa)
            keys.append(P.residue_key(r))
    return "".join(seq), {k: i for i, k in enumerate(keys)}


def align_index(seq_a, seq_b):
    """Map index in seq_a -> index in seq_b via global alignment (free end gaps)."""
    from Bio import Align
    al = Align.PairwiseAligner()
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = "global"
    al.end_gap_score = 0.0
    try:
        aln = al.align(seq_a, seq_b)[0]
    except Exception:                                            # noqa: BLE001
        return {}
    m = {}
    ia = ib = 0
    for ca, cb in zip(aln[0], aln[1]):
        if ca != "-" and cb != "-":
            m[ia] = ib
        if ca != "-":
            ia += 1
        if cb != "-":
            ib += 1
    return m


if __name__ == "__main__":
    idx = antigen_uniprot()
    matched = {u: v for u, v in idx.items() if v.get("antibody") and v.get("ppi")}
    print(f"{len(matched)} proteins bound by BOTH an antibody and a natural partner")

    contacts = {"antibody": defaultdict(set), "ppi": defaultdict(set)}
    for arm in C.ARMS:
        for r in C.dedup_rows(arm, C.load_contacts(arm)):
            contacts[arm][(r["pdb_id"], r["ag_chain"])].add(r["ag_resi"].strip())

    rows = []
    for u, v in sorted(matched.items()):
        ab = v["antibody"][0]
        pp = v["ppi"][0]
        sa, ka = chain_seq_and_positions(ab[0], ab[1], "antibody")
        sb, kb = chain_seq_and_positions(pp[0], pp[1], "ppi")
        if len(sa) < 40 or len(sb) < 40:
            continue
        m = align_index(sa, sb)
        if not m:
            continue
        ident = sum(1 for i, j in m.items() if sa[i] == sb[j]) / max(1, len(m))
        if ident < 0.7:
            continue                       # not really the same protein as solved

        ab_pos = {m[ka[x]] for x in contacts["antibody"].get(ab, set())
                  if x in ka and ka[x] in m}
        pp_pos = {kb[x] for x in contacts["ppi"].get(pp, set()) if x in kb}
        if len(ab_pos) < 5 or len(pp_pos) < 5:
            continue
        inter = len(ab_pos & pp_pos)
        rows.append({
            "uniprot": u,
            "antibody_pdb": ab[0], "ppi_pdb": pp[0],
            "seq_identity": round(ident, 3),
            "n_epitope": len(ab_pos), "n_ppi_site": len(pp_pos),
            "n_shared": inter,
            "jaccard": round(inter / len(ab_pos | pp_pos), 3),
            "frac_epitope_in_ppi_site": round(inter / len(ab_pos), 3),
            "n_antibody_structures": len(v["antibody"]),
            "n_ppi_structures": len(v["ppi"]),
        })

    rows.sort(key=lambda r: -r["jaccard"])
    jac = np.array([r["jaccard"] for r in rows]) if rows else np.array([])
    fe = np.array([r["frac_epitope_in_ppi_site"] for r in rows]) if rows else np.array([])

    out = {
        "n_matched_proteins": len(matched),
        "n_analysed": len(rows),
        "median_jaccard_overlap": float(np.median(jac)) if len(jac) else None,
        "median_frac_epitope_in_ppi_site": float(np.median(fe)) if len(fe) else None,
        "frac_pairs_with_no_overlap": float(np.mean(jac == 0)) if len(jac) else None,
        "pairs": rows,
        "note": ("Matching is on UniProt accession, not description, so 'Green "
                 "fluorescent protein' is not merged with 'GFP-like chromoprotein'. "
                 "Residue numbering is reconciled by pairwise alignment between the "
                 "two structures' copies of the protein; pairs below 70% identity "
                 "are dropped. One representative structure per arm per protein."),
        "interpretation": (
            "Low overlap means antibodies and natural partners bind DIFFERENT "
            "surfaces of the same protein, so the corpus-wide chemical differences "
            "cannot be attributed to target identity alone -- the same protein is "
            "read in two different places, in two different ways."),
    }
    (C.OUT_MET / "matched_targets.json").write_text(json.dumps(out, indent=1))

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
    ax = axes[0]
    if len(jac):
        ax.hist(jac, bins=np.linspace(0, 1, 21), color="#762a83", alpha=0.8)
        ax.axvline(np.median(jac), color="k", ls="--", lw=1.3,
                   label=f"median {np.median(jac):.2f}")
        ax.legend(fontsize=9, frameon=False)
    ax.set_xlabel("Jaccard overlap: epitope vs natural-partner site")
    ax.set_ylabel("matched proteins")
    ax.set_title(f"Do antibodies bind where the natural partner binds?\n"
                 f"{len(rows)} proteins solved with both", fontsize=10)

    ax = axes[1]
    top = rows[:16]
    y = np.arange(len(top))
    ax.barh(y, [r["jaccard"] for r in top], color="#762a83")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['uniprot']}  {r['antibody_pdb']}/{r['ppi_pdb']}"
                        for r in top], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Jaccard overlap")
    ax.set_title("Highest-overlap matched pairs", fontsize=10)

    fig.suptitle("Matched targets - the same protein read by an antibody and by a "
                 "natural partner", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.savefig(C.OUT_FIG / "case_matched_targets.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"analysed {len(rows)} matched pairs")
    if rows:
        print(f"median Jaccard overlap = {np.median(jac):.3f}; "
              f"{np.mean(jac == 0):.0%} of pairs share NO contacted position")
        print(f"median fraction of epitope inside the natural-partner site = "
              f"{np.median(fe):.3f}")
        print("\ntop overlaps:")
        for r in rows[:8]:
            print(f"   {r['uniprot']:10s} {r['antibody_pdb']}/{r['ppi_pdb']}  "
                  f"J={r['jaccard']:.2f}  epitope={r['n_epitope']:3d} "
                  f"site={r['n_ppi_site']:3d} shared={r['n_shared']:3d}")
