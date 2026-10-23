"""grammar_score.py - score any antibody-antigen complex against the recognition grammar.

A command-line tool built on the study's released matrices. Given a PDB/mmCIF antibody-antigen
complex it extracts the interface (4.5 A heavy-atom residue-pair contacts), then reports:
  - grammar score: mean antibody-matrix log-odds over the interface contacts (how well the
    interface's residue pairing matches the learned antibody recognition grammar);
  - a transfer score against the general protein-protein matrix, for comparison;
  - the side-chain vs main-chain readout fractions;
  - the number of tyrosine-hydroxyl -> antigen main-chain hydrogen bonds (the antibody
    signature this study localises), with the contacted PDB positions.

Usage:
    python3 code/grammar_score.py complex.pdb
    python3 code/grammar_score.py complex.cif --json

Matrices are read from interactions/analysis/matrix_{antibody,ppi}_cluster30.csv (log2_or).
"""
from __future__ import annotations
import sys, csv, json, re, pathlib, argparse
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder

ROOT = pathlib.Path(__file__).resolve().parent.parent
AAS = list("ARNDCEQGHILKMFPSTWYV"); IDX = {a: i for i, a in enumerate(AAS)}
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
             'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
             'TYR':'Y','VAL':'V','MSE':'M'}
BB = {'N', 'O', 'C', 'CA', 'OXT'}
# Ig J-region (FR4) motifs: heavy chain WGxG(T), light chain FGxG(T) -- catches both so a
# light chain is never mistaken for the antigen.
AB_MOTIF = re.compile(r"WG[QRKEGA]G|FG[A-Z]G[TS]")
ppb = PPBuilder()


def load_matrix(arm):
    M = np.full((20, 20), np.nan)
    for r in csv.DictReader(open(ROOT / f"interactions/analysis/matrix_{arm}_cluster30.csv")):
        if r["probe_aa"] in IDX and r["target_aa"] in IDX:
            M[IDX[r["probe_aa"]], IDX[r["target_aa"]]] = float(r["log2_or"])
    return M


def parse(path):
    path = str(path)
    p = MMCIFParser(QUIET=True) if path.lower().endswith(".cif") else PDBParser(QUIET=True)
    return next(p.get_structure("x", path).get_models())


def chain_seqs(model):
    out = {}
    for ch in model:
        s = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(ch))
        if len(s) >= 25:
            out[ch.id] = s
    return out


def score(path):
    model = parse(path)
    seqs = chain_seqs(model)
    ab = {c for c, s in seqs.items() if AB_MOTIF.search(s[:130])}
    ag = {c for c, s in seqs.items() if c not in ab and len(s) >= 50}
    if not ab or not ag:
        raise SystemExit("Could not identify both an antibody (Ig framework motif) and a "
                         ">=50-residue protein antigen in this structure.")
    ag = {max(ag, key=lambda c: len(seqs[c]))}            # largest non-antibody chain
    atoms = [a for a in model.get_atoms()
             if a.get_parent().get_parent().id in (ab | ag) and a.element != "H"]
    ns = NeighborSearch(atoms)
    side = lambda a: "ab" if a.get_parent().get_parent().id in ab else "ag"

    Mab, Mppi = load_matrix("antibody"), load_matrix("ppi")
    seen = set(); pair_lo = []; n_side = n_main = 0; tyr_bb = []
    for a in atoms:
        if side(a) != "ab":
            continue
        ra = a.get_parent()
        for b in ns.search(a.coord, 4.5):
            if side(b) != "ag":
                continue
            rb = b.get_parent()
            key = (id(ra), id(rb), a.name, b.name)
            if key in seen:
                continue
            seen.add(key)
            aa_ab = THREE2ONE.get(ra.get_resname()); aa_ag = THREE2ONE.get(rb.get_resname())
            if not aa_ab or not aa_ag:
                continue
            if b.name in BB: n_main += 1
            else: n_side += 1
            lo = Mab[IDX[aa_ab], IDX[aa_ag]], Mppi[IDX[aa_ab], IDX[aa_ag]]
            pair_lo.append(lo)
        if ra.get_resname() == "TYR" and a.name == "OH":
            for b in ns.search(a.coord, 3.5):
                rb = b.get_parent()
                # antigen main-chain N/O of a standard amino acid only (exclude water/ligand O)
                if side(b) == "ag" and b.name in ("N", "O") and rb.get_resname() in THREE2ONE:
                    tyr_bb.append(f"{a.get_parent().get_parent().id}:Tyr{ra.id[1]}-OH..{rb.get_parent().id}:"
                                  f"{rb.get_resname()}{rb.id[1]}:{b.name}")
    arr = np.array(pair_lo)
    tot = n_side + n_main
    return {
        "antibody_chains": sorted(ab), "antigen_chain": sorted(ag),
        "n_interface_contacts": len(arr),
        "grammar_score_antibody_matrix": round(float(np.nanmean(arr[:, 0])), 3) if len(arr) else None,
        "transfer_score_ppi_matrix": round(float(np.nanmean(arr[:, 1])), 3) if len(arr) else None,
        "side_chain_fraction": round(n_side / tot, 3) if tot else None,
        "main_chain_fraction": round(n_main / tot, 3) if tot else None,
        "n_tyrosine_hydroxyl_to_backbone_bonds": len(tyr_bb),
        "tyrosine_backbone_bonds": sorted(set(tyr_bb)),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("structure", help="antibody-antigen complex (.pdb or .cif)")
    ap.add_argument("--json", action="store_true", help="emit JSON only")
    a = ap.parse_args()
    out = score(a.structure)
    if a.json:
        print(json.dumps(out, indent=1))
    else:
        print(f"Antibody {out['antibody_chains']}  vs  antigen {out['antigen_chain']}")
        print(f"  interface contacts            : {out['n_interface_contacts']}")
        print(f"  grammar score (antibody matrix): {out['grammar_score_antibody_matrix']}  "
              f"(higher = better matches the antibody grammar)")
        print(f"  transfer score (PPI matrix)    : {out['transfer_score_ppi_matrix']}")
        print(f"  side-chain / main-chain readout: {out['side_chain_fraction']} / {out['main_chain_fraction']}")
        print(f"  Tyr-OH -> backbone H-bonds     : {out['n_tyrosine_hydroxyl_to_backbone_bonds']}")
        for t in out["tyrosine_backbone_bonds"]:
            print(f"      {t}")
