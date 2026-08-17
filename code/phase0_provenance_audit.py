"""Phase 0 - provenance audit of the inherited contact tables.

This study did not compute `antibody_contacts.csv` / `ppi_contacts.csv`; it inherited
them from a working area that also contains `_superseded_contaminated/` and sits in a
project containing `_FABRICATED_removed/`. Nothing is built on those tables until they
are shown to be re-derivable from the raw coordinates.

Four independent checks:

  1. RE-EXTRACTION. For a random sample of structures, recompute every
     residue-residue contact straight from the PDB with an extractor written here,
     from the criteria in the extraction report, without reading the inherited code.
     Report precision / recall / distance concordance. Gate: >=99% F1.
  2. OVERLAP. Any PDB present in both arms? An antibody complex sitting in the PPI
     control would blunt exactly the comparison this study exists to make.
  3. REDUNDANCY. How concentrated is each corpus? The unit weighting exists to absorb
     hyper-solved antigens (lysozyme, spike, HIV Env); quantify what it is absorbing.
  4. LINEAGE. Do the inherited tables share provenance with the quarantined
     directories?

Writes metadata/analysis/provenance_audit.json. Prints a PASS/FAIL verdict.
"""
from __future__ import annotations

import csv
import json
import pathlib
import random
import sys
from collections import Counter, defaultdict

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

CONTACT_CUTOFF = 4.5          # from metadata/antibody_extraction_report.json
N_SAMPLE_AB = 25
N_SAMPLE_PPI = 15
SEED = 20260816

AA3 = set(C.AA3_TO_1)
_pdb = PDBParser(QUIET=True)
_cif = MMCIFParser(QUIET=True)


def _structure_path(arm: str, pdb_id: str) -> pathlib.Path | None:
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    for name in (f"{pdb_id}.pdb", f"{pdb_id}.cif",
                 f"{pdb_id.lower()}.pdb", f"{pdb_id.lower()}.cif"):
        p = d / name
        if p.exists():
            return p
    return None


def _load_model(path: pathlib.Path):
    parser = _cif if path.suffix == ".cif" else _pdb
    return parser.get_structure(path.stem, str(path))[0]


def _residue_key(res):
    """Residue id as the CSV writes it: resseq with any insertion code appended.

    The inherited tables DO carry insertion codes (e.g. '30A'), which matters --
    Chothia-numbered antibody CDR loops are full of them, and an extractor that
    dropped them would silently merge distinct loop residues. It did not.
    """
    het, resseq, icode = res.get_id()
    return f"{resseq}{icode.strip()}"


def reextract(path: pathlib.Path, chains_a: set[str], chains_b: set[str]):
    """Every residue pair across the two chain sets within 4.5 A, heavy atoms only.

    Written against the stated criteria, not against the inherited implementation --
    an audit that reuses the code it audits proves nothing.
    """
    model = _load_model(path)

    def heavy_atoms(chain_ids):
        out = []
        for ch in model:
            if ch.id not in chain_ids:
                continue
            for res in ch:
                if res.get_resname().strip() not in AA3:
                    continue                      # waters, ligands, nucleotides
                for at in res:
                    if at.element == "H":
                        continue
                    if at.is_disordered():
                        at = at.disordered_get()  # highest-occupancy altloc
                    out.append((at, ch.id, res))
        return out

    A = heavy_atoms(chains_a)
    B = heavy_atoms(chains_b)
    if not A or not B:
        return {}, {"empty_side": True}

    ns = NeighborSearch([a[0] for a in B])
    owner = {id(a[0]): (a[1], a[2]) for a in B}

    pairs: dict[tuple, float] = {}
    icode_seen = set()
    for at, ch_a, res_a in A:
        for hit in ns.search(at.coord, CONTACT_CUTOFF):
            ch_b, res_b = owner[id(hit)]
            ra = _residue_key(res_a)
            rb = _residue_key(res_b)
            if not ra[-1].isdigit():
                icode_seen.add((ch_a, ra))
            if not rb[-1].isdigit():
                icode_seen.add((ch_b, rb))
            key = (ch_a, ra, ch_b, rb)
            d = float(np.linalg.norm(at.coord - hit.coord))
            if key not in pairs or d < pairs[key]:
                pairs[key] = d
    return pairs, {"n_insertion_coded_residues": len(icode_seen)}


def audit_arm(arm: str, n_sample: int, rng: random.Random) -> dict:
    # group the inherited table by structure, without the unit-weighting filters --
    # this is an audit of the FILE, not of the analysis view of it
    by_pdb: dict[str, dict] = defaultdict(dict)
    chains: dict[str, tuple[set, set]] = defaultdict(lambda: (set(), set()))
    with (C.INT / f"{arm}_contacts.csv").open() as fh:
        for row in csv.DictReader(fh):
            pdb = row["pdb_id"]
            # PPI rows are emitted in both directions; audit one canonical direction
            if arm == "ppi" and row.get("direction") != "a_to_b":
                continue
            key = (row["ab_chain"], row["ab_resi"].strip(),
                   row["ag_chain"], row["ag_resi"].strip())
            by_pdb[pdb][key] = float(row["min_dist"])
            chains[pdb][0].add(row["ab_chain"])
            chains[pdb][1].add(row["ag_chain"])

    all_pdbs = sorted(by_pdb)
    sample = rng.sample(all_pdbs, min(n_sample, len(all_pdbs)))

    per_struct, missing = [], []
    tot_tp = tot_fp = tot_fn = 0
    dist_dev = []
    icode_total = 0

    for pdb in sample:
        path = _structure_path(arm, pdb)
        if path is None:
            missing.append(pdb)
            continue
        ca, cb = chains[pdb]
        try:
            mine, info = reextract(path, ca, cb)
        except Exception as e:                     # noqa: BLE001 - report, don't crash
            per_struct.append({"pdb": pdb, "error": f"{type(e).__name__}: {e}"})
            continue
        icode_total += info.get("n_insertion_coded_residues", 0)

        theirs = by_pdb[pdb]
        tp = set(mine) & set(theirs)
        fp = set(mine) - set(theirs)               # I find it, the CSV does not
        fn = set(theirs) - set(mine)               # the CSV claims it, I cannot
        tot_tp += len(tp); tot_fp += len(fp); tot_fn += len(fn)
        for k in tp:
            dist_dev.append(abs(mine[k] - theirs[k]))

        per_struct.append({
            "pdb": pdb, "n_mine": len(mine), "n_theirs": len(theirs),
            "tp": len(tp), "fp": len(fp), "fn": len(fn),
            "recall": round(len(tp) / max(1, len(theirs)), 4),
        })

    prec = tot_tp / max(1, tot_tp + tot_fp)
    rec = tot_tp / max(1, tot_tp + tot_fn)
    f1 = 2 * prec * rec / max(1e-12, prec + rec)
    dev = np.array(dist_dev) if dist_dev else np.array([np.nan])

    return {
        "arm": arm,
        "n_structures_in_file": len(all_pdbs),
        "n_sampled": len(sample),
        "n_structure_files_missing": len(missing),
        "missing": missing[:10],
        "precision": round(prec, 5),
        "recall": round(rec, 5),
        "f1": round(f1, 5),
        "n_pairs_compared": tot_tp + tot_fp + tot_fn,
        "min_dist_max_abs_dev": float(np.nanmax(dev)),
        "min_dist_median_abs_dev": float(np.nanmedian(dev)),
        "n_insertion_coded_residues_in_sample": icode_total,
        "per_structure": per_struct,
    }


def overlap_and_redundancy() -> dict:
    def pdbs(arm):
        with (C.INT / f"{arm}_contacts.csv").open() as fh:
            return {r["pdb_id"] for r in csv.DictReader(fh)}

    ab, ppi = pdbs("antibody"), pdbs("ppi")
    shared = sorted(ab & ppi)

    red = {}
    for arm in C.ARMS:
        u = C.units(arm)
        if arm == "antibody":
            counts = Counter(v.get("cluster30") for v in u.values())
            names = {v.get("cluster30"): v.get("antigen_description") for v in u.values()}
        else:
            counts = Counter(d["cluster30"] for v in u.values()
                             for d in (v.get("a_to_b"), v.get("b_to_a")) if d)
            names = {}
        top = counts.most_common(8)
        red[arm] = {
            "n_units": len(counts),
            "n_structures": len(u),
            "largest_cluster_frac": round(top[0][1] / max(1, sum(counts.values())), 4),
            "top_clusters": [{"cluster30": c, "n": n, "antigen": names.get(c, "")}
                             for c, n in top],
            "n_singleton_units": sum(1 for v in counts.values() if v == 1),
        }
    return {"n_shared_pdbs": len(shared), "shared_pdbs": shared[:25],
            "redundancy": red}


def lineage() -> dict:
    """Do the inherited tables come from the quarantined directories?"""
    hv = pathlib.Path.home() / "Documents" / "heatmap_validation"
    out = {}
    for q in ("_superseded_contaminated", "_superseded", "_archive"):
        d = hv / q
        if not d.exists():
            out[q] = "absent"
            continue
        hits = [str(p.relative_to(d)) for p in d.rglob("*")
                if p.is_file() and ("antibody" in p.name or "ppi" in p.name)]
        out[q] = {"exists": True, "n_matching_files": len(hits), "examples": hits[:10]}
    real = {}
    for arm in C.ARMS:
        p = (C.INT / f"{arm}_contacts.csv").resolve()
        real[arm] = {"resolved_path": str(p),
                     "inside_quarantine": any(q in p.parts for q in
                                              ("_superseded_contaminated", "_superseded",
                                               "_archive", "_FABRICATED_removed"))}
    return {"quarantine_dirs": out, "table_paths": real}


if __name__ == "__main__":
    rng = random.Random(SEED)
    report = {
        "phase": 0,
        "criteria": {"contact_cutoff_A": CONTACT_CUTOFF,
                     "gate_f1": 0.99,
                     "note": "extractor written from the stated criteria, "
                             "independent of the inherited implementation"},
        "reextraction": {
            "antibody": audit_arm("antibody", N_SAMPLE_AB, rng),
            "ppi": audit_arm("ppi", N_SAMPLE_PPI, rng),
        },
    }
    report.update(overlap_and_redundancy())
    report["lineage"] = lineage()

    gates = {
        "antibody_f1": report["reextraction"]["antibody"]["f1"] >= 0.99,
        "ppi_f1": report["reextraction"]["ppi"]["f1"] >= 0.99,
        "no_quarantine_lineage": not any(
            v["inside_quarantine"] for v in report["lineage"]["table_paths"].values()),
    }
    report["gates"] = gates
    report["verdict"] = "PASS" if all(gates.values()) else "FAIL"

    (C.OUT_MET / "provenance_audit.json").write_text(json.dumps(report, indent=1))

    for arm in C.ARMS:
        r = report["reextraction"][arm]
        print(f"{arm:9s} n={r['n_sampled']:3d}  precision={r['precision']:.4f}  "
              f"recall={r['recall']:.4f}  F1={r['f1']:.4f}  "
              f"pairs={r['n_pairs_compared']:,}  "
              f"max|d dev|={r['min_dist_max_abs_dev']:.3f}A  "
              f"missing_files={r['n_structure_files_missing']}")
    print(f"shared PDBs between arms: {report['n_shared_pdbs']}")
    for arm in C.ARMS:
        rd = report["redundancy"][arm]
        print(f"{arm:9s} units={rd['n_units']:5d}  largest cluster "
              f"{rd['largest_cluster_frac']:.1%}  singletons={rd['n_singleton_units']}")
    print("GATES:", gates)
    print("VERDICT:", report["verdict"])
