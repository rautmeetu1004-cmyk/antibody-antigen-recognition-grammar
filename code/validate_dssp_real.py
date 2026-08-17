"""Validate the in-house assignment against the real DSSP program.

Peer review asked for a check against an established implementation. The pydssp
cross-check answered whether our code implements the algorithm correctly, but it
shares our per-chain constraint and so could not bound the inter-chain gap. PDB-REDO
publishes DSSP annotations computed with the actual program on the complete
structure, which removes that limitation.

  https://pdb-redo.eu/dssp/db/<pdbid>/mmcif

Residues are matched between the two files by CA coordinate rather than by residue
number, because the DSSP mmCIF uses label_asym_id and label_seq_id while our tables
use author numbering, and coordinate identity is exact where a numbering map would
need assumptions.

Two quantities are reported. The first is three-state agreement, which is the
validation the reviewer requested. The second is the direction of disagreement,
which is the more informative number here: our implementation runs per chain, so
strands whose bridge partner lies in another chain should appear as residues that
real DSSP calls strand and we call loop. Measuring that rate turns a stated
limitation into a quantified one.

Cached downloads live in dssp_reference/.
"""
from __future__ import annotations

import json
import pathlib
import random
import sys
import time
import urllib.error
import urllib.request

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aa_common as C
import phase6a_structural as P

CACHE = C.ROOT / "dssp_reference"
CACHE.mkdir(exist_ok=True)
N_PER_ARM = 60
SEED = 20260817

# DSSP one-letter codes to three states
DSSP3 = {"H": "helix", "G": "helix", "I": "helix", "P": "loop",
         "E": "strand", "B": "strand", "T": "loop", "S": "loop"}


def fetch(pdb_id):
    path = CACHE / f"{pdb_id.lower()}.cif"
    if path.exists() and path.stat().st_size > 1000:
        return path
    url = f"https://pdb-redo.eu/dssp/db/{pdb_id.lower()}/mmcif"
    try:
        with urllib.request.urlopen(url, timeout=90) as fh:
            data = fh.read()
        if len(data) < 1000:
            return None
        path.write_bytes(data)
        return path
    except Exception:                                            # noqa: BLE001
        return None


def parse_dssp(path):
    """{(x,y,z) rounded -> ss3} from the _dssp_struct_summary loop."""
    text = path.read_text(errors="ignore")
    key = "_dssp_struct_summary."
    if key not in text:
        return {}
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith(key))
    cols, i = [], start
    while i < len(lines) and lines[i].strip().startswith(key):
        cols.append(lines[i].strip().split(".")[1])
        i += 1
    idx = {c: k for k, c in enumerate(cols)}
    need = ("secondary_structure", "x_ca", "y_ca", "z_ca")
    if not all(n in idx for n in need):
        return {}

    out = {}
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#") or line.startswith("loop_") \
                or line.startswith("_"):
            break
        parts = line.split()
        if len(parts) >= len(cols):
            ss = parts[idx["secondary_structure"]]
            try:
                xyz = tuple(round(float(parts[idx[c]]), 1)
                            for c in ("x_ca", "y_ca", "z_ca"))
            except ValueError:
                i += 1
                continue
            out[xyz] = DSSP3.get(ss, "loop")
        i += 1
    return out


if __name__ == "__main__":
    rng = random.Random(SEED)
    plan = []
    for arm in C.ARMS:
        ids = sorted(P.targets(arm))
        rng.shuffle(ids)
        plan += [(arm, p) for p in ids[:N_PER_ARM]]

    states = ["helix", "strand", "loop"]
    si = {s: i for i, s in enumerate(states)}
    conf = np.zeros((3, 3), dtype=int)
    n_struct = n_res = n_missing = 0
    t0 = time.time()

    for k, (arm, pdb) in enumerate(plan, 1):
        ref = fetch(pdb)
        if ref is None:
            n_missing += 1
            continue
        dssp = parse_dssp(ref)
        if not dssp:
            n_missing += 1
            continue

        d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
        path = next((d / nm for nm in (f"{pdb}.pdb", f"{pdb}.cif")
                     if (d / nm).exists()), None)
        if path is None:
            continue
        try:
            model = P.load_model(path)
        except Exception:                                        # noqa: BLE001
            continue

        matched = 0
        for ch in model:
            res = [r for r in ch if r.get_resname().strip() in C.AA3_TO_1
                   and "CA" in r]
            if len(res) < 12:
                continue
            ours8 = P.ks_secondary_structure(res)
            ours3 = [P.SS3.get(c, "loop") for c in ours8]
            for r, o in zip(res, ours3):
                a = r["CA"]
                if a.is_disordered():
                    a = a.disordered_get()
                key = tuple(round(float(v), 1) for v in a.coord)
                t = dssp.get(key)
                if t is None:
                    continue
                conf[si[o], si[t]] += 1
                matched += 1
        if matched:
            n_struct += 1
            n_res += matched
        if k % 20 == 0:
            print(f"  {k}/{len(plan)}  matched so far {n_res:,}  "
                  f"{time.time()-t0:4.0f}s", flush=True)

    total = conf.sum()
    if total == 0:
        raise SystemExit("no residues matched; check the coordinate key")
    po = np.trace(conf) / total
    pe = sum((conf[i].sum() / total) * (conf[:, i].sum() / total) for i in range(3))
    kappa = (po - pe) / (1 - pe)

    # the directional question: strands we miss because their partner is elsewhere
    strand_missed = int(conf[si["loop"], si["strand"]])
    strand_total = int(conf[:, si["strand"]].sum())
    helix_missed = int(conf[si["loop"], si["helix"]])
    helix_total = int(conf[:, si["helix"]].sum())

    out = {
        "reference": "DSSP via PDB-REDO precomputed mmCIF, whole structure",
        "n_structures_compared": n_struct,
        "n_structures_unavailable": n_missing,
        "n_residues_compared": int(total),
        "three_state_agreement": float(po),
        "cohens_kappa": float(kappa),
        "confusion_ours_rows_dssp_cols": conf.tolist(),
        "states": states,
        "strand_recall": float(1 - strand_missed / strand_total)
        if strand_total else None,
        "strand_called_loop_by_us": strand_missed,
        "strand_total_in_dssp": strand_total,
        "helix_recall": float(1 - helix_missed / helix_total)
        if helix_total else None,
        "interpretation": (
            "Residues that DSSP calls strand and we call loop are the expected "
            "signature of our per-chain implementation, since a bridge partner in "
            "another chain is invisible to it. The strand recall reported here "
            "bounds that gap directly, replacing the qualitative statement in the "
            "manuscript's limitations."),
    }
    (C.OUT_MET / "dssp_real_validation.json").write_text(json.dumps(out, indent=1))

    print(f"\nstructures compared : {n_struct} ({n_missing} unavailable)")
    print(f"residues compared   : {total:,}")
    print(f"three-state agreement: {po:.4f}")
    print(f"Cohen's kappa        : {kappa:.4f}")
    print("\nconfusion (rows ours, cols real DSSP):")
    print("            " + "".join(f"{s:>10s}" for s in states))
    for s, i in si.items():
        print(f"{s:>10s}  " + "".join(f"{conf[i, j]:10d}" for j in range(3)))
    print(f"\nstrand recall vs DSSP : {out['strand_recall']:.4f} "
          f"({strand_missed:,} of {strand_total:,} DSSP strand residues called loop)")
    print(f"helix recall vs DSSP  : {out['helix_recall']:.4f}")
