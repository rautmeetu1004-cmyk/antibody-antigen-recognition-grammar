"""Recompute phi/psi/chi1/chi2 after the dihedral-convention bug, in place.

The Phase 6a pass computed torsions with a dihedral routine that returned
-180 - theta. Secondary structure (Kabsch-Sander: H-bond energies + CA angles),
SASA and burial never touch a torsion, so those columns are correct and are kept.
Only the four torsion columns and the rotamer bin are recomputed, which avoids
redoing the expensive Shrake-Rupley step.

The corrected routine is asserted against Bio.PDB.vectors.calc_dihedral before any
structure is processed; if that assertion fails, nothing is written.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import phase6a_structural as P

worst = P.verify_dihedral()
print(f"dihedral verified against Biopython: max deviation {worst:.2e} deg")

src = C.OUT_INT / "structural_residues.csv"
rows = list(csv.DictReader(src.open()))
print(f"{len(rows):,} residue rows to patch")

# group row indices by structure so each file is parsed once
by_struct: dict[tuple, list[int]] = {}
for i, r in enumerate(rows):
    by_struct.setdefault((r["arm"], r["pdb_id"]), []).append(i)
print(f"{len(by_struct):,} structures")

t0 = time.time()
n_done = n_fail = 0
for k, ((arm, pdb), idxs) in enumerate(by_struct.items(), 1):
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    path = next((d / nm for nm in (f"{pdb}.pdb", f"{pdb}.cif") if (d / nm).exists()),
                None)
    if path is None:
        n_fail += 1
        continue
    try:
        model = P.load_model(path)
    except Exception:                                            # noqa: BLE001
        n_fail += 1
        continue

    tor: dict[tuple, tuple] = {}
    for ch in model:
        res = [r for r in ch if r.get_resname().strip() in C.AA3_TO_1]
        for i, r in enumerate(res):
            phi = psi = chi1 = chi2 = None
            try:
                if i > 0 and "C" in res[i - 1] and all(a in r for a in ("N", "CA", "C")):
                    phi = P.dihedral(res[i - 1]["C"].coord, r["N"].coord,
                                     r["CA"].coord, r["C"].coord)
                if i < len(res) - 1 and "N" in res[i + 1] and all(
                        a in r for a in ("N", "CA", "C")):
                    psi = P.dihedral(r["N"].coord, r["CA"].coord, r["C"].coord,
                                     res[i + 1]["N"].coord)
                rn = r.get_resname().strip()
                g = P.CHI1.get(rn)
                if g and all(a in r for a in ("N", "CA", "CB")) and g in r:
                    chi1 = P.dihedral(r["N"].coord, r["CA"].coord, r["CB"].coord,
                                      r[g].coord)
                t = P.CHI2.get(rn)
                if t and "CA" in r and all(a in r for a in t):
                    chi2 = P.dihedral(r["CA"].coord, r[t[0]].coord, r[t[1]].coord,
                                      r[t[2]].coord)
            except Exception:                                    # noqa: BLE001
                pass
            tor[(ch.id, P.residue_key(r))] = (phi, psi, chi1, chi2)

    for i in idxs:
        r = rows[i]
        v = tor.get((r["chain"], r["resi"]))
        if v is None:
            continue
        phi, psi, chi1, chi2 = v
        r["phi"] = "" if phi is None else round(phi, 1)
        r["psi"] = "" if psi is None else round(psi, 1)
        r["chi1"] = "" if chi1 is None else round(chi1, 1)
        r["chi2"] = "" if chi2 is None else round(chi2, 1)
        r["rotamer"] = P.rotamer_bin(chi1)
    n_done += 1
    if k % 100 == 0:
        el = time.time() - t0
        print(f"{k}/{len(by_struct)}  {el:5.0f}s  "
              f"eta {el/k*(len(by_struct)-k)/60:4.1f} min", flush=True)

with src.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

chi = np.array([float(r["chi1"]) for r in rows if r["chi1"] != ""])
h, e = np.histogram(chi, bins=36, range=(-180, 180))
top = sorted(int((e[i] + e[i + 1]) / 2) for i in np.argsort(h)[-3:])
print(f"\npatched {n_done} structures ({n_fail} failed)")
print(f"chi1 modal bins now at {top} (canonical: -60, +60, 180)")
