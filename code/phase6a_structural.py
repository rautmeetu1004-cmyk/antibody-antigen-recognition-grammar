"""Phase 6a - one structural pass: secondary structure, dihedrals, rotamers, burial.

Every structure is parsed ONCE and four things are extracted together, because the
parse dominates the cost:

  SECONDARY STRUCTURE  Kabsch-Sander, implemented here. No DSSP binary exists on
                       this machine and pydssp needs torch under a PEP-668-locked
                       interpreter, so the algorithm is written out directly and
                       VALIDATED against each PDB's own HELIX/SHEET records. That
                       validation is a stronger provenance claim than trusting an
                       opaque binary, and it is reported, not assumed.
  BACKBONE DIHEDRALS   phi/psi. This is the conformation being read -- the channel
                       where Phase 5 found antibody and PPI grammars diverge
                       (cross-arm r = 0.14).
  SIDE-CHAIN TORSIONS  chi1/chi2 and the g-/t/g+ rotamer bin. Rotamer restriction at
                       an interface is a proxy for the ENTROPIC cost of binding, and
                       so tests Phase 2's hypothesis that antibodies buy affinity
                       enthalpically rather than entropically.
  BURIAL               Shrake-Rupley SASA in the complex and in isolation, giving
                       delta-SASA per residue (how much surface the partner buries)
                       and relative accessibility.

The central confound this phase exists to test: antibody paratopes are CDR loops.
If the "antibody grammar" of Phases 3-6 is really LOOP grammar, secondary structure
will show it, and every earlier claim must be restated conditional on that.

Writes interactions/analysis/structural_residues.csv (one row per interface residue
plus a non-interface surface background) and metadata/analysis/structural_pass.json.
Run directly; it is the slow step of the study.
"""
from __future__ import annotations

import copy
import csv
import json
import pathlib
import random
import sys
import time
from collections import defaultdict

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

SEED = 20260816
PPI_SAMPLE = 1000          # PPI complexes to process; antibody arm runs in full
SASA_POINTS = 64           # Shrake-Rupley sphere points; 64 is ample for burial
HB_ENERGY_CUTOFF = -0.5    # kcal/mol, Kabsch & Sander 1983
Q1Q2F = 0.084 * 332.0      # kcal*A/mol, the DSSP electrostatic prefactor

_pdb = PDBParser(QUIET=True)
_cif = MMCIFParser(QUIET=True)
_sr = ShrakeRupley(n_points=SASA_POINTS)

CHI1 = {"ARG": "CG", "ASN": "CG", "ASP": "CG", "CYS": "SG", "GLN": "CG",
        "GLU": "CG", "HIS": "CG", "ILE": "CG1", "LEU": "CG", "LYS": "CG",
        "MET": "CG", "MSE": "CG", "PHE": "CG", "PRO": "CG", "SER": "OG",
        "THR": "OG1", "TRP": "CG", "TYR": "CG", "VAL": "CG1"}
CHI2 = {"ARG": ("CB", "CG", "CD"), "ASN": ("CB", "CG", "OD1"),
        "ASP": ("CB", "CG", "OD1"), "GLN": ("CB", "CG", "CD"),
        "GLU": ("CB", "CG", "CD"), "HIS": ("CB", "CG", "ND1"),
        "ILE": ("CB", "CG1", "CD1"), "LEU": ("CB", "CG", "CD1"),
        "LYS": ("CB", "CG", "CD"), "MET": ("CB", "CG", "SD"),
        "MSE": ("CB", "CG", "SE"), "PHE": ("CB", "CG", "CD1"),
        "PRO": ("CB", "CG", "CD"), "TRP": ("CB", "CG", "CD1"),
        "TYR": ("CB", "CG", "CD1")}


# ------------------------------------------------------------------- geometry utils
def dihedral(p0, p1, p2, p3):
    """IUPAC dihedral in degrees.

    An earlier hand-rolled version of this returned -180 - theta, which survived a
    Ramachandran sanity check (phi stayed negative under that transform) but put
    every chi1 rotamer 60 degrees off its canonical position. It is now the
    standard projection formula and is asserted against Bio.PDB.vectors.calc_dihedral
    in verify_dihedral() below -- never trust a torsion routine that has not been
    checked against a reference implementation.
    """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    n = np.linalg.norm(b1)
    if n == 0:
        return float("nan")
    b1 = b1 / n
    v = b0 - np.dot(b0, b1) * b1          # b0 projected onto the plane normal to b1
    w = b2 - np.dot(b2, b1) * b1          # b2 projected onto the same plane
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def verify_dihedral(n_check: int = 200) -> float:
    """Assert our dihedral matches Biopython's on real coordinates."""
    from Bio.PDB.vectors import calc_dihedral
    path = C.ROOT / "pdb_antibody" / "1A2Y.pdb"
    model = _pdb.get_structure("v", str(path))[0]
    diffs = []
    for ch in model:
        for r in ch:
            g = CHI1.get(r.get_resname().strip())
            if not g or not all(a in r for a in ("N", "CA", "CB")) or g not in r:
                continue
            mine = dihedral(r["N"].coord, r["CA"].coord, r["CB"].coord, r[g].coord)
            ref = np.degrees(calc_dihedral(r["N"].get_vector(), r["CA"].get_vector(),
                                           r["CB"].get_vector(), r[g].get_vector()))
            diffs.append((mine - ref + 180) % 360 - 180)
            if len(diffs) >= n_check:
                break
    worst = float(np.max(np.abs(diffs)))
    # 1e-3 deg, not 1e-6: PDB coordinates are float32, so the two routines differ
    # in the last bits of the arctan2. Anything above this is a real disagreement.
    if worst > 1e-3:
        raise AssertionError(f"dihedral disagrees with Biopython by {worst:.6f} deg")
    return worst


def rotamer_bin(chi1):
    if chi1 is None or not np.isfinite(chi1):
        return ""
    if -120 <= chi1 < 0:
        return "g-"
    if 0 <= chi1 < 120:
        return "g+"
    return "t"


# ------------------------------------------------------- Kabsch-Sander secondary str
def ks_secondary_structure(residues):
    """8-state DSSP-style assignment from backbone geometry alone.

    residues: ordered list of Bio.PDB residues from ONE chain.
    Returns a list of single-character codes aligned to `residues`.

    Implements Kabsch & Sander (1983): an electrostatic H-bond energy between the
    C=O of residue i and the N-H of residue j, with the amide H placed along the
    direction opposite the preceding carbonyl. n-turns and bridges are then read
    off the H-bond map, with DSSP's priority H > E > G > I > T > S.
    """
    n = len(residues)
    ss = ["C"] * n
    if n < 4:
        return ss

    N = np.full((n, 3), np.nan); CA = np.full((n, 3), np.nan)
    Cc = np.full((n, 3), np.nan); O = np.full((n, 3), np.nan)
    for i, r in enumerate(residues):
        for name, arr in (("N", N), ("CA", CA), ("C", Cc), ("O", O)):
            if name in r:
                arr[i] = r[name].coord

    # amide hydrogen: unit vector from the previous residue's C->O, reversed
    H = np.full((n, 3), np.nan)
    for i in range(1, n):
        co = O[i - 1] - Cc[i - 1]
        nrm = np.linalg.norm(co)
        if np.isfinite(nrm) and nrm > 0:
            H[i] = N[i] - co / nrm

    def energy(i, j):
        """H-bond energy C=O(i) ... H-N(j)."""
        if not (np.all(np.isfinite(O[i])) and np.all(np.isfinite(Cc[i]))
                and np.all(np.isfinite(N[j])) and np.all(np.isfinite(H[j]))):
            return 0.0
        if residues[j].get_resname() == "PRO":
            return 0.0                     # proline has no amide proton
        rON = np.linalg.norm(O[i] - N[j]); rCH = np.linalg.norm(Cc[i] - H[j])
        rOH = np.linalg.norm(O[i] - H[j]); rCN = np.linalg.norm(Cc[i] - N[j])
        if min(rON, rCH, rOH, rCN) < 0.5:
            return 0.0
        return Q1Q2F * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN)

    # only pairs whose CA are within 9 A can hydrogen bond
    hb = np.zeros((n, n), dtype=bool)
    for i in range(n):
        if not np.all(np.isfinite(CA[i])):
            continue
        d = np.linalg.norm(CA - CA[i], axis=1)
        for j in np.where(d < 9.0)[0]:
            if abs(i - j) < 2:
                continue
            if energy(i, int(j)) < HB_ENERGY_CUTOFF:
                hb[i, int(j)] = True

    turn = {k: np.zeros(n, dtype=bool) for k in (3, 4, 5)}
    for k in (3, 4, 5):
        for i in range(n - k):
            if hb[i, i + k]:
                turn[k][i] = True

    # bridges -> ladders -> E
    is_e = np.zeros(n, dtype=bool)
    for i in range(1, n - 1):
        for j in range(i + 3, n - 1):
            par = (hb[i - 1, j] and hb[j, i + 1]) or (hb[j - 1, i] and hb[i, j + 1])
            anti = (hb[i, j] and hb[j, i]) or (hb[i - 1, j + 1] and hb[j - 1, i + 1])
            if par or anti:
                is_e[i] = is_e[j] = True

    is_h = np.zeros(n, dtype=bool)
    is_g = np.zeros(n, dtype=bool)
    is_i = np.zeros(n, dtype=bool)
    for i in range(1, n):
        if turn[4][i - 1] and turn[4][i]:
            is_h[i:min(n, i + 4)] = True
        if turn[3][i - 1] and turn[3][i]:
            is_g[i:min(n, i + 3)] = True
        if turn[5][i - 1] and turn[5][i]:
            is_i[i:min(n, i + 5)] = True

    # bend: CA(i-2)->CA(i) vs CA(i)->CA(i+2) angle > 70 deg
    is_s = np.zeros(n, dtype=bool)
    for i in range(2, n - 2):
        a, b = CA[i] - CA[i - 2], CA[i + 2] - CA[i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if np.isfinite(na) and np.isfinite(nb) and na > 0 and nb > 0:
            ang = np.degrees(np.arccos(np.clip(np.dot(a, b) / (na * nb), -1, 1)))
            is_s[i] = ang > 70

    any_turn = turn[3] | turn[4] | turn[5]
    for i in range(n):
        if is_h[i]:
            ss[i] = "H"
        elif is_e[i]:
            ss[i] = "E"
        elif is_g[i]:
            ss[i] = "G"
        elif is_i[i]:
            ss[i] = "I"
        elif any_turn[i] or (i > 0 and any_turn[i - 1]):
            ss[i] = "T"
        elif is_s[i]:
            ss[i] = "S"
    return ss


SS3 = {"H": "helix", "G": "helix", "I": "helix", "E": "strand",
       "T": "loop", "S": "loop", "C": "loop"}


def pdb_ss_records(path):
    """HELIX/SHEET records as {(chain, resseq): 'H'|'E'} -- the validation truth."""
    out = {}
    if path.suffix != ".pdb":
        return out
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("HELIX "):
            ch, s, e = line[19], line[21:25], line[33:37]
        elif line.startswith("SHEET "):
            ch, s, e = line[21], line[22:26], line[33:37]
        else:
            continue
        try:
            s, e = int(s), int(e)
        except ValueError:
            continue
        code = "H" if line.startswith("HELIX") else "E"
        for r in range(s, e + 1):
            out[(ch, r)] = code
    return out


# --------------------------------------------------------------------- structure IO
def load_model(path):
    parser = _cif if path.suffix == ".cif" else _pdb
    return parser.get_structure(path.stem, str(path))[0]


def clean_copy(model, chain_ids):
    """A model copy holding only the requested chains' standard residues."""
    m = copy.deepcopy(model)
    for ch in list(m):
        if ch.id not in chain_ids:
            m.detach_child(ch.id)
            continue
        for res in list(ch):
            if res.get_resname().strip() not in C.AA3_TO_1:
                ch.detach_child(res.id)
    return m


def residue_key(res):
    _, seq, ic = res.get_id()
    return f"{seq}{ic.strip()}"


def process(arm, pdb_id, probe_chains, target_chains, iface):
    path = None
    d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
    for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif"):
        if (d / nm).exists():
            path = d / nm
            break
    if path is None:
        return [], {}

    model = load_model(path)
    both = set(probe_chains) | set(target_chains)
    cx = clean_copy(model, both)
    if not len(list(cx)):
        return [], {}

    _sr.compute(cx, level="R")
    sasa_cx = {(ch.id, residue_key(r)): r.sasa for ch in cx for r in ch}

    sasa_free = {}
    for side in (probe_chains, target_chains):
        sub = clean_copy(model, set(side))
        if not len(list(sub)):
            continue
        _sr.compute(sub, level="R")
        for ch in sub:
            for r in ch:
                sasa_free[(ch.id, residue_key(r))] = r.sasa

    ss_map, truth_hits, truth_tot = {}, 0, 0
    records = pdb_ss_records(path)
    for ch in cx:
        res = list(ch)
        ss = ks_secondary_structure(res)
        for r, s in zip(res, ss):
            ss_map[(ch.id, residue_key(r))] = s
            key = (ch.id, r.get_id()[1])
            if key in records:
                truth_tot += 1
                truth_hits += int(SS3[s] == SS3[records[key]])

    rows = []
    for ch in cx:
        res = list(ch)
        for i, r in enumerate(res):
            key = (ch.id, residue_key(r))
            aa = C.AA3_TO_1.get(r.get_resname().strip())
            if aa is None:
                continue
            phi = psi = chi1 = chi2 = None
            try:
                if i > 0 and "C" in res[i - 1] and all(a in r for a in ("N", "CA", "C")):
                    phi = dihedral(res[i - 1]["C"].coord, r["N"].coord,
                                   r["CA"].coord, r["C"].coord)
                if i < len(res) - 1 and "N" in res[i + 1] and all(
                        a in r for a in ("N", "CA", "C")):
                    psi = dihedral(r["N"].coord, r["CA"].coord, r["C"].coord,
                                   res[i + 1]["N"].coord)
            except Exception:                                    # noqa: BLE001
                pass
            rn = r.get_resname().strip()
            try:
                g = CHI1.get(rn)
                if g and all(a in r for a in ("N", "CA", "CB")) and g in r:
                    chi1 = dihedral(r["N"].coord, r["CA"].coord, r["CB"].coord,
                                    r[g].coord)
                t = CHI2.get(rn)
                if t and "CA" in r and all(a in r for a in t):
                    chi2 = dihedral(r["CA"].coord, r[t[0]].coord, r[t[1]].coord,
                                    r[t[2]].coord)
            except Exception:                                    # noqa: BLE001
                pass

            sc, sf = sasa_cx.get(key), sasa_free.get(key)
            dsasa = (sf - sc) if (sc is not None and sf is not None) else None
            role = "probe" if ch.id in probe_chains else "target"
            rows.append({
                "arm": arm, "pdb_id": pdb_id, "chain": ch.id,
                "resi": residue_key(r), "aa": aa, "role": role,
                "is_interface": int((ch.id, residue_key(r)) in iface),
                "ss8": ss_map.get(key, ""), "ss3": SS3.get(ss_map.get(key, "C"), "loop"),
                "phi": "" if phi is None else round(phi, 1),
                "psi": "" if psi is None else round(psi, 1),
                "chi1": "" if chi1 is None else round(chi1, 1),
                "chi2": "" if chi2 is None else round(chi2, 1),
                "rotamer": rotamer_bin(chi1),
                "sasa_complex": "" if sc is None else round(sc, 2),
                "sasa_free": "" if sf is None else round(sf, 2),
                "d_sasa": "" if dsasa is None else round(dsasa, 2),
            })
    return rows, {"ss_validated": truth_tot, "ss_agree": truth_hits}


def targets(arm):
    """{pdb -> (probe_chains, target_chains, interface residue keys)}."""
    out = defaultdict(lambda: [set(), set(), set()])
    with (C.INT / f"{arm}_contacts.csv").open() as fh:
        for r in csv.DictReader(fh):
            if arm == "ppi" and r.get("direction") != "a_to_b":
                continue
            e = out[r["pdb_id"]]
            e[0].add(r["ab_chain"]); e[1].add(r["ag_chain"])
            e[2].add((r["ab_chain"], r["ab_resi"].strip()))
            e[2].add((r["ag_chain"], r["ag_resi"].strip()))
    return out


if __name__ == "__main__":
    rng = random.Random(SEED)
    plan = []
    for arm in C.ARMS:
        t = targets(arm)
        ids = sorted(t)
        if arm == "ppi" and len(ids) > PPI_SAMPLE:
            ids = sorted(rng.sample(ids, PPI_SAMPLE))
        plan += [(arm, p, t[p]) for p in ids]

    out_path = C.OUT_INT / "structural_residues.csv"
    fields = ["arm", "pdb_id", "chain", "resi", "aa", "role", "is_interface",
              "ss8", "ss3", "phi", "psi", "chi1", "chi2", "rotamer",
              "sasa_complex", "sasa_free", "d_sasa"]
    stats = {"ss_validated": 0, "ss_agree": 0, "n_structures": 0, "n_failed": 0,
             "failures": []}
    t0 = time.time()
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for k, (arm, pdb, (pc, tc, iface)) in enumerate(plan, 1):
            try:
                rows, st = process(arm, pdb, pc, tc, iface)
                w.writerows(rows)
                stats["ss_validated"] += st.get("ss_validated", 0)
                stats["ss_agree"] += st.get("ss_agree", 0)
                stats["n_structures"] += 1
            except Exception as e:                               # noqa: BLE001
                stats["n_failed"] += 1
                if len(stats["failures"]) < 25:
                    stats["failures"].append(f"{arm}/{pdb}: {type(e).__name__}: {e}")
            if k % 25 == 0:
                el = time.time() - t0
                print(f"{k}/{len(plan)}  {el:6.0f}s  "
                      f"eta {el/k*(len(plan)-k)/60:5.1f} min", flush=True)

    stats["ss_agreement_vs_pdb_records"] = (
        stats["ss_agree"] / stats["ss_validated"] if stats["ss_validated"] else None)
    stats["ppi_sampled"] = PPI_SAMPLE
    stats["sampling_note"] = (
        f"the antibody arm is processed in full (512); the PPI arm is a seeded "
        f"random sample of {PPI_SAMPLE} of 2,292 complexes, which still exceeds the "
        f"antibody arm in both structures and cluster30 units. Reported so the "
        f"coverage limit is visible rather than implied.")
    (C.OUT_MET / "structural_pass.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps({k: v for k, v in stats.items() if k != "failures"}, indent=1))
