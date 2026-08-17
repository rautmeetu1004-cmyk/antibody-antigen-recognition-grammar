"""Shared foundation for the antibody-antigen / protein-protein recognition study.

Two arms, both already extracted (see metadata/*_extraction_report.json):

  antibody  curated antibody-antigen complexes, folded protein antigens only
            (>=50-residue antigen floor)                             512 structures
  ppi       general (non-antibody) protein-protein complexes,
            each extracted in BOTH directions and flagged            2,310 complexes

This module is the protein-protein analogue of the nucleic-acid arm's
`strata_common.py`, and deliberately mirrors its statistics so the two arms are
comparable:

  * loading an arm's residue-pair contact table
  * the independence-unit machinery (cluster30 primary, uniprot sensitivity)
  * the 20x20 statistics: unit-weighted log2 odds ratio, bootstrap CI over UNITS,
    Fisher-exact + BH, per-unit sign replication

Key differences from the nucleic-acid arm, which are real and must not be
glossed over:

  * the partner alphabet is 20 amino acids, not 4 bases, so the matrix is 20x20
    (400 cells) rather than 20x4 (80). Power per cell is far lower and the BH
    correction is over 5x more tests.
  * the ANTIBODY matrix is ASYMMETRIC and directional by construction: rows are
    paratope (antibody) residues, columns are epitope (antigen) residues. The two
    sides are not exchangeable -- one is a germline-encoded, somatically matured
    binder, the other is an arbitrary protein surface.
  * the PPI matrix is SYMMETRIC in principle: two co-evolved partners have no
    probe/target asymmetry. The extractor emitted both directions, so we
    symmetrise explicitly rather than pretending a direction is meaningful.
  * salt-bridge counts are NOT comparable to the nucleic-acid arm (basic N to
    carboxylate O here, basic N to phosphate O there). Same 4 A Barlow-Thornton
    ion-pair geometry, different acceptor chemistry.

Nothing here writes figures; the analysis scripts own their outputs.
"""
from __future__ import annotations

import csv
import json
import pathlib
from collections import defaultdict

import numpy as np
from scipy.stats import fisher_exact

ROOT = pathlib.Path(__file__).resolve().parent.parent
INT = ROOT / "interactions"
MET = ROOT / "metadata"
FIG = ROOT / "figures"

OUT_INT = INT / "analysis"
OUT_MET = MET / "analysis"
OUT_FIG = FIG
for _d in (OUT_INT, OUT_MET, OUT_FIG):
    _d.mkdir(parents=True, exist_ok=True)

ARMS = ("antibody", "ppi")
ARM_LABEL = {
    "antibody": "Antibody-antigen",
    "ppi": "Protein-protein (non-antibody)",
}

# amino-acid ordering / colouring shared with the nucleic-acid arm's figures so
# the two studies' heatmaps can be read side by side without re-learning the axis
AAS = list("RKHWYFQNSTPGAMILVCDE")
AA_GROUP = {**{a: "basic" for a in "RKH"}, **{a: "aromatic" for a in "WYF"},
            **{a: "polar" for a in "QNSTC"}, **{a: "special" for a in "GPA"},
            **{a: "hydrophobic" for a in "MILV"}, **{a: "acidic" for a in "DE"}}
GROUP_COLOR = {"basic": "#2166ac", "aromatic": "#762a83", "polar": "#1b7837",
               "special": "#7f7f7f", "hydrophobic": "#b8860b", "acidic": "#b2182b"}
GROUPS = ("basic", "aromatic", "polar", "special", "hydrophobic", "acidic")

# selenomethionine is a crystallographic substitution for MET, not a biological
# modification; the extractor kept it in the PPI arm. Fold it into MET so it does
# not masquerade as a 21st residue. Genuine PTMs never reached these tables --
# see aa_modified.py, which goes back to the raw PDBs for that question.
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

BOND_COLS = ("salt", "hbond", "vdw", "stack")

UNIT_LEVELS = ("cluster30", "uniprot")     # cluster30 primary, uniprot sensitivity
N_BOOT = 2000
SEED = 20260816


# ------------------------------------------------------------------ unit metadata
_UNITS_CACHE: dict[str, dict] = {}


def units(arm: str) -> dict:
    """Raw independence-unit metadata for an arm, keyed by PDB id.

    antibody: {pdb -> {cluster30, uniprot, antigen_length, resolution, ...}}
    ppi:      {pdb -> {a_to_b: {cluster30, uniprot, ...}, b_to_a: {...},
                       resolution}}
    """
    if arm not in _UNITS_CACHE:
        _UNITS_CACHE[arm] = json.loads((MET / f"{arm}_units.json").read_text())
    return _UNITS_CACHE[arm]


def _unit_of(arm: str, row: dict, level: str) -> str | None:
    """The independence unit a single contact row belongs to.

    For PPI this depends on the row's `direction`: a complex contributes two
    half-interfaces, and each half is clustered by the identity of the protein
    playing the 'probe' role in that direction.
    """
    u = units(arm).get(row["pdb_id"])
    if u is None:
        return None
    if arm == "ppi":
        u = u.get(row.get("direction"))
        if u is None:
            return None
    val = u.get(level)
    # a chain with no UniProt mapping cannot join the sensitivity analysis; it is
    # dropped there rather than silently pooled into one giant 'NA' unit
    if not val or val == "NA":
        return None
    return str(val)


# ----------------------------------------------------------------------- loading
_CONTACT_CACHE: dict[str, list] = {}


def load_contacts(arm: str, level: str = "cluster30") -> list[dict]:
    """Residue-pair contacts for an arm, annotated with `unit` and `weight`.

    One row per protein-residue <-> protein-residue contact. Rows are annotated,
    never aggregated, so downstream scripts can slice by bond type, moiety, or
    structure without re-reading the CSV.

    Weighting: every contact carries weight 1/(structures in its unit), so a
    cluster represented by 45 lysozyme complexes counts once, not 45 times. This
    is the same redundancy control the nucleic-acid arm uses, and it matters more
    here -- the antibody corpus is dominated by a handful of hyper-solved
    antigens (lysozyme, HIV Env, SARS-CoV-2 spike).
    """
    key = f"{arm}:{level}"
    if key in _CONTACT_CACHE:
        return _CONTACT_CACHE[key]

    rows = []
    with (INT / f"{arm}_contacts.csv").open() as fh:
        for row in csv.DictReader(fh):
            a1 = AA3_TO_1.get(row["ab_aa3"])
            a2 = AA3_TO_1.get(row["ag_aa3"])
            if a1 is None or a2 is None:
                continue                      # non-standard residue, not in the 20x20
            unit = _unit_of(arm, row, level)
            if unit is None:
                continue
            row["aa_probe"] = a1              # antibody: paratope. ppi: 'probe' side
            row["aa_target"] = a2             # antibody: epitope.  ppi: 'target' side
            row["unit"] = unit
            for c in BOND_COLS:
                row[c] = int(row[c]) if row.get(c) not in (None, "") else 0
            row["min_dist"] = float(row["min_dist"]) if row["min_dist"] else np.nan
            rows.append(row)

    # 1/n_structures_in_unit, computed on the rows that actually survived
    per_unit_structs = defaultdict(set)
    for r in rows:
        per_unit_structs[r["unit"]].add(r["pdb_id"])
    for r in rows:
        r["weight"] = 1.0 / len(per_unit_structs[r["unit"]])

    _CONTACT_CACHE[key] = rows
    return rows


def dedup_rows(arm: str, rows: list[dict]) -> list[dict]:
    """One row per physical contact, for anything that needs a true sample size.

    The PPI extractor emitted every complex in both directions, so each physical
    residue-residue contact appears TWICE (verified exactly: 406,397 a_to_b rows
    and 406,397 mirrored b_to_a rows). That is correct for the unit-weighted
    matrices -- each half-interface is clustered by its own protein and both
    deserve a vote -- but it is wrong for Fisher exact, which would be handed
    double the evidence it actually has and return anti-conservative p-values.

    Effect sizes (log2 odds ratios) are ratios and are unaffected either way.
    """
    if arm != "ppi":
        return rows
    return [r for r in rows if r.get("direction") == "a_to_b"]


def corpus_summary(arm: str, level: str = "cluster30") -> dict:
    rows = load_contacts(arm, level)
    return {
        "arm": arm,
        "label": ARM_LABEL[arm],
        "n_contacts": len(rows),
        "n_structures": len({r["pdb_id"] for r in rows}),
        "n_units": len({r["unit"] for r in rows}),
        "unit_level": level,
    }


# -------------------------------------------------------------------- the matrix
def contact_matrix(rows, weighted: bool = True, symmetrise: bool = False,
                   bond: str | None = None) -> np.ndarray:
    """20x20 counts, rows = probe residue, cols = target residue.

    bond=None counts every contact; bond='hbond' etc. counts only pairs making
    that bond type (a pair can make several, so the bond matrices overlap).
    symmetrise=True adds the transpose and halves -- correct for PPI, where the
    two partners are exchangeable, and wrong for antibody, where they are not.
    """
    idx = {a: i for i, a in enumerate(AAS)}
    M = np.zeros((20, 20))
    for r in rows:
        if bond is not None and not r[bond]:
            continue
        M[idx[r["aa_probe"]], idx[r["aa_target"]]] += r["weight"] if weighted else 1.0
    if symmetrise:
        M = (M + M.T) / 2.0
    return M


def log2_odds(M: np.ndarray) -> np.ndarray:
    """log2(observed / expected), expected = independent marginals.

    This is the recognition-grammar quantity: how much more often a residue pair
    is in contact than the two residues' interface abundances alone predict. It
    is blind to the fact that some residues are simply commoner at interfaces,
    which is exactly what we want when asking about *preference*.
    """
    total = M.sum()
    if total <= 0:
        return np.full_like(M, np.nan)
    exp = np.outer(M.sum(1), M.sum(0)) / total
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.log2(M / exp)
    out[~np.isfinite(out)] = np.nan
    return out


def unit_matrices(rows, symmetrise: bool = False, bond: str | None = None):
    """Per-unit 20x20 weighted count matrices, stacked as (n_units, 20, 20).

    The bootstrap is a sum over resampled units, and a unit's contribution is
    fixed. Precomputing each unit's matrix once turns 2,000 resamples from 2,000
    passes over 800k rows into a single matrix multiply.
    """
    idx = {a: i for i, a in enumerate(AAS)}
    by_unit = defaultdict(lambda: np.zeros((20, 20)))
    for r in rows:
        if bond is not None and not r[bond]:
            continue
        by_unit[r["unit"]][idx[r["aa_probe"]], idx[r["aa_target"]]] += r["weight"]
    ids = list(by_unit)
    U = np.stack([by_unit[u] for u in ids]) if ids else np.zeros((0, 20, 20))
    if symmetrise:
        U = (U + U.transpose(0, 2, 1)) / 2.0
    return ids, U


def bootstrap_ci(rows, symmetrise: bool = False, bond: str | None = None,
                 n_boot: int = N_BOOT, seed: int = SEED):
    """Percentile CI on the log2 odds ratio, resampling UNITS with replacement.

    Resampling units (not contacts) is the whole point: contacts within one
    structure are massively correlated, so a contact-level bootstrap would report
    absurdly tight intervals. Returns (point, lo, hi, frac_same_sign).

    Implemented as a multinomial draw over units followed by one matmul, which is
    mathematically identical to resampling unit ids and re-accumulating.
    """
    rng = np.random.default_rng(seed)
    ids, U = unit_matrices(rows, symmetrise=symmetrise, bond=bond)
    n = len(ids)
    point = log2_odds(U.sum(0)) if n else np.full((20, 20), np.nan)
    if n == 0:
        return point, point, point, np.zeros((20, 20))

    flat = U.reshape(n, 400)
    # each row of `counts` is one bootstrap replicate's multiplicity per unit
    counts = rng.multinomial(n, np.full(n, 1.0 / n), size=n_boot).astype(float)
    draws = (counts @ flat).reshape(n_boot, 20, 20)

    tot = draws.sum((1, 2))[:, None, None]
    exp = (draws.sum(2)[:, :, None] * draws.sum(1)[:, None, :]) / np.where(tot > 0, tot, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        lo2 = np.log2(draws / exp)
    lo2[~np.isfinite(lo2)] = np.nan

    with np.errstate(invalid="ignore"):
        lo = np.nanpercentile(lo2, 2.5, axis=0)
        hi = np.nanpercentile(lo2, 97.5, axis=0)
        same = np.nanmean(np.sign(lo2) == np.sign(point)[None, :, :], axis=0)
    return point, lo, hi, same


def fisher_bh(M: np.ndarray):
    """Per-cell Fisher exact (cell vs rest of table) + Benjamini-Hochberg.

    Run on RAW counts -- Fisher needs integers and a real sample size, so feeding
    it unit-weighted fractions would invent significance. The weighted matrix
    carries the effect size; this carries the evidence. Both are reported.
    """
    N = M.sum()
    p = np.ones((20, 20))
    for i in range(20):
        for j in range(20):
            a = M[i, j]
            table = [[a, M[i].sum() - a],
                     [M[:, j].sum() - a, N - M[i].sum() - M[:, j].sum() + a]]
            p[i, j] = fisher_exact(np.rint(table).astype(int))[1]

    flat = p.ravel()
    order = np.argsort(flat)
    m = flat.size
    q = np.empty(m)
    prev = 1.0
    for rank, k in enumerate(order[::-1], start=1):
        prev = min(prev, flat[k] * m / (m - rank + 1))
        q[k] = prev
    return p, q.reshape(20, 20)


def group_matrix(M: np.ndarray) -> np.ndarray:
    """Collapse the 20x20 to a 6x6 chemical-class matrix.

    The per-residue matrix is under-powered in its rarest cells (C, W); the class
    matrix is where the shared-chemistry question is actually answerable.
    """
    gi = {g: i for i, g in enumerate(GROUPS)}
    G = np.zeros((6, 6))
    for i, a in enumerate(AAS):
        for j, b in enumerate(AAS):
            G[gi[AA_GROUP[a]], gi[AA_GROUP[b]]] += M[i, j]
    return G


def matrix_similarity(A: np.ndarray, B: np.ndarray) -> float:
    """Pearson r between two log2-odds maps over cells finite in both."""
    ok = np.isfinite(A) & np.isfinite(B)
    if ok.sum() < 3:
        return float("nan")
    return float(np.corrcoef(A[ok], B[ok])[0, 1])
