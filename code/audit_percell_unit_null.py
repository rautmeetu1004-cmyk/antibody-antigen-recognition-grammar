"""Independence-respecting per-cell significance for the recognition matrices.

Addresses two review points about the per-cell map (Fig. 2 crosses; the 'universal'
shared cells of phase4):

  M1  aa_common.fisher_bh() tests each cell on RAW, non-declustered contact counts, so
      the 45-lysozyme and 39-spike clusters are treated as independent draws and the
      q-values are anti-conservative. Here significance comes instead from the UNIT
      bootstrap the study already uses everywhere else: a cell's two-sided p is how
      often its log2-odds crosses zero across 2,000 cluster30-unit resamples.

  M2  Benjamini-Hochberg was applied over all 400 cells, but the PPI matrix is
      symmetrised, so ~190 off-diagonal cells are exact duplicates. Here BH runs over
      the UNIQUE hypotheses only: 400 for the asymmetric antibody matrix, 210
      (upper triangle + diagonal) for the symmetric PPI matrix.

This does not touch aa_common; it writes an independent, corrected table that can be
compared cell-by-cell with the Fisher/BH result. Run after the contact tables exist.

Output: metadata/analysis/percell_significance_unitnull.json
"""
from __future__ import annotations
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

N_BOOT = C.N_BOOT          # 2000
SEED = C.SEED


def bh_over(pvals, mask):
    """Benjamini-Hochberg across only the cells flagged True in `mask` (the unique
    hypotheses). Returns a full 20x20 q array; masked-out cells are NaN."""
    q = np.full(pvals.shape, np.nan)
    idx = np.argwhere(mask)
    ps = np.array([pvals[i, j] for i, j in idx])
    m = len(ps)
    order = np.argsort(ps)
    prev = 1.0
    qs = np.empty(m)
    for rank, k in enumerate(order[::-1], start=1):          # largest p first
        prev = min(prev, ps[k] * m / (m - rank + 1))
        qs[k] = prev
    for (i, j), qv in zip(idx, qs):
        q[i, j] = qv
    return q, m


def per_cell(arm, symmetrise, seed=SEED):
    rows = C.load_contacts(arm)
    ids, U = C.unit_matrices(rows, symmetrise=symmetrise)     # (n_units, 20, 20)
    n = len(ids)
    point = C.log2_odds(U.sum(0))

    rng = np.random.default_rng(seed)
    flat = U.reshape(n, 400)
    counts = rng.multinomial(n, np.full(n, 1.0 / n), size=N_BOOT).astype(float)
    draws = (counts @ flat).reshape(N_BOOT, 20, 20)
    with np.errstate(divide="ignore", invalid="ignore"):
        lo2 = np.stack([C.log2_odds(draws[b]) for b in range(N_BOOT)])
    # two-sided p: fraction of resamples whose sign disagrees with the point estimate
    with np.errstate(invalid="ignore"):
        frac_cross = np.nanmean(np.sign(lo2) != np.sign(point)[None], axis=0)
    p = np.clip(2.0 * frac_cross, 1.0 / (N_BOOT + 1), 1.0)
    p[~np.isfinite(point)] = np.nan

    # unique-hypothesis mask (M2): full for asymmetric, upper+diagonal for symmetric
    mask = np.isfinite(point)
    if symmetrise:
        mask &= np.triu(np.ones((20, 20), bool))
    q, m = bh_over(np.where(np.isfinite(p), p, 1.0), mask)
    sig = np.isfinite(q) & (q < 0.05)
    return {
        "n_units": n,
        "n_unique_hypotheses": int(m),
        "n_significant_q<0.05": int(np.nansum(sig)),
        "cells": [
            {"probe": C.AAS[i], "target": C.AAS[j],
             "log2_or": round(float(point[i, j]), 3),
             "p_unit_boot": round(float(p[i, j]), 4),
             "q_bh_unique": round(float(q[i, j]), 4)}
            for i in range(20) for j in range(20)
            if mask[i, j] and np.isfinite(q[i, j]) and q[i, j] < 0.05
        ],
    }


if __name__ == "__main__":
    out = {
        "method": ("per-cell two-sided p from 2000 cluster30-unit bootstrap resamples "
                   "(sign-crossing); BH over unique hypotheses only (400 antibody, 210 PPI)"),
        "antibody": per_cell("antibody", symmetrise=False),
        "ppi": per_cell("ppi", symmetrise=True),
    }
    (C.OUT_MET / "percell_significance_unitnull.json").write_text(json.dumps(out, indent=1))
    for arm in ("antibody", "ppi"):
        a = out[arm]
        print(f"{arm:9s} units={a['n_units']} unique_tests={a['n_unique_hypotheses']} "
              f"significant(q<0.05)={a['n_significant_q<0.05']}")
