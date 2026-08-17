"""TOST equivalence test for Fab versus single-domain antibody grammar.

Section 3.6 argued that the two antibody formats share one grammar because their
cross-format correlation (0.18) falls inside the Fab split-half reliability ceiling
(0.25, interval 0.12 to 0.37). That is an overlap argument, and the rest of the paper
holds claims of no difference to the stricter standard of an equivalence test. This
applies that standard: the cross-format correlation and the within-Fab reliability
correlation are each bootstrapped over units, and the difference between them is
tested against an equivalence margin.
"""
import csv, json, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, 'code')
import aa_common as C

rng = np.random.default_rng(C.SEED)
ct = {r['pdb_id']: r['antibody_type'] for r in
      csv.DictReader(open('metadata/corpus_manifest.csv')) if r['group'] == 'antibody'}
ab = C.dedup_rows('antibody', C.load_contacts('antibody'))
fab = [r for r in ab if ct.get(r['pdb_id']) == 'FAB']
vhh = [r for r in ab if ct.get(r['pdb_id']) == 'SD-H']

def M(rows):
    return C.log2_odds(C.contact_matrix(rows, symmetrise=False))

byu_f, byu_v = defaultdict(list), defaultdict(list)
for r in fab: byu_f[r['unit']].append(r)
for r in vhh: byu_v[r['unit']].append(r)
uf, uv = list(byu_f), list(byu_v)
nv = len(uv)

obs_cross = C.matrix_similarity(M(fab), M(vhh))

# bootstrap the cross-format correlation, resampling units in both formats
cross = []
for _ in range(400):
    a = [r for u in rng.choice(uf, size=len(uf), replace=True) for r in byu_f[u]]
    b = [r for u in rng.choice(uv, size=len(uv), replace=True) for r in byu_v[u]]
    v = C.matrix_similarity(M(a), M(b))
    if np.isfinite(v): cross.append(v)

# within-Fab reliability at the same sample size: two disjoint halves of nv units
split = []
for _ in range(400):
    perm = rng.permutation(uf)
    h1 = [r for u in perm[:nv] for r in byu_f[u]]
    h2 = [r for u in perm[nv:2*nv] for r in byu_f[u]]
    if len(h1) < 500 or len(h2) < 500: continue
    v = C.matrix_similarity(M(h1), M(h2))
    if np.isfinite(v): split.append(v)

cross, split = np.array(cross), np.array(split)

# margin: 20 percent of the antibody-versus-PPI grammar gap, the difference this
# study treats as substantive
hm = json.load(open('metadata/analysis/heatmap.json'))
gap = abs(hm['chemical_classes']['class_similarity_antibody_vs_ppi']
          - hm['cluster30']['grammar_similarity_antibody_vs_ppi'])
margin = 0.20 * gap

from scipy.stats import t as tdist
na, nb = len(cross), len(split)
ma, mb = cross.mean(), split.mean()
va, vb = cross.var(ddof=1)/na, split.var(ddof=1)/nb
se = np.sqrt(va+vb)
df = (va+vb)**2 / (va**2/(na-1) + vb**2/(nb-1))
p_low = 1 - tdist.cdf((ma-mb+margin)/se, df)
p_up = tdist.cdf((ma-mb-margin)/se, df)
p_tost = max(p_low, p_up)

diffs = cross[:min(na,nb)] - split[:min(na,nb)]
out = {
 'observed_cross_format_r': float(obs_cross),
 'cross_format_bootstrap_mean': float(ma),
 'within_fab_reliability_mean': float(mb),
 'difference': float(ma-mb),
 'difference_ci': [float(np.percentile(diffs,2.5)), float(np.percentile(diffs,97.5))],
 'equivalence_margin': float(margin),
 'margin_basis': ('20 percent of the class-level minus residue-level grammar gap '
                  f'({gap:.3f}), the difference this study treats as substantive'),
 'tost_p': float(p_tost),
 'equivalent': bool(p_tost < 0.05),
 'n_units_fab': len(uf), 'n_units_nanobody': nv,
}
json.dump(out, open('metadata/analysis/format_equivalence.json','w'), indent=1)
print(f"observed cross-format r        = {obs_cross:.3f}")
print(f"cross-format bootstrap mean    = {ma:.3f}")
print(f"within-Fab reliability mean    = {mb:.3f}")
print(f"difference                     = {ma-mb:+.3f}  CI {out['difference_ci']}")
print(f"equivalence margin             = {margin:.3f}")
print(f"TOST p                         = {p_tost:.4f}  -> "
      f"{'EQUIVALENT' if p_tost < 0.05 else 'not established'}")
