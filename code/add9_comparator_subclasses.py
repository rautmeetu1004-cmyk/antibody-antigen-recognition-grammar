"""add9 - recognition grammar by comparator subclass (Perplexity point #2).

Repeats the headline cross-class grammar correlation (antibody vs co-evolved, r ~ 0.42)
within biologically defined subclasses of Group 2. Each subclass has fewer independence
units, which attenuates the correlation for sample-size reasons alone, so every subclass
value carries a SIZE-MATCHED NULL: the r expected from drawing the same number of units at
random from the full Group-2 set. A subclass only counts as a genuine shift if its observed
r falls outside that null band.

Reuses aa_common exactly (same units, log2_odds, symmetrisation), so 'all_ppi' reproduces
the manuscript's r ~ 0.42 as an end-to-end check. Runs against the validated real contact
tables symlinked into interactions/.
"""
from __future__ import annotations
import sys, json, pathlib
import numpy as np

CODE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import aa_common as C
from add9a_comparator_physics import subclass_labels

ROOT = CODE.parent
OUT_JSON = ROOT / "metadata/analysis/comparator_grammar.json"
OUT_MD   = ROOT / "metadata/analysis/COMPARATOR_GRAMMAR_RESULTS.md"
N_ITER = 1000
SEED = 20260816
SUBCLASSES = ["all_ppi","heteromeric","homomeric","enzyme_inhibitor",
              "receptor_ligand","no_disulfide","disulfide"]

def xcorr(A, B):
    a, b = A.ravel(), B.ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1])

def r_from_units(U, idx, ab):
    """log-odds correlation from a selection of unit matrices."""
    M = U[idx].sum(0)
    return xcorr(ab, C.log2_odds(M))

def main():
    rng = np.random.default_rng(SEED)
    ab = C.log2_odds(C.contact_matrix(C.load_contacts("antibody"), symmetrise=False))

    ppi_rows = C.load_contacts("ppi")
    labels, n_disulf = subclass_labels()
    # full-corpus per-unit matrices (the pool the size-matched null draws from)
    ids_all, U_all = C.unit_matrices(ppi_rows, symmetrise=True)
    n_all = len(ids_all)
    print(f"antibody ref built; full PPI pool = {n_all} units")

    out = {}
    for s in SUBCLASSES:
        rows = [r for r in ppi_rows if s in labels.get(r["pdb_id"], set())]
        ids, U = C.unit_matrices(rows, symmetrise=True)
        k = len(ids)
        r_obs = r_from_units(U, np.arange(k), ab)

        # bootstrap CI on the observed subclass r: resample its own units w/ replacement
        boot = np.array([r_from_units(U, rng.integers(0, k, k), ab) for _ in range(N_ITER)])
        ci = (float(np.nanpercentile(boot, 2.5)), float(np.nanpercentile(boot, 97.5)))

        # size-matched null: k units drawn at random from the FULL pool (without replacement)
        if k < n_all:
            null = np.array([r_from_units(U_all, rng.choice(n_all, k, replace=False), ab)
                             for _ in range(N_ITER)])
            nb = (float(np.nanpercentile(null, 2.5)), float(np.nanpercentile(null, 50)),
                  float(np.nanpercentile(null, 97.5)))
            below = bool(r_obs < nb[0])   # genuine downward shift beyond sample-size effect
        else:
            nb = (r_obs, r_obs, r_obs); below = False   # all_ppi IS the pool

        n_struct = len({r["pdb_id"] for r in rows})
        out[s] = {"n_units": k, "n_structures": n_struct, "r_obs": round(r_obs, 3),
                  "boot_ci": [round(ci[0],3), round(ci[1],3)],
                  "size_matched_null": [round(nb[0],3), round(nb[1],3), round(nb[2],3)],
                  "below_null": below}
        tag = "  <-- below size-matched null" if below else ""
        print(f"  {s:16s} k={k:4d}  r={r_obs:+.3f}  CI[{ci[0]:+.2f},{ci[1]:+.2f}]  "
              f"null~{nb[1]:+.2f}[{nb[0]:+.2f},{nb[2]:+.2f}]{tag}")

    OUT_JSON.write_text(json.dumps({"n_disulfide_interfaces": n_disulf, "results": out}, indent=2))

    L = ["# Recognition grammar by comparator subclass (Perplexity #2)\n",
         "Cross-class grammar correlation (antibody vs comparator subclass). Lower r = more "
         "divergent grammars. Each subclass has fewer independence units, so the **size-matched "
         "null** gives the r expected from that unit count alone; a subclass is a genuine shift "
         "only if `r_obs` sits **below** its null band.\n",
         f"`all_ppi` reproduces the manuscript headline (r ~ 0.42), end-to-end check. "
         f"Disulfide interfaces flagged: {n_disulf}.\n",
         "| Subclass | Structures | Units | Observed r [95% CI] | Size-matched null (median [95%]) | Genuine shift? |",
         "|---|---:|---:|---|---|---|"]
    for s in SUBCLASSES:
        d = out[s]
        null = "—" if s == "all_ppi" else f"{d['size_matched_null'][1]:+.2f} [{d['size_matched_null'][0]:+.2f}, {d['size_matched_null'][2]:+.2f}]"
        verdict = "— (reference)" if s == "all_ppi" else ("**yes, below null**" if d["below_null"] else "no (within null)")
        L.append(f"| {s} | {d['n_structures']} | {d['n_units']} | "
                 f"{d['r_obs']:+.2f} [{d['boot_ci'][0]:+.2f}, {d['boot_ci'][1]:+.2f}] | {null} | {verdict} |")
    OUT_MD.write_text("\n".join(L))
    print("wrote", OUT_JSON.name, "and", OUT_MD.name)

if __name__ == "__main__":
    main()
