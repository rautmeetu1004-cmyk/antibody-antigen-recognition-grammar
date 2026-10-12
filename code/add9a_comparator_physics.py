"""add9a - physics by comparator subclass (Perplexity point #2, down-payment).

Does the antibody-vs-general physics signal (Fig. 1) depend on WHICH slice of the
protein-protein set it is compared against? For each interface feature we recompute the
antibody-vs-comparator rank-biserial effect against biologically defined subclasses of
Group 2, on unit-level means (one value per cluster30 independence unit).

Self-check: the 'all_ppi' column must reproduce the manuscript's Fig. 1 effect sizes
(apolar ~ -0.52, aromatic ~ +0.53, H-bond ~ +0.44, cation-pi ~ +0.34).

Runs entirely from committed tables plus the validated raw contact table (real data:
heatmap_validation/interactions/ppi_contacts.csv, used only to flag disulfide interfaces).
"""
from __future__ import annotations
import csv, json, pathlib, collections
import numpy as np
from scipy.stats import mannwhitneyu

ROOT = pathlib.Path(__file__).resolve().parent.parent
INV  = ROOT / "interactions/analysis/interface_inventory.csv"
THERMO = ROOT / "interactions/analysis/thermodynamics.csv"
MAN  = ROOT / "metadata/corpus_manifest.csv"
RAW_PPI = pathlib.Path("/Users/meetrayuraut/Documents/heatmap_validation/interactions/ppi_contacts.csv")
OUT_JSON = ROOT / "metadata/analysis/comparator_physics.json"
OUT_MD   = ROOT / "metadata/analysis/COMPARATOR_PHYSICS_RESULTS.md"

FEATURES = [
    ("f_apolar_apolar", "apolar-apolar fraction"),
    ("f_aromatic",      "aromatic-residue fraction"),
    ("f_hbond",         "hydrogen-bond fraction"),
    ("f_cation_pi",     "cation-pi fraction"),
    ("f_salt",          "salt-bridge fraction"),
    ("f_stack",         "aromatic-stacking fraction"),
    ("f_target_sidechain_read", "side-chain-read fraction"),
    ("n_pairs",         "interface size (residue pairs)"),
    ("mean_min_dist",   "mean closest-approach (A)"),
]
THERMO_FEATURES = [
    ("dCp_per_1000A2",     "area-derived dCp per 1000 A^2"),
    ("frac_apolar_buried", "apolar fraction of buried area"),
]

def norm(s): return "".join(ch for ch in s.lower() if ch.isalnum() or ch==" ").strip()

def subclass_labels():
    """Group-2 pdb_id -> set of subclass tags, from manifest + disulfide flag."""
    man = {r["pdb_id"]: r for r in csv.DictReader(open(MAN))
           if r["group"] != "antibody" and r["in_analysis_set"] == "1"}
    # disulfide interfaces: any cross-interface CYS-CYS contact in the real table
    disulf = set()
    with open(RAW_PPI) as fh:
        for row in csv.DictReader(fh):
            if row["ab_aa3"] == "CYS" and row["ag_aa3"] == "CYS":
                disulf.add(row["pdb_id"])
    labels = collections.defaultdict(set)
    for pid, r in man.items():
        labels[pid].add("all_ppi")
        a, b = norm(r["side_a_description"]), norm(r["side_b_description"])
        labels[pid].add("homomeric" if (a and a == b) else "heteromeric")
        both = a + " | " + b
        if "inhibitor" in both:                labels[pid].add("enzyme_inhibitor")
        if "receptor" in both:                 labels[pid].add("receptor_ligand")
        labels[pid].add("disulfide" if pid in disulf else "no_disulfide")
    return labels, len(disulf)

def unit_means(rows, feat):
    """one value per cluster30 unit = mean of the feature over that unit's interfaces."""
    by = collections.defaultdict(list)
    for r in rows:
        v = r.get(feat, "")
        if v not in ("", None):
            by[r["unit"]].append(float(v))
    return np.array([np.mean(v) for v in by.values()])

def rank_biserial(x, y):
    """rank-biserial for x vs y; positive => x (antibody) tends to exceed y."""
    if len(x) < 3 or len(y) < 3: return np.nan, np.nan
    U, p = mannwhitneyu(x, y, alternative="two-sided")
    return 2.0 * U / (len(x) * len(y)) - 1.0, p

def main():
    inv = list(csv.DictReader(open(INV)))
    ab  = [r for r in inv if r["arm"] == "antibody"]
    ppi = [r for r in inv if r["arm"] == "ppi"]           # already deduped to a_to_b
    # attach thermo features onto inventory rows by (arm,pdb_id)
    thermo = {(r["arm"], r["pdb_id"]): r for r in csv.DictReader(open(THERMO))}
    for r in inv:
        t = thermo.get((r["arm"], r["pdb_id"]))
        for f, _ in THERMO_FEATURES:
            r[f] = t[f] if t else ""

    labels, n_disulf = subclass_labels()
    subclasses = ["all_ppi","heteromeric","homomeric","enzyme_inhibitor",
                  "receptor_ligand","no_disulfide","disulfide"]
    # subclass -> its ppi rows
    sub_rows = {s: [r for r in ppi if s in labels.get(r["pdb_id"], set())] for s in subclasses}

    results = {}
    for feat, label in FEATURES + THERMO_FEATURES:
        ab_vals = unit_means(ab, feat)
        row = {"label": label, "antibody_median": float(np.median(ab_vals)),
               "antibody_n_units": int(len(ab_vals)), "subclasses": {}}
        for s in subclasses:
            pv = unit_means(sub_rows[s], feat)
            rb, p = rank_biserial(ab_vals, pv)
            row["subclasses"][s] = {"n_units": int(len(pv)),
                                    "comparator_median": (float(np.median(pv)) if len(pv) else None),
                                    "rank_biserial": (None if np.isnan(rb) else round(float(rb),3)),
                                    "p": (None if np.isnan(p) else float(p))}
        results[feat] = row

    meta = {"n_structures": {s: len({r["pdb_id"] for r in sub_rows[s]}) for s in subclasses},
            "n_disulfide_interfaces": n_disulf,
            "note": "PPI thermo features cover the 1,000-structure structural subsample only."}
    OUT_JSON.write_text(json.dumps({"meta": meta, "results": results}, indent=2))

    # markdown
    L = ["# Physics by comparator subclass (Perplexity #2, down-payment)\n",
         "Antibody-vs-comparator **rank-biserial effect** on unit-level means "
         "(positive = antibody higher). The `all_ppi` column is the manuscript's Fig. 1 "
         "value and serves as the self-check.\n",
         "Structures per subclass: " +
         ", ".join(f"{s} {meta['n_structures'][s]}" for s in subclasses) +
         f". Disulfide interfaces flagged: {n_disulf}.\n"]
    hdr = ["Feature"] + subclasses
    L.append("| " + " | ".join(hdr) + " |")
    L.append("|" + "|".join(["---"]*len(hdr)) + "|")
    for feat, label in FEATURES + THERMO_FEATURES:
        r = results[feat]
        cells = [label]
        for s in subclasses:
            d = r["subclasses"][s]; rb = d["rank_biserial"]
            cells.append("—" if rb is None else f"{rb:+.2f}")
        L.append("| " + " | ".join(cells) + " |")
    L.append("\n*n units per subclass:* " +
             ", ".join(f"{s}={results['f_hbond']['subclasses'][s]['n_units']}" for s in subclasses))
    OUT_MD.write_text("\n".join(L))
    print("wrote", OUT_JSON.name, "and", OUT_MD.name)
    # console self-check
    print("\nSELF-CHECK (all_ppi, should match manuscript Fig.1):")
    for feat in ["f_apolar_apolar","f_aromatic","f_hbond","f_cation_pi"]:
        print(f"  {feat:20s} {results[feat]['subclasses']['all_ppi']['rank_biserial']:+.3f}")

if __name__ == "__main__":
    main()
