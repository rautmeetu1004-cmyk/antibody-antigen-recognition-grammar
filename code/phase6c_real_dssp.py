"""Re-run the loop-confound test using secondary structure from the DSSP program.

Section 3.6 rests on secondary structure, and that assignment was made by an
in-house implementation that operates per chain. Validation against the DSSP program
put three-state agreement at 97.67 percent and showed the one systematic gap: strand
recall of 0.950, because a bridge partner in another chain is invisible to a
per-chain method.

That is small, but Section 3.6 is load-bearing, so the honest check is to redo the
test with authoritative assignments rather than argue that 5 percent cannot matter.
DSSP annotations computed on the complete structure were retrieved for 2,767 of the
2,804 corpus entries, and residues are matched by alpha-carbon coordinate.

The test itself is unchanged: the cross-arm grammar correlation is recomputed within
matched structural context, and compared against a null drawn at the same sample
size, because restricting data attenuates correlation on its own.

Outputs
  metadata/analysis/loop_confound_real_dssp.json
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import defaultdict

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aa_common as C
import phase6a_structural as P
import validate_dssp_real as V


def dssp_labels():
    """{(arm, pdb, chain, resi) -> ss3} from real DSSP, matched by CA coordinate."""
    out = {}
    n_struct = n_missing = 0
    for arm in C.ARMS:
        for pdb in sorted(P.targets(arm)):
            ref = V.CACHE / f"{pdb.lower()}.cif"
            if not ref.exists() or ref.stat().st_size < 1000:
                n_missing += 1
                continue
            table = V.parse_dssp(ref)
            if not table:
                n_missing += 1
                continue
            d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
            path = next((d / nm for nm in (f"{pdb}.pdb", f"{pdb}.cif")
                         if (d / nm).exists()), None)
            if path is None:
                continue
            try:
                model = P.load_model(path)
            except Exception:                                    # noqa: BLE001
                continue
            hit = False
            for ch in model:
                for r in ch:
                    if r.get_resname().strip() not in C.AA3_TO_1 or "CA" not in r:
                        continue
                    a = r["CA"]
                    if a.is_disordered():
                        a = a.disordered_get()
                    key = tuple(round(float(v), 1) for v in a.coord)
                    ss = table.get(key)
                    if ss:
                        out[(arm, pdb, ch.id, P.residue_key(r))] = ss
                        hit = True
            n_struct += int(hit)
    return out, n_struct, n_missing


def annotate(arm, labels, level="cluster30"):
    rows = []
    for r in C.dedup_rows(arm, C.load_contacts(arm, level)):
        a = labels.get((arm, r["pdb_id"], r["ab_chain"], r["ab_resi"].strip()))
        b = labels.get((arm, r["pdb_id"], r["ag_chain"], r["ag_resi"].strip()))
        if a is None or b is None:
            continue
        r = dict(r)
        r["ss_probe"], r["ss_target"] = a, b
        rows.append(r)
    return rows


def grammar_r(ab, pp):
    return C.matrix_similarity(
        C.log2_odds(C.contact_matrix(ab, symmetrise=False)),
        C.log2_odds(C.contact_matrix(pp, symmetrise=True)))


def size_matched_null(ab_all, pp_all, n_a, n_p, n_draw=60, seed=C.SEED):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_draw):
        ia = rng.choice(len(ab_all), size=min(n_a, len(ab_all)), replace=False)
        ip = rng.choice(len(pp_all), size=min(n_p, len(pp_all)), replace=False)
        vals.append(grammar_r([ab_all[i] for i in ia], [pp_all[i] for i in ip]))
    v = np.array(vals)
    return float(np.median(v)), [float(np.percentile(v, 2.5)),
                                 float(np.percentile(v, 97.5))]


if __name__ == "__main__":
    print("mapping real DSSP labels onto the corpus...", flush=True)
    labels, n_struct, n_missing = dssp_labels()
    print(f"  {len(labels):,} residues labelled across {n_struct} structures "
          f"({n_missing} without DSSP)", flush=True)

    ann = {arm: annotate(arm, labels) for arm in C.ARMS}
    for arm in C.ARMS:
        print(f"  {arm}: {len(ann[arm]):,} contacts annotated", flush=True)

    # secondary-structure composition per side
    comp = {}
    for arm in C.ARMS:
        for side, key in (("probe", "ss_probe"), ("target", "ss_target")):
            w = defaultdict(float)
            for r in ann[arm]:
                w[r[key]] += r["weight"]
            tot = sum(w.values()) or 1.0
            comp[f"{arm}_{side}"] = {s: round(w.get(s, 0.0) / tot, 4)
                                     for s in ("helix", "strand", "loop")}

    strata = {
        "all": lambda r: True,
        "loop_probe": lambda r: r["ss_probe"] == "loop",
        "loop_both": lambda r: r["ss_probe"] == "loop" and r["ss_target"] == "loop",
        "target_loop": lambda r: r["ss_target"] == "loop",
        "target_strand": lambda r: r["ss_target"] == "strand",
        "target_helix": lambda r: r["ss_target"] == "helix",
    }
    res = {}
    for name, pred in strata.items():
        a = [r for r in ann["antibody"] if pred(r)]
        p = [r for r in ann["ppi"] if pred(r)]
        if len(a) < 300 or len(p) < 300:
            continue
        r_obs = grammar_r(a, p)
        entry = {"r": r_obs, "n_antibody": len(a), "n_ppi": len(p)}
        if name != "all":
            med, ci = size_matched_null(ann["antibody"], ann["ppi"], len(a), len(p))
            entry["null_median"] = med
            entry["null_ci"] = ci
            entry["verdict"] = ("below null, matching context does not explain the "
                                "difference" if r_obs < ci[0] else
                                "within null, consistent with sample size"
                                if r_obs <= ci[1] else "above null")
        res[name] = entry

    out = {
        "source": "DSSP program, whole structure, via PDB-REDO",
        "n_structures_with_dssp": n_struct,
        "n_structures_without_dssp": n_missing,
        "n_residues_labelled": len(labels),
        "ss_composition": comp,
        "loop_confound": res,
        "comparison_note": (
            "This repeats the Section 3.6 test with authoritative secondary "
            "structure. The in-house assignment gave all-contacts 0.387 and "
            "loop-to-loop 0.185, the latter below its size-matched null. If the "
            "conclusion holds here, the per-chain limitation of that assignment does "
            "not drive it."),
    }
    (C.OUT_MET / "loop_confound_real_dssp.json").write_text(json.dumps(out, indent=1))

    print("\nsecondary structure at the interface (real DSSP, unit-weighted):")
    for k, v in comp.items():
        print(f"  {k:18s} " + "  ".join(f"{s}={v[s]:.3f}"
                                        for s in ("helix", "strand", "loop")))
    print("\nloop confound with real DSSP:")
    for k, v in res.items():
        if "null_median" in v:
            print(f"  {k:14s} r={v['r']:+.3f}  null={v['null_median']:+.3f} "
                  f"[{v['null_ci'][0]:+.3f},{v['null_ci'][1]:+.3f}]  "
                  f"n_ab={v['n_antibody']:,}  -> {v['verdict']}")
        else:
            print(f"  {k:14s} r={v['r']:+.3f}  n_ab={v['n_antibody']:,}")
