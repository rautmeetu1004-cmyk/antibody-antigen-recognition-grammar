"""Cross-validate the in-house Kabsch-Sander assignment against pydssp.

Peer review asked for the secondary-structure implementation to be checked against
an established program rather than only against deposited HELIX and SHEET records,
which are author annotations and a coarser reference. No DSSP executable is
available in this environment, and pydssp declares torch as a dependency. However
pydssp's numpy backend (`pydssp_numpy.py`) imports only numpy and einops; the torch
import lives in the package's `main.py`. Loading that module directly by file path
gives an independent implementation of the same algorithm with no torch and no
third-party build.

What this does and does not settle. It is a genuine implementation cross-check:
pydssp was written independently, from the same 1983 definition, so agreement tests
whether our code implements the algorithm faithfully. It does NOT address the
per-chain limitation, because pydssp also assigns from a single ordered backbone
array and so likewise cannot see beta bridges formed between chains. That limitation
requires a full DSSP binary and remains stated in the manuscript.

pydssp returns three states (helix, strand, loop), which matches the three-state
form in which every secondary-structure claim in this study is made.

Run with the isolated interpreter that has einops:
  /tmp/dsspenv/bin/python code/validate_dssp.py
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import random
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

PYDSSP_PATH = "/tmp/pydssp_check/pd/pydssp/pydssp_numpy.py"
N_STRUCTURES = 120
SEED = 20260817


def load_pydssp():
    spec = importlib.util.spec_from_file_location("pydssp_numpy", PYDSSP_PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    import aa_common as C
    import phase6a_structural as P

    pd = load_pydssp()
    rng = random.Random(SEED)

    plan = []
    for arm in C.ARMS:
        ids = sorted(P.targets(arm))
        rng.shuffle(ids)
        plan += [(arm, p) for p in ids[:N_STRUCTURES // 2]]

    # three-state confusion matrix, ours (rows) against pydssp (cols)
    states = ["helix", "strand", "loop"]
    conf = np.zeros((3, 3), dtype=int)
    si = {s: i for i, s in enumerate(states)}
    n_chains = n_res = n_fail = 0

    for arm, pdb in plan:
        d = C.ROOT / ("pdb_antibody" if arm == "antibody" else "pdb_ppi")
        path = next((d / nm for nm in (f"{pdb}.pdb", f"{pdb}.cif")
                     if (d / nm).exists()), None)
        if path is None:
            continue
        try:
            model = P.load_model(path)
        except Exception:                                        # noqa: BLE001
            n_fail += 1
            continue

        for ch in model:
            res = [r for r in ch if r.get_resname().strip() in C.AA3_TO_1]
            # pydssp needs a complete N, CA, C, O for every position
            res = [r for r in res if all(a in r for a in ("N", "CA", "C", "O"))]
            if len(res) < 12:
                continue
            coord = np.array([[r["N"].coord, r["CA"].coord, r["C"].coord,
                               r["O"].coord] for r in res], dtype=float)
            try:
                # the numpy backend takes only (coord, donor_mask) and returns a
                # one-hot array ordered [loop, helix, strand]; the out_type argument
                # belongs to the package wrapper in main.py, which imports torch
                onehot = pd.assign(coord)
            except Exception as e:                               # noqa: BLE001
                n_fail += 1
                if n_fail <= 3:
                    print(f"  pydssp failed on {pdb} chain {ch.id}: "
                          f"{type(e).__name__}: {e}")
                continue
            theirs = ["loop", "helix", "strand"]
            theirs = [theirs[i] for i in np.argmax(onehot, axis=-1)]

            ours8 = P.ks_secondary_structure(res)
            ours3 = [P.SS3.get(c, "loop") for c in ours8]
            their3 = theirs
            if len(ours3) != len(their3):
                n_fail += 1
                continue
            for a, b in zip(ours3, their3):
                conf[si[a], si[b]] += 1
            n_chains += 1
            n_res += len(res)

    total = conf.sum()
    agree = np.trace(conf)
    per_state = {}
    for s, i in si.items():
        row, col = conf[i].sum(), conf[:, i].sum()
        per_state[s] = {
            "ours_n": int(row), "pydssp_n": int(col),
            "recall_vs_pydssp": float(conf[i, i] / col) if col else None,
            "precision_vs_pydssp": float(conf[i, i] / row) if row else None,
        }

    # Cohen's kappa, which corrects for agreement expected by chance
    po = agree / total
    pe = sum((conf[i].sum() / total) * (conf[:, i].sum() / total) for i in range(3))
    kappa = (po - pe) / (1 - pe) if pe < 1 else float("nan")

    out = {
        "reference": "pydssp 0.9.1 numpy backend, loaded directly (no torch)",
        "n_structures_attempted": len(plan),
        "n_chains_compared": n_chains,
        "n_residues_compared": int(total),
        "n_failed": n_fail,
        "three_state_agreement": float(po),
        "cohens_kappa": float(kappa),
        "confusion_ours_rows_pydssp_cols": conf.tolist(),
        "states": states,
        "per_state": per_state,
        "scope_note": (
            "This is an implementation cross-check. pydssp assigns from a single "
            "ordered backbone array, as our implementation does, so it shares the "
            "per-chain limitation and cannot detect beta bridges between chains. "
            "Agreement here tests fidelity to the algorithm, not coverage of "
            "inter-chain sheets."),
    }
    (C.OUT_MET / "dssp_validation.json").write_text(json.dumps(out, indent=1))

    print(f"chains compared      : {n_chains}")
    print(f"residues compared    : {total:,}")
    print(f"three-state agreement: {po:.4f}")
    print(f"Cohen's kappa        : {kappa:.4f}")
    print(f"failures             : {n_fail}")
    print("\nconfusion (rows ours, cols pydssp):")
    print("            " + "".join(f"{s:>10s}" for s in states))
    for s, i in si.items():
        print(f"{s:>10s}  " + "".join(f"{conf[i, j]:10d}" for j in range(3)))
    print()
    for s, v in per_state.items():
        r_ = v['recall_vs_pydssp']; p_ = v['precision_vs_pydssp']
        print(f"{s:>7s}  ours {v['ours_n']:7d}  pydssp {v['pydssp_n']:7d}  "
              f"recall {r_ if r_ is None else round(r_, 3)}  "
              f"precision {p_ if p_ is None else round(p_, 3)}")


if __name__ == "__main__":
    main()
