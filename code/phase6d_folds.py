"""Phase 6d - fold composition (CATH / SCOP) and the fold-matched control.

The remaining composition confound. Phase 6c showed the antibody signature is not
loop context. This asks the other half: is it FOLD context? If antibody antigens are
mainly beta-rich viral surface proteins while PPI partners are mainly alpha/beta
enzymes, then a difference in recognition grammar could be a difference in what is
being recognised, not in how.

This is the direct analogue of the nucleic-acid study's target-composition control,
which is the analysis that reframed that paper's headline claim.

Fold assignments come from RCSB (data.rcsb.org GraphQL, instance level -- CATH and
SCOP are chain annotations, not entity annotations). Coverage is partial: recently
deposited chains are often unclassified, and that fraction is reported rather than
quietly dropped.

The control: restrict BOTH arms to target chains carrying the same CATH class, and
recompute the cross-arm grammar correlation. As in Phase 6c, any change is compared
against a SIZE-MATCHED null, because restricting the data shrinks it and small
samples attenuate correlations on their own.

Outputs
  metadata/analysis/fold_annotations.json   the cached RCSB lookup
  metadata/analysis/folds.json              composition + control results
  figures/phase6d_folds.png
"""
from __future__ import annotations

import json
import pathlib
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

GRAPHQL = "https://data.rcsb.org/graphql"
BATCH = 80
CACHE = C.OUT_MET / "fold_annotations.json"

CATH_CLASS = {"1": "mainly alpha", "2": "mainly beta", "3": "alpha/beta",
              "4": "few secondary structures", "6": "special"}


def target_instances():
    """{arm -> {pdb -> set(chain)}} for the side being RECOGNISED.

    antibody: the antigen chains. ppi: partner B chains (the a_to_b direction's
    target), so the two arms are compared on the same footing -- what is being read.
    """
    out = {arm: defaultdict(set) for arm in C.ARMS}
    for arm in C.ARMS:
        for r in C.dedup_rows(arm, C.load_contacts(arm)):
            out[arm][r["pdb_id"]].add(r["ag_chain"])
    return out


def fetch(instance_ids: list[str]) -> dict:
    """RCSB instance annotations, batched, cached, resumable."""
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    todo = [i for i in instance_ids if i not in cache]
    print(f"{len(cache):,} cached, {len(todo):,} to fetch")

    for k in range(0, len(todo), BATCH):
        chunk = todo[k:k + BATCH]
        ids = ",".join(f'"{i}"' for i in chunk)
        q = ('{polymer_entity_instances(instance_ids:[' + ids + ']){rcsb_id '
             'rcsb_polymer_instance_annotation{type annotation_id name}}}')
        body = json.dumps({"query": q}).encode()
        req = urllib.request.Request(GRAPHQL, data=body,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as fh:
                data = json.load(fh)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            print(f"  batch {k//BATCH} failed: {type(e).__name__}", flush=True)
            time.sleep(2)
            continue
        got = (data.get("data") or {}).get("polymer_entity_instances") or []
        for e in got:
            if not e:
                continue
            ann = e.get("rcsb_polymer_instance_annotation") or []
            cache[e["rcsb_id"]] = [
                {"type": a["type"], "id": a["annotation_id"], "name": a.get("name")}
                for a in ann
                if a["type"] in ("CATH", "SCOP", "SCOP2", "ECOD")]
        for i in chunk:                       # record misses so we do not refetch
            cache.setdefault(i, [])
        if (k // BATCH) % 10 == 0:
            CACHE.write_text(json.dumps(cache))
            print(f"  {k + len(chunk):,}/{len(todo):,}", flush=True)
    CACHE.write_text(json.dumps(cache))
    return cache


def cath_of(entries):
    for a in entries:
        if a["type"] == "CATH":
            return a["id"]
    return None


def scop_of(entries):
    for a in entries:
        if a["type"] in ("SCOP2", "SCOP"):
            return a.get("name") or a["id"]
    return None


def size_matched_null(ann_a, ann_p, n_a, n_p, n_draw=40, seed=C.SEED):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_draw):
        ia = rng.choice(len(ann_a), size=min(n_a, len(ann_a)), replace=False)
        ip = rng.choice(len(ann_p), size=min(n_p, len(ann_p)), replace=False)
        a = C.log2_odds(C.contact_matrix([ann_a[i] for i in ia], symmetrise=False))
        b = C.log2_odds(C.contact_matrix([ann_p[i] for i in ip], symmetrise=True))
        out.append(C.matrix_similarity(a, b))
    o = np.array(out)
    return float(np.median(o)), [float(np.percentile(o, 2.5)),
                                 float(np.percentile(o, 97.5))]


if __name__ == "__main__":
    tgt = target_instances()
    inst = sorted({f"{p}.{c}" for arm in C.ARMS for p, chs in tgt[arm].items()
                   for c in chs})
    print(f"{len(inst):,} target chain instances to annotate")
    cache = fetch(inst)

    # per-structure fold label: the CATH class of its target chains
    fold = {arm: {} for arm in C.ARMS}
    for arm in C.ARMS:
        for p, chs in tgt[arm].items():
            cs = [cath_of(cache.get(f"{p}.{c}", [])) for c in chs]
            cs = [c for c in cs if c]
            fold[arm][p] = cs[0] if cs else None

    comp = {}
    for arm in C.ARMS:
        got = [v for v in fold[arm].values() if v]
        comp[arm] = {
            "n_structures": len(fold[arm]),
            "n_with_cath": len(got),
            "coverage": round(len(got) / max(1, len(fold[arm])), 3),
            "class_fractions": {CATH_CLASS.get(k, k): round(v / max(1, len(got)), 4)
                                for k, v in Counter(c.split(".")[0]
                                                    for c in got).most_common()},
            "top_topologies": Counter(".".join(c.split(".")[:3])
                                      for c in got).most_common(10),
        }

    # ---- the fold-matched control
    ann = {arm: [r for r in C.dedup_rows(arm, C.load_contacts(arm))] for arm in C.ARMS}
    for arm in C.ARMS:
        for r in ann[arm]:
            r["cath_class"] = (fold[arm].get(r["pdb_id"]) or ".").split(".")[0]

    control = {"all": {}}
    a_all = ann["antibody"]; p_all = ann["ppi"]
    control["all"] = {
        "r": C.matrix_similarity(
            C.log2_odds(C.contact_matrix(a_all, symmetrise=False)),
            C.log2_odds(C.contact_matrix(p_all, symmetrise=True))),
        "n_antibody": len(a_all), "n_ppi": len(p_all)}

    for cls in ("1", "2", "3"):
        a = [r for r in a_all if r["cath_class"] == cls]
        p = [r for r in p_all if r["cath_class"] == cls]
        if len(a) < 500 or len(p) < 500:
            continue
        r_obs = C.matrix_similarity(
            C.log2_odds(C.contact_matrix(a, symmetrise=False)),
            C.log2_odds(C.contact_matrix(p, symmetrise=True)))
        med, ci = size_matched_null(a_all, p_all, len(a), len(p))
        control[f"cath_class_{cls}"] = {
            "label": CATH_CLASS[cls], "r": r_obs,
            "n_antibody": len(a), "n_ppi": len(p),
            "size_matched_null_median": med, "size_matched_null_ci": ci,
            "verdict": ("below null - fold matching does NOT explain the difference"
                        if r_obs < ci[0] else
                        "within null - difference is consistent with sample size"
                        if r_obs <= ci[1] else "above null - fold matching raises r")}

    out = {"composition": comp, "fold_matched_control": control,
           "note": ("CATH/SCOP coverage is partial because recently deposited chains "
                    "are not yet classified; the covered fraction is reported per arm. "
                    "The control compares observed r against a SIZE-MATCHED null "
                    "drawn from the unrestricted pool, because restriction shrinks "
                    "the data and small samples attenuate correlations by themselves.")}
    (C.OUT_MET / "folds.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
    classes = ["mainly alpha", "mainly beta", "alpha/beta",
               "few secondary structures"]
    x = np.arange(len(classes))
    for arm, off, col in (("antibody", -0.2, "#b2182b"), ("ppi", 0.2, "#2166ac")):
        v = [comp[arm]["class_fractions"].get(c, 0.0) for c in classes]
        axes[0].bar(x + off, v, width=0.38, color=col, label=C.ARM_LABEL[arm])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(classes, fontsize=8, rotation=15, ha="right")
    axes[0].set_ylabel("fraction of targets")
    axes[0].set_title("CATH class of the recognised protein", fontsize=10)
    axes[0].legend(fontsize=8, frameon=False)

    keys = [k for k in control if k != "all"]
    if keys:
        ax = axes[1]
        ys = np.arange(len(keys))
        for k, y in zip(keys, ys):
            d = control[k]
            ax.plot(d["size_matched_null_ci"], [y, y], c="0.6", lw=6,
                    solid_capstyle="butt")
            ax.plot(d["r"], y, "o", ms=9, c="#b2182b")
        ax.axvline(control["all"]["r"], ls="--", c="0.4", lw=1,
                   label=f"all contacts r={control['all']['r']:.2f}")
        ax.set_yticks(ys)
        ax.set_yticklabels([control[k]["label"] for k in keys], fontsize=9)
        ax.set_xlabel("antibody vs PPI grammar correlation (r)")
        ax.set_title("Fold-matched control\n"
                     "grey = size-matched null, red = observed", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Phase 6d - fold composition and the fold-matched control",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(C.OUT_FIG / "phase6d_folds.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for arm in C.ARMS:
        c = comp[arm]
        print(f"\n{arm}: CATH coverage {c['coverage']:.2f} "
              f"({c['n_with_cath']}/{c['n_structures']})")
        print("   ", c["class_fractions"])
    print("\nfold-matched control:")
    for k, d in control.items():
        if k == "all":
            print(f"   all              r={d['r']:+.3f} (n_ab={d['n_antibody']:,})")
        else:
            print(f"   {d['label']:22s} r={d['r']:+.3f}  null="
                  f"{d['size_matched_null_median']:+.3f} "
                  f"[{d['size_matched_null_ci'][0]:+.3f},"
                  f"{d['size_matched_null_ci'][1]:+.3f}]  -> {d['verdict']}")
