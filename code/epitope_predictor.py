"""Epitope predictor: can structure predict which EXPOSED antigen residues are epitope?

Tests the paper's central claim directly. Among solvent-exposed antigen residues, predict
whether each is contacted by the antibody (epitope) from features computable on the UNBOUND
antigen only: amino-acid identity, three-state secondary structure, and relative solvent
accessibility. Leaky features (d_sasa, sasa_complex) are excluded because they encode the
interface. Cross-validation is grouped by cluster30 so no antigen family appears in both
train and test. A near-chance AUC confirms 'the epitope is ordinary exposed surface'; a high
AUC would refute it.

Input:  structural_residues.csv (per-residue structural pass; set STRUCT_CSV or pass as argv[1])
        metadata/corpus_manifest.csv (pdb_id -> cluster30)
Output: metadata/analysis/epitope_predictor.json
"""
from __future__ import annotations
import csv, json, os, sys, pathlib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = pathlib.Path(__file__).resolve().parent.parent
STRUCT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else
                      os.environ.get("STRUCT_CSV", str(ROOT / "interactions/analysis/structural_residues.csv")))
OUT = ROOT / "metadata/analysis/epitope_predictor.json"

AAS = list("ARNDCEQGHILKMFPSTWYV")
AA_IDX = {a: i for i, a in enumerate(AAS)}
# Tien et al. 2013 theoretical max ASA (A^2)
MAXASA = {"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,"I":197,
          "L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,"Y":263,"V":174}
SS3 = {"loop":0, "helix":1, "strand":2}
RSA_EXPOSED = 0.25          # "surface" = relative accessibility >= 25% (standard cutoff)


def load_cluster30():
    m = {}
    for r in csv.DictReader(open(ROOT / "metadata/corpus_manifest.csv")):
        m[r["pdb_id"].upper()] = r.get("cluster30") or r["pdb_id"]
    return m


def build():
    clu = load_cluster30()
    X, y, groups, rsa_only = [], [], [], []
    for r in csv.DictReader(open(STRUCT)):
        if r["arm"] != "antibody" or r["role"] != "target":
            continue
        aa = r["aa"]
        if aa not in AA_IDX:
            continue
        try:
            sasa_free = float(r["sasa_free"])
        except (ValueError, KeyError):
            continue
        rsa = sasa_free / MAXASA[aa]
        if rsa < RSA_EXPOSED:                 # restrict to exposed surface residues
            continue
        feat = [0.0] * 20
        feat[AA_IDX[aa]] = 1.0                # amino-acid one-hot
        ss = [0.0, 0.0, 0.0]
        ss[SS3.get(r.get("ss3", "loop"), 0)] = 1.0
        X.append(feat + ss + [min(rsa, 1.5)])  # +ss3 one-hot + RSA
        rsa_only.append([min(rsa, 1.5)])
        y.append(int(r["is_interface"]))
        groups.append(clu.get(r["pdb_id"].upper(), r["pdb_id"]))
    return np.array(X), np.array(y), np.array(groups), np.array(rsa_only)


def grouped_auc(X, y, groups, model_fn, scale=True):
    gkf = GroupKFold(n_splits=5)
    aucs = []
    for tr, te in gkf.split(X, y, groups):
        Xtr, Xte = X[tr], X[te]
        if scale:
            sc = StandardScaler().fit(Xtr); Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
        m = model_fn().fit(Xtr, y[tr])
        p = m.predict_proba(Xte)[:, 1]
        if len(np.unique(y[te])) > 1:
            aucs.append(roc_auc_score(y[te], p))
    return float(np.mean(aucs)), float(np.std(aucs))


if __name__ == "__main__":
    X, y, groups, rsa = build()
    base = y.mean()
    print(f"exposed antigen residues: {len(y):,}  epitope rate: {base:.3f}  "
          f"antigen clusters: {len(set(groups))}")
    lr_full = grouped_auc(X, y, groups, lambda: LogisticRegression(max_iter=1000, C=1.0))
    lr_rsa  = grouped_auc(rsa, y, groups, lambda: LogisticRegression(max_iter=1000))
    rf_full = grouped_auc(X, y, groups, lambda: RandomForestClassifier(
        n_estimators=300, max_depth=8, n_jobs=-1, random_state=20260816), scale=False)
    out = {
        "task": "predict epitope among exposed antigen residues (RSA>=0.25), cluster30-grouped 5-fold CV",
        "n_residues": int(len(y)), "epitope_rate": round(float(base), 4),
        "n_antigen_clusters": len(set(groups)),
        "features_full": "aa one-hot (20) + ss3 one-hot (3) + relative accessibility (1)",
        "auc_logreg_full":   [round(lr_full[0], 3), round(lr_full[1], 3)],
        "auc_logreg_rsa_only": [round(lr_rsa[0], 3), round(lr_rsa[1], 3)],
        "auc_randomforest_full": [round(rf_full[0], 3), round(rf_full[1], 3)],
        "interpretation": ("AUC near 0.5-0.65 indicates exposed epitopes are weakly "
                           "distinguishable from other exposed surface, consistent with the "
                           "study's finding that epitopes are ordinary accessible surface; "
                           "compare full vs accessibility-only to see what composition/SS add."),
    }
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
