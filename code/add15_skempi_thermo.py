"""add15 - SKEMPI 2.0 enthalpy/entropy partition by interaction class (manuscript section 3.2).

Reconstructs metadata/analysis/skempi_thermo.json. The area-derived heat-capacity proxy
predicts antibody association is less hydrophobically / entropically driven than co-evolved
PPI. This tests that partition directly against measured wild-type calorimetry in SKEMPI 2.0,
which tabulates dH_wt and dS_wt for the ITC subset.

Method (one wild-type record per complex, classified by SKEMPI's own Hold_out_type):
  AB/AG -> antibody; Pr/PI -> protease-inhibitor; TCR/pMHC -> TCR-pMHC;
  '' -> general PPI; anything else (e.g. 'AB/AG,Pr/PI') -> other.
  dG   = RT ln(Kd),  Kd = Affinity_wt_parsed, T = Temperature (default 298 K)
  dH   = dH_wt (kcal/mol)
  TdS  = T * dS_wt / 1000 (dS_wt is cal/mol/K)
Antibody is tested against the pooled non-antibody set (dH and TdS) by Mann-Whitney U,
with the rank-biserial effect size reported alongside, matching the rest of the study.

SKEMPI 2.0 is a public download (https://life.bsc.es/pid/skempi2); not redistributed here.
Set SKEMPI_CSV or drop skempi_v2.csv into metadata/skempi/.
"""
from __future__ import annotations
import csv, json, math, os, re, pathlib
import numpy as np
from scipy.stats import mannwhitneyu

ROOT = pathlib.Path(__file__).resolve().parent.parent
SKEMPI = pathlib.Path(os.environ.get("SKEMPI_CSV", str(ROOT / "metadata" / "skempi" / "skempi_v2.csv")))
OUT_JSON = ROOT / "metadata" / "analysis" / "skempi_thermo.json"
R = 1.987e-3  # kcal/mol/K

CLASS = {"AB/AG": "antibody", "Pr/PI": "protease-inhibitor", "TCR/pMHC": "TCR-pMHC", "": "general PPI"}


def num(v):
    if v in ("", "n.d.", None):
        return None
    m = re.search(r"-?\d+\.?\d*(?:[eE]-?\d+)?", str(v))
    return float(m.group()) if m else None


def rank_biserial(a, b):
    """Rank-biserial correlation from Mann-Whitney U (a vs b): 2U/(n_a n_b) - 1."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    U, _ = mannwhitneyu(a, b, alternative="two-sided")
    return float(2.0 * U / (len(a) * len(b)) - 1.0)


def main():
    rows = list(csv.DictReader(open(SKEMPI), delimiter=";"))
    seen = {}
    for r in rows:                                   # one wild-type record per complex
        pid = r["#Pdb"].split("_")[0].upper()
        seen.setdefault(pid, r)

    classes = {"antibody": {"dG": [], "dH": [], "TdS": []},
               "general PPI": {"dG": [], "dH": [], "TdS": []},
               "protease-inhibitor": {"dG": [], "dH": [], "TdS": []},
               "TCR-pMHC": {"dG": [], "dH": [], "TdS": []},
               "other": {"dG": [], "dH": [], "TdS": []}}
    for r in seen.values():
        cls = CLASS.get(r["Hold_out_type"].strip(), "other")
        T = num(r["Temperature"]) or 298.0
        kd = num(r["Affinity_wt_parsed"])
        if kd and kd > 0:
            classes[cls]["dG"].append(round(R * T * math.log(kd), 3))
        dh = num(r["dH_wt (kcal mol^(-1))"])
        if dh is not None:
            classes[cls]["dH"].append(dh)
        ds = num(r["dS_wt (cal mol^(-1) K^(-1))"])
        if ds is not None:
            classes[cls]["TdS"].append(round(T * ds / 1000.0, 3))

    json.dump(classes, open(OUT_JSON, "w"), indent=1)

    # antibody vs pooled non-antibody (dH, TdS)
    def pool(field):
        out = []
        for c, v in classes.items():
            if c != "antibody":
                out += v[field]
        return out
    ab_dh, nab_dh = classes["antibody"]["dH"], pool("dH")
    ab_ts, nab_ts = classes["antibody"]["TdS"], pool("TdS")
    pH = mannwhitneyu(ab_dh, nab_dh, alternative="two-sided")[1]
    pS = mannwhitneyu(ab_ts, nab_ts, alternative="two-sided")[1]
    summary = {
        "antibody_dH_median": float(np.median(ab_dh)), "n_antibody_dH": len(ab_dh),
        "nonantibody_dH_median": float(np.median(nab_dh)), "n_nonantibody_dH": len(nab_dh),
        "dH_MWU_p": round(float(pH), 4), "dH_rank_biserial": round(rank_biserial(ab_dh, nab_dh), 3),
        "antibody_TdS_median": float(np.median(ab_ts)),
        "nonantibody_TdS_median": float(np.median(nab_ts)),
        "TdS_MWU_p": round(float(pS), 4), "TdS_rank_biserial": round(rank_biserial(ab_ts, nab_ts), 3),
        "antibody_dG_median": float(np.median(classes["antibody"]["dG"])),
        "generalPPI_dG_median": float(np.median(classes["general PPI"]["dG"])),
    }
    json.dump(summary, open(ROOT / "metadata" / "analysis" / "skempi_thermo_summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
