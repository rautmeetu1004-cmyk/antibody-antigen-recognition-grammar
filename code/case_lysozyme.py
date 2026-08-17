"""Case study - lysozyme C: what does the immune system actually target?

Forty-two independently raised antibodies against one 129-residue antigen. This asks
three questions of that natural experiment which the corpus-wide analyses cannot:

  CONSERVATION. Are immunodominant positions conserved or variable? Corpus-wide
  conservation is NOT computable here -- ConSurf-DB is unreachable and no HMMER or
  BLAST is installed -- so this is deliberately case-study scale: UniRef50 homologues
  are pulled from UniProt and aligned pairwise to the reference with Biopython, and
  per-position identity is the conservation measure. Reported as a case study, not
  generalised to the corpus.

  FUNCTIONAL SITE. Lysozyme's catalytic residues are Glu35 and Asp52 (mature
  numbering). Do antibodies preferentially target the active site -- which would bear
  on neutralisation -- or merely its surroundings?

  IMMUNODOMINANCE vs SURFACE. Is antibody targeting explained by accessibility alone?

Numbering. UniProt P00698 includes an 18-residue signal peptide, so UniProt position
= mature position + 18. This offset is applied explicitly; getting it wrong would
misplace the catalytic residues by 18 positions.

Outputs
  metadata/analysis/case_lysozyme.json
  figures/case_lysozyme.png
"""
from __future__ import annotations

import json
import pathlib
import sys
import urllib.request
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
import immunodominance as IM

ACC = "P00698"
SIGNAL_PEPTIDE = 18          # UniProt position = mature position + 18


def fetch_json(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return json.load(fh)


def functional_sites():
    d = fetch_json(f"https://rest.uniprot.org/uniprotkb/{ACC}.json"
                   "?fields=ft_act_site,ft_binding,ft_site,ft_disulfid")
    out = []
    for f in d.get("features", []):
        s = f["location"]["start"]["value"]
        e = f["location"]["end"]["value"]
        if s is None or e is None:
            continue
        for p in range(int(s), int(e) + 1):
            mature = p - SIGNAL_PEPTIDE
            if 1 <= mature <= len(IM.HEWL):
                out.append({"type": f["type"], "position": mature,
                            "residue": IM.HEWL[mature - 1],
                            "description": f.get("description", "")})
    return out


def conservation(max_homologues=250):
    """Per-position identity across UniRef50 homologues, aligned to HEWL.

    The UniRef record lists member ACCESSIONS but not their sequences, so those are
    batch-fetched from UniProtKB. Homologues are full-length (they carry an
    18-residue signal peptide that mature HEWL lacks), so end gaps are free --
    otherwise the terminal penalty drags the alignment out of frame.
    """
    from Bio import Align
    d = fetch_json(f"https://rest.uniprot.org/uniref/UniRef50_{ACC}.json")
    accs = []
    for m in d.get("members", []):
        a = (m.get("accessions") or [None])[0]
        if a and 80 <= (m.get("sequenceLength") or 0) <= 400:
            accs.append(a)
    accs = list(dict.fromkeys(accs))[:max_homologues]

    seqs = []
    for i in range(0, len(accs), 90):
        chunk = ",".join(accs[i:i + 90])
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
               + chunk + "&format=fasta")
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                txt = fh.read().decode()
        except Exception:                                        # noqa: BLE001
            continue
        cur = []
        for line in txt.splitlines():
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
    rep = (d.get("representativeMember") or {}).get("sequence", {}).get("value")
    if rep:
        seqs.append(rep)
    seqs = [s for s in dict.fromkeys(seqs) if 80 <= len(s) <= 400]

    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    match = np.zeros(len(IM.HEWL))
    total = np.zeros(len(IM.HEWL))
    used = 0
    for s in seqs:
        try:
            aln = aligner.align(IM.HEWL, s)[0]
        except Exception:                                        # noqa: BLE001
            continue
        a, b = aln[0], aln[1]
        qi = 0
        for ca, cb in zip(a, b):
            if ca != "-":
                if qi < len(IM.HEWL):
                    total[qi] += 1
                    if ca == cb:
                        match[qi] += 1
                qi += 1
        used += 1
    cons = np.divide(match, np.maximum(total, 1))
    return cons, used, len(seqs)


if __name__ == "__main__":
    # ---- immunodominance (recomputed here so the case study is self-contained)
    units = C.units("antibody")
    rows = C.dedup_rows("antibody", C.load_contacts("antibody"))
    cluster = max(
        ((c, sum(1 for p in {r["pdb_id"] for r in rows}
                 if units.get(p, {}).get("cluster30") == c))
         for c in {units[p]["cluster30"] for p in {r["pdb_id"] for r in rows}
                   if p in units}), key=lambda kv: kv[1])[0]

    by_pdb = defaultdict(lambda: defaultdict(set))
    for r in rows:
        if units.get(r["pdb_id"], {}).get("cluster30") != cluster:
            continue
        by_pdb[r["pdb_id"]][r["ag_chain"]].add(r["ag_resi"].strip())

    hits = np.zeros(len(IM.HEWL))
    n_used = 0
    for pdb, chains in sorted(by_pdb.items()):
        ok = False
        for ch, resis in chains.items():
            rm = IM.chain_residues(pdb, ch)
            off, score = IM.best_offset(rm)
            if off is None or score < 0.8:
                continue
            for x in resis:
                v = rm.get(x)
                if v:
                    pos = v[0] - off + 1
                    if 1 <= pos <= len(IM.HEWL):
                        hits[pos - 1] += 1
                        ok = True
        n_used += int(ok)
    frac_ab = hits / max(1, n_used)

    print(f"immunodominance from {n_used} antibody complexes", flush=True)
    sites = functional_sites()
    cons, n_aln, n_seq = conservation()
    print(f"conservation from {n_aln}/{n_seq} UniRef50 homologues", flush=True)

    contacted = frac_ab > 0
    act = sorted({s["position"] for s in sites
                  if s["type"] in ("Active site", "Binding site")})

    from scipy.stats import mannwhitneyu, spearmanr
    rho, prho = spearmanr(cons, frac_ab)
    hot = frac_ab >= 0.5
    cold = (frac_ab > 0) & (frac_ab < 0.5)
    try:
        u, pmw = mannwhitneyu(cons[hot], cons[cold], alternative="two-sided")
        eff = 2 * u / (hot.sum() * cold.sum()) - 1
    except ValueError:
        pmw, eff = float("nan"), float("nan")

    out = {
        "antigen": "Lysozyme C (P00698)",
        "n_antibody_complexes": n_used,
        "conservation": {"n_homologues_aligned": n_aln,
                         "mean_identity": float(cons.mean()),
                         "source": "UniRef50 via UniProt REST, pairwise BLOSUM62 "
                                   "global alignment (case-study scale only)"},
        "functional_sites": sites,
        "catalytic_positions_mature": act,
        "catalytic_targeting": [
            {"position": p, "residue": IM.HEWL[p - 1],
             "frac_antibodies_contacting": round(float(frac_ab[p - 1]), 3),
             "conservation": round(float(cons[p - 1]), 3)} for p in act],
        "immunodominant_top": [
            {"position": int(i + 1), "residue": IM.HEWL[i],
             "frac_antibodies": round(float(frac_ab[i]), 3),
             "conservation": round(float(cons[i]), 3)}
            for i in np.argsort(-frac_ab)[:12]],
        "conservation_vs_immunodominance": {
            "spearman_rho": float(rho), "p": float(prho),
            "hot_vs_cold_effect": float(eff), "hot_vs_cold_p": float(pmw),
            "mean_conservation_hot": float(cons[hot].mean()) if hot.any() else None,
            "mean_conservation_cold": float(cons[cold].mean()) if cold.any() else None,
        },
        "note": ("Conservation here is case-study scale: corpus-wide per-residue "
                 "conservation is not computable in this environment (ConSurf-DB "
                 "unreachable, no HMMER/BLAST installed) and is not claimed."),
    }
    (C.OUT_MET / "case_lysozyme.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(2, 1, figsize=(14, 7.4), sharex=True)
    ax = axes[0]
    ax.bar(range(1, len(IM.HEWL) + 1), frac_ab, color="#2166ac", width=0.9)
    for p in act:
        ax.axvline(p, color="#b2182b", lw=1.6, ls="--", alpha=0.9)
        ax.text(p, ax.get_ylim()[1] * 0.94, f"{IM.HEWL[p-1]}{p}", color="#b2182b",
                fontsize=8, ha="center")
    for i in np.argsort(-frac_ab)[:6]:
        ax.text(i + 1, frac_ab[i] + 0.02, f"{IM.HEWL[i]}{i+1}", ha="center",
                fontsize=8)
    ax.set_ylabel(f"fraction of {n_used} antibodies")
    ax.set_title("Lysozyme C: immunodominance (bars) vs catalytic residues "
                 "(red dashed)", fontsize=11)

    ax = axes[1]
    ax.bar(range(1, len(IM.HEWL) + 1), cons, color="#1b7837", width=0.9)
    for p in act:
        ax.axvline(p, color="#b2182b", lw=1.6, ls="--", alpha=0.9)
    ax.set_ylabel(f"conservation\n({n_aln} UniRef50 homologues)")
    ax.set_xlabel("lysozyme residue (mature numbering)")
    ax.set_title(f"Sequence conservation   "
                 f"(Spearman rho vs immunodominance = {rho:+.2f}, p={prho:.3f})",
                 fontsize=11)
    ax.set_xlim(0, len(IM.HEWL) + 1)
    fig.tight_layout()
    fig.savefig(C.OUT_FIG / "case_lysozyme.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\ncatalytic residues (mature): {act}")
    for t in out["catalytic_targeting"]:
        print(f"   {t['residue']}{t['position']}: contacted by "
              f"{t['frac_antibodies_contacting']:.0%} of antibodies, "
              f"conservation {t['conservation']:.2f}")
    print(f"\nconservation vs immunodominance: rho={rho:+.3f} (p={prho:.3g})")
    cvi = out["conservation_vs_immunodominance"]
    print(f"mean conservation: hot positions {cvi['mean_conservation_hot']}, "
          f"cold {cvi['mean_conservation_cold']}, effect "
          f"{cvi['hot_vs_cold_effect']:+.3f} p={cvi['hot_vs_cold_p']:.3g}")
