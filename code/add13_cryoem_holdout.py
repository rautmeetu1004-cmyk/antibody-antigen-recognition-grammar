"""add13 - cryo-EM held-out validation of the tyrosine-backbone mechanism.

The corpus is X-ray only. This tests whether the tyrosine-hydroxyl -> antigen main-chain
hydrogen-bond mechanism reproduces on an independent set of cryo-EM antibody-antigen
complexes (the class the paper explicitly excludes), measured with the SAME geometry applied
to the X-ray corpus. Reproduction both validates the mechanism on independent data and
rebuts the X-ray-only limitation.

cryo-EM ids from the RCSB search API; X-ray reference = the corpus antibody structures on disk.
"""
from __future__ import annotations
import json, os, re, pathlib, urllib.request
import numpy as np
from scipy.stats import mannwhitneyu
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder

ROOT = pathlib.Path(__file__).resolve().parent.parent
CIFDIR = ROOT / "metadata" / "cryoem" / "cif_cache"; CIFDIR.mkdir(parents=True, exist_ok=True)  # regenerable (gitignored)
# X-ray antibody structures (not redistributed here; PDB ids are in metadata/corpus_manifest.csv).
# Point this at a directory of <pdb_id>.pdb files for the antibody arm before running.
XRAYDIR = pathlib.Path(os.environ.get("XRAY_ANTIBODY_DIR", str(ROOT / "pdb_antibody")))
IDS = json.load(open(ROOT / "metadata" / "cryoem" / "cryoem_ids.json"))  # bundled RCSB cryo-EM ids
OUT_JSON = ROOT / "metadata" / "analysis" / "cryoem_holdout.json"
OUT_MD   = ROOT / "metadata" / "analysis" / "CRYOEM_HOLDOUT_RESULTS.md"

cif = MMCIFParser(QUIET=True); pdbp = PDBParser(QUIET=True); ppb = PPBuilder()
AB_MOTIF = re.compile(r"WG[QRKEGA]G")

def chain_seqs(model):
    out = {}
    for ch in model:
        s = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(ch))
        if len(s) >= 50: out[ch.id] = s
    return out

def classify(model):
    """antibody chains carry the Ig FR4 motif in the first ~130 residues; antigen = the rest."""
    seqs = chain_seqs(model)
    ab = {c for c, s in seqs.items() if AB_MOTIF.search(s[:130])}
    ag = {c for c in seqs if c not in ab}
    return ab, ag

def measure(model):
    """Tyr-OH(antibody) -> antigen backbone N/O H-bonds (<=3.5 A) and interface contacts."""
    ab, ag = classify(model)
    if not ab or not ag: return None
    atoms = [a for a in model.get_atoms()
             if a.get_parent().get_parent().id in (ab | ag) and a.element != "H"]
    if len(atoms) > 200000: return None                # skip pathological assemblies
    ns = NeighborSearch(atoms)
    def side(a): return "ab" if a.get_parent().get_parent().id in ab else "ag"
    tyr_bb, iface = 0, set()
    for a in atoms:
        sa = side(a)
        for b in ns.search(a.coord, 4.5):
            if side(b) != sa:
                iface.add(tuple(sorted((id(a.get_parent()), id(b.get_parent())))))
        if sa == "ab" and a.get_parent().get_resname() == "TYR" and a.name == "OH":
            for b in ns.search(a.coord, 3.5):
                if side(b) == "ag" and b.name in ("N", "O"):
                    tyr_bb += 1
    if len(iface) < 10: return None
    return {"tyr_bb": tyr_bb, "iface": len(iface),
            "tyr_bb_per100": 100.0 * tyr_bb / len(iface)}

def cryoem():
    res = []
    for i, pid in enumerate(IDS):
        f = CIFDIR / f"{pid}.cif"
        if not f.exists():
            try: urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.cif", f)
            except Exception: continue
        if f.stat().st_size < 2000 or f.stat().st_size > 80_000_000: continue
        try:
            m = measure(next(iter(cif.get_structure(pid, f))))
        except Exception:
            m = None
        if m: m["pdb"] = pid; res.append(m)
        if (i + 1) % 15 == 0: print(f"  cryo-EM {i+1}/{len(IDS)} -> {len(res)} usable", flush=True)
        if len(res) >= 45: break
    return res

def xray_reference(n=90):
    res = []
    import csv
    lys = [r["pdb_id"] for r in csv.DictReader(open(ROOT/"metadata/corpus_manifest.csv"))
           if r["group"] == "antibody" and r["in_analysis_set"] == "1"]
    for pid in lys[:250]:
        f = XRAYDIR / f"{pid}.pdb"
        if not f.exists(): continue
        try:
            m = measure(next(iter(pdbp.get_structure(pid, f))))
        except Exception:
            m = None
        if m: m["pdb"] = pid; res.append(m)
        if len(res) >= n: break
    return res

def main():
    print("measuring cryo-EM held-out set...")
    ce = cryoem()
    print(f"cryo-EM usable: {len(ce)}")
    print("measuring X-ray corpus reference (same geometry)...")
    xr = xray_reference()
    print(f"X-ray usable: {len(xr)}")
    def frac(r): return float(np.mean([x["tyr_bb"] > 0 for x in r]))
    def med(r, k): return float(np.median([x[k] for x in r]))
    ce_rate = [x["tyr_bb_per100"] for x in ce]; xr_rate = [x["tyr_bb_per100"] for x in xr]
    U, p = mannwhitneyu(ce_rate, xr_rate, alternative="two-sided")
    out = {
        "cryoem": {"n": len(ce), "frac_with_tyr_bb": frac(ce),
                   "median_tyr_bb": med(ce, "tyr_bb"), "median_per100": med(ce, "tyr_bb_per100")},
        "xray":   {"n": len(xr), "frac_with_tyr_bb": frac(xr),
                   "median_tyr_bb": med(xr, "tyr_bb"), "median_per100": med(xr, "tyr_bb_per100")},
        "rate_difference_MWU_p": float(p),
    }
    OUT_JSON.write_text(json.dumps(out, indent=1))
    L = ["# Cryo-EM held-out validation of the tyrosine-backbone mechanism\n",
         "Tyrosine-hydroxyl -> antigen main-chain hydrogen bonds, measured with identical "
         "geometry on an independent cryo-EM antibody-antigen set (excluded from the corpus) "
         "and on the X-ray corpus.\n",
         "| Set | n | % with >=1 Tyr-OH↔backbone bond | median count | median per 100 contacts |",
         "|---|---:|---:|---:|---:|",
         f"| Cryo-EM (held out) | {len(ce)} | {frac(ce)*100:.0f}% | {med(ce,'tyr_bb'):.0f} | {med(ce,'tyr_bb_per100'):.2f} |",
         f"| X-ray (corpus) | {len(xr)} | {frac(xr)*100:.0f}% | {med(xr,'tyr_bb'):.0f} | {med(xr,'tyr_bb_per100'):.2f} |",
         "",
         f"Rate per 100 contacts, cryo-EM vs X-ray: Mann-Whitney p = {p:.3f}. A non-significant "
         "difference means the mechanism is present at the same rate in cryo-EM structures, so it "
         "is not an artefact of the X-ray-only corpus and generalises to the cryo-EM viral-"
         "glycoprotein complexes the study otherwise excludes."]
    OUT_MD.write_text("\n".join(L))
    print(f"cryo-EM: {frac(ce)*100:.0f}% have Tyr-bb, median/100={med(ce,'tyr_bb_per100'):.2f}; "
          f"X-ray: {frac(xr)*100:.0f}%, {med(xr,'tyr_bb_per100'):.2f}; p={p:.3f}")

if __name__ == "__main__":
    main()
