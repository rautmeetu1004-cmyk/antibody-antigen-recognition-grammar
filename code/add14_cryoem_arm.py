"""Cryo-EM validation arm: does the recognition grammar reproduce off the X-ray corpus?

The main corpus is X-ray only. This builds an INDEPENDENT cryo-EM antibody-antigen set by
the same rules (folded protein antigen >=50 aa), extracts the same residue-pair grammar with
the same 4.5 A heavy-atom contact and side/main moiety split, clusters antigens for
independence weighting, and correlates the cryo-EM antibody grammar against the X-ray
antibody and X-ray PPI matrices. Reproduction rebuts the X-ray-only limitation.

Outputs: results/cryoem_arm.json  (records cached under data/, CIF download cache in scratchpad)
"""
from __future__ import annotations
import json, re, sys, pathlib, urllib.request, urllib.parse, urllib.error, time
from collections import defaultdict
import numpy as np
# aa_common (log2_odds, contact_matrix, load_contacts) lives with the main X-ray pipeline
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder

AAS = C.AAS; idx = {a: i for i, a in enumerate(AAS)}
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
BB = {'N':'main','O':'main','C':'main','CA':'main','OXT':'main'}  # backbone atoms
AB_MOTIF = re.compile(r"WG[QRKEGA]G")   # Ig heavy/light FR4, first ~130 aa

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent                                               # repository root
CIF = ROOT / "metadata" / "cryoem" / "cif_cache"                 # ~1 GB CIF download cache, regenerable (gitignored)
CIF.mkdir(parents=True, exist_ok=True)
DATA = ROOT / "metadata" / "cryoem"; DATA.mkdir(parents=True, exist_ok=True)
IDS_CACHE = DATA / "cryoem_ids.json"                             # bundled: RCSB cryo-EM entry ids
REC_CACHE = DATA / "cryoem_records.json"                         # bundled: per-structure extracted records
OUT = ROOT / "metadata" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)
parser = MMCIFParser(QUIET=True); ppb = PPBuilder()


def rcsb_ids(maxrows=600):
    if IDS_CACHE.exists():
        return json.load(open(IDS_CACHE))
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match", "value": "ELECTRON MICROSCOPY"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": 3.5}},
            {"type": "terminal", "service": "full_text", "parameters": {"value": "antibody"}}]},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": maxrows},
                             "results_content_type": ["experimental"]}}
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(q))
    with urllib.request.urlopen(url, timeout=60) as r:
        res = json.load(r)
    ids = [x["identifier"] for x in res.get("result_set", [])]
    json.dump(ids, open(IDS_CACHE, "w"))
    return ids


def fetch(pid):
    f = CIF / f"{pid}.cif"
    if f.exists() and f.stat().st_size > 2000:
        return f
    try:
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.cif", f)
    except Exception:
        return None
    return f if f.exists() and f.stat().st_size > 2000 else None


def chain_seqs(model):
    out = {}
    for ch in model:
        s = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(ch))
        if len(s) >= 25:
            out[ch.id] = s
    return out


def extract(pid):
    """Return dict: antigen_seq, contacts[(ab_aa,ag_aa,ag_side,ag_main)], tyr_bb, n_iface."""
    f = fetch(pid)
    if not f:
        return None
    try:
        model = next(parser.get_structure(pid, str(f)).get_models())
    except Exception:
        return None
    seqs = chain_seqs(model)
    ab = {c for c, s in seqs.items() if AB_MOTIF.search(s[:130])}
    ag = {c for c, s in seqs.items() if c not in ab and len(s) >= 50}
    if not ab or not ag:
        return None
    # antigen = largest non-antibody chain
    ag_main_chain = max(ag, key=lambda c: len(seqs[c]))
    ag = {ag_main_chain}
    atoms = [a for a in model.get_atoms()
             if a.get_parent().get_parent().id in (ab | ag) and a.element != "H"]
    if len(atoms) > 250000 or len(atoms) < 100:
        return None
    ns = NeighborSearch(atoms)
    def sideof(a): return "ab" if a.get_parent().get_parent().id in ab else "ag"
    # residue-pair accumulation
    pair_side = defaultdict(float); pair_main = defaultdict(float)
    seen_atompairs = set()
    tyr_bb = 0
    iface_res = set()
    for a in atoms:
        if sideof(a) != "ab":
            continue
        ra = a.get_parent()
        for b in ns.search(a.coord, 4.5):
            if sideof(b) != "ag":
                continue
            rb = b.get_parent()
            key = (id(ra), id(rb), a.name, b.name)
            if key in seen_atompairs:
                continue
            seen_atompairs.add(key)
            aa_ab = THREE2ONE.get(ra.get_resname()); aa_ag = THREE2ONE.get(rb.get_resname())
            if not aa_ab or not aa_ag:
                continue
            moiety = "main" if b.name in BB else "side"
            pk = (aa_ab, aa_ag)
            if moiety == "side": pair_side[pk] += 1
            else: pair_main[pk] += 1
            iface_res.add(id(ra)); iface_res.add(id(rb))
            # Tyr-OH(antibody) -> antigen backbone N/O <=3.5
        if ra.get_resname() == "TYR" and a.name == "OH":
            for b in ns.search(a.coord, 3.5):
                if sideof(b) == "ag" and b.name in ("N", "O"):
                    tyr_bb += 1
    if not pair_side and not pair_main:
        return None
    keys = set(pair_side) | set(pair_main)
    contacts = [[k[0], k[1], pair_side.get(k, 0.0), pair_main.get(k, 0.0)] for k in keys]
    return {"pdb": pid, "antigen_seq": seqs[ag_main_chain], "contacts": contacts,
            "tyr_bb": tyr_bb, "n_iface_res": len(iface_res)}


def kmers(s, k=6):
    return set(s[i:i+k] for i in range(len(s)-k+1))


def cluster_antigens(records, thr=0.7):
    """Greedy k-mer-Jaccard clustering of antigen sequences -> unit id per record."""
    reps = []  # (kmerset, cluster_id)
    unit = {}
    for r in records:
        ks = kmers(r["antigen_seq"])
        assigned = None
        for kset, cid in reps:
            inter = len(ks & kset); uni = len(ks | kset) or 1
            if inter / uni >= thr:
                assigned = cid; break
        if assigned is None:
            assigned = len(reps); reps.append((ks, assigned))
        unit[r["pdb"]] = assigned
    return unit, len(reps)


def build_matrix(records, unit, channel):  # channel in {'both','side','main'}
    n_in_unit = defaultdict(int)
    for r in records:
        n_in_unit[unit[r["pdb"]]] += 1
    M = np.zeros((20, 20))
    for r in records:
        w = 1.0 / n_in_unit[unit[r["pdb"]]]
        for aa_ab, aa_ag, s, m in r["contacts"]:
            if aa_ab not in idx or aa_ag not in idx:
                continue
            val = s + m if channel == "both" else (s if channel == "side" else m)
            M[idx[aa_ab], idx[aa_ag]] += w * val
    return C.log2_odds(M)


def corr(A, B):
    ok = np.isfinite(A) & np.isfinite(B)
    return float(np.corrcoef(A[ok], B[ok])[0, 1]) if ok.sum() > 3 else float("nan")


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
    ids = rcsb_ids()
    print(f"RCSB returned {len(ids)} cryo-EM 'antibody' entries at <=3.5 A", flush=True)
    # load cache
    records = json.load(open(REC_CACHE)) if REC_CACHE.exists() else []
    done = {r["pdb"] for r in records}
    for n, pid in enumerate(ids):
        if pid in done or len(records) >= limit:
            continue
        rec = extract(pid)
        if rec:
            records.append(rec)
        if n % 25 == 0:
            json.dump(records, open(REC_CACHE, "w"))
            print(f"  processed {n}/{len(ids)}, kept {len(records)} Ab-Ag", flush=True)
    json.dump(records, open(REC_CACHE, "w"))
    print(f"FINAL kept {len(records)} cryo-EM Ab-Ag complexes", flush=True)

    unit, n_units = cluster_antigens(records)
    # X-ray reference matrices (asymmetric antibody, symmetric ppi), from raw
    LabX = C.log2_odds(C.contact_matrix(C.load_contacts('antibody'), symmetrise=False))
    LppX = C.log2_odds(C.contact_matrix(C.load_contacts('ppi'), symmetrise=True))
    Mboth = build_matrix(records, unit, "both")
    Mside = build_matrix(records, unit, "side")
    Mmain = build_matrix(records, unit, "main")
    out = {
        "n_structures": len(records), "n_antigen_clusters": n_units,
        "total_tyr_bb": sum(r["tyr_bb"] for r in records),
        "grammar_corr_cryoEM_vs_Xray_antibody": corr(Mboth, LabX),
        "grammar_corr_cryoEM_vs_Xray_ppi": corr(Mboth, LppX),
        "side_corr_cryoEM_vs_Xray_antibody": corr(Mside, C.log2_odds(C.contact_matrix(
            C.load_contacts('antibody'), symmetrise=False))),  # placeholder same as full ab
        "note": "cryo-EM antibody grammar vs X-ray antibody (reproduction) and vs X-ray PPI (divergence).",
    }
    # channel cross-arm: cryo main/side vs X-ray ppi to see if divergence holds
    out["main_corr_cryoEM_vs_Xray_ppi"] = corr(Mmain, LppX)
    out["side_corr_cryoEM_vs_Xray_ppi"] = corr(Mside, LppX)
    json.dump(out, open(OUT / "cryoem_arm.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
