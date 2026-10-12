"""add12 - does the tyrosine-backbone mechanism carry energetic signal? (validation).

Kills "the grammar is geometry, not energetics." For the antibody-antigen complexes in
SKEMPI 2.0 with a measured wild-type affinity, we download the structure, count the
tyrosine-hydroxyl -> antigen main-chain hydrogen bonds (the mechanism this paper localises),
and test whether that count tracks the measured binding free energy. If more tyrosine-backbone
bonds go with tighter binding, the grammar has energetic meaning, not merely geometric.

Independent of the structural corpus. SKEMPI affinity + RCSB structures.
"""
from __future__ import annotations
import csv, json, math, os, re, pathlib, urllib.request
import numpy as np
from scipy.stats import spearmanr
from Bio.PDB import PDBParser, NeighborSearch

ROOT = pathlib.Path(__file__).resolve().parent.parent
# SKEMPI 2.0 is a public download (https://life.bsc.es/pid/skempi2); not redistributed here.
# Set SKEMPI_CSV to its path, or drop skempi_v2.csv into metadata/skempi/.
SKEMPI = pathlib.Path(os.environ.get("SKEMPI_CSV", str(ROOT / "metadata" / "skempi" / "skempi_v2.csv")))
PDBDIR = ROOT / "metadata" / "skempi" / "pdb_cache"; PDBDIR.mkdir(parents=True, exist_ok=True)  # regenerable (gitignored)
OUT_JSON = ROOT / "metadata" / "analysis" / "affinity_mechanism.json"
OUT_MD   = ROOT / "metadata" / "analysis" / "AFFINITY_MECHANISM_RESULTS.md"
R = 1.987e-3
parser = PDBParser(QUIET=True)
AB_KW = ("antibody","fab","immunoglobulin","heavy chain","light chain","scfv","nanobody","vhh","igg","fv ")

def num(v):
    if v in ('','n.d.',None): return None
    m=re.search(r'-?\d+\.?\d*(?:[eE]-?\d+)?', str(v)); return float(m.group()) if m else None

def load_abag():
    rows=list(csv.DictReader(open(SKEMPI), delimiter=';'))
    seen={}
    for r in rows:
        if r['Hold_out_type']!='AB/AG': continue
        pid=r['#Pdb'].split('_')[0].upper()
        if pid in seen: continue
        kd=num(r['Affinity_wt_parsed']); T=num(r['Temperature']) or 298.0
        if not kd or kd<=0: continue
        parts=r['#Pdb'].split('_')
        if len(parts)!=3: continue
        c1,c2=set(parts[1]),set(parts[2])
        p1,p2=r['Protein 1'].lower(),r['Protein 2'].lower()
        # antibody side = the protein whose name matches antibody keywords
        if any(k in p1 for k in AB_KW):   ab,ag=c1,c2
        elif any(k in p2 for k in AB_KW): ab,ag=c2,c1
        else: ab,ag=(c1,c2) if len(c1)>=len(c2) else (c2,c1)  # fallback: 2-chain side
        seen[pid]={'pdb':pid,'ab':ab,'ag':ag,'dG':R*T*math.log(kd)}
    return list(seen.values())

def fetch(pid):
    f=PDBDIR/f"{pid}.pdb"
    if not f.exists():
        try: urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.pdb", f)
        except Exception: return None
    return f if f.exists() and f.stat().st_size>1000 else None

def analyse(entry):
    f=fetch(entry['pdb'])
    if not f: return None
    try: st=parser.get_structure(entry['pdb'], f)
    except Exception: return None
    model=next(iter(st))
    chains={c.id for c in model}
    ab=entry['ab']&chains; ag=entry['ag']&chains
    if not ab or not ag: return None
    atoms=[a for a in model.get_atoms() if a.get_parent().get_parent().id in (ab|ag)]
    ns=NeighborSearch(atoms)
    def side(a): return 'ab' if a.get_parent().get_parent().id in ab else 'ag'
    # tyrosine OH (antibody) -> antigen backbone N/O within 3.5 A
    tyr_bb=0; iface=set()
    for a in atoms:
        for b in ns.search(a.coord, 4.5):
            if side(a)==side(b): continue
            iface.add(tuple(sorted([id(a.get_parent()),id(b.get_parent())])))
        res=a.get_parent()
        if side(a)=='ab' and res.get_resname()=='TYR' and a.name=='OH':
            for b in ns.search(a.coord, 3.5):
                if side(b)=='ag' and b.name in ('N','O'):
                    tyr_bb+=1
    return {'pdb':entry['pdb'],'dG':entry['dG'],'tyr_bb':tyr_bb,'iface_pairs':len(iface),
            'tyr_bb_per100':100.0*tyr_bb/max(1,len(iface))}

def main():
    ents=load_abag()
    print(f"{len(ents)} AB/AG complexes with affinity; downloading + measuring...")
    res=[]
    for i,e in enumerate(ents):
        r=analyse(e)
        if r: res.append(r)
        if (i+1)%10==0: print(f"  {i+1}/{len(ents)}")
    dG=np.array([r['dG'] for r in res])
    tb=np.array([r['tyr_bb'] for r in res])
    tbn=np.array([r['tyr_bb_per100'] for r in res])
    rho1,p1=spearmanr(tb,dG); rho2,p2=spearmanr(tbn,dG)
    out={'n':len(res),
         'spearman_tyrBB_vs_dG':[round(rho1,3),round(p1,4)],
         'spearman_tyrBBper100_vs_dG':[round(rho2,3),round(p2,4)],
         'median_tyr_bb':float(np.median(tb)),
         'frac_with_tyr_bb':float(np.mean(tb>0)),
         'records':res}
    OUT_JSON.write_text(json.dumps(out,indent=1))
    L=["# Affinity vs the tyrosine-backbone mechanism (validation)\n",
       f"For {len(res)} antibody-antigen complexes in SKEMPI 2.0 with a measured wild-type "
       "affinity, the structure was fetched from the PDB and the number of tyrosine-hydroxyl "
       "-> antigen main-chain hydrogen bonds counted, then correlated with binding free "
       "energy (ΔG = RT ln Kd; more negative = tighter).\n",
       f"- Tyr-OH↔backbone count vs ΔG: Spearman rho = {rho1:+.2f} (p = {p1:.3f}).",
       f"- Count per 100 interface contacts vs ΔG: rho = {rho2:+.2f} (p = {p2:.3f}).",
       f"- {out['frac_with_tyr_bb']*100:.0f}% of complexes have >=1 such bond; median count "
       f"{np.median(tb):.0f}.\n",
       "A **negative** rho means more tyrosine-backbone hydrogen bonds accompany tighter "
       "binding — the mechanism the grammar localises carries measurable energetic signal, "
       "not merely a geometric preference."]
    OUT_MD.write_text("\n".join(L))
    print(f"n={len(res)}  rho(tyrBB,dG)={rho1:+.3f} p={p1:.3f}  rho(per100,dG)={rho2:+.3f} p={p2:.3f}")
    print("wrote", OUT_JSON.name)

if __name__=="__main__":
    main()
