"""Download real DSSP annotations for every structure in the structural pass."""
import pathlib, sys, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C, phase6a_structural as P
import validate_dssp_real as V

ids = []
for arm in C.ARMS:
    ids += [(arm, p) for p in sorted(P.targets(arm))]
print(f"{len(ids)} structures", flush=True)
ok = miss = cached = 0
t0 = time.time()
for k, (arm, pdb) in enumerate(ids, 1):
    p = V.CACHE / f"{pdb.lower()}.cif"
    if p.exists() and p.stat().st_size > 1000:
        cached += 1
    elif V.fetch(pdb) is not None:
        ok += 1
    else:
        miss += 1
    if k % 200 == 0:
        el = time.time() - t0
        print(f"  {k}/{len(ids)}  new={ok} cached={cached} miss={miss}  "
              f"{el:5.0f}s  eta {el/k*(len(ids)-k)/60:4.1f} min", flush=True)
print(f"done: {ok} downloaded, {cached} cached, {miss} unavailable "
      f"({miss/len(ids):.1%})")
