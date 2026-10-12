"""Regenerate the matplotlib publication figures at 400 dpi with panel labels.

Reuses make_publication_figures' panel-letter/no-title savefig patch, then forces every
savefig to 400 dpi. Runs only the scripts that produce embedded figures. Deterministic:
recomputation reproduces the committed results.
"""
import runpy, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.figure
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_publication_figures as MPF

MPF.install_patches()                     # adds (a),(b)... panel letters; no-ops titles
_patched = matplotlib.figure.Figure.savefig
def savefig_400(self, *a, **k):
    k["dpi"] = 400
    return _patched(self, *a, **k)
matplotlib.figure.Figure.savefig = savefig_400

# only the scripts behind embedded Fig 1,2,3,5,7,8,10
SCRIPTS = ["phase2_inventory.py","aa_heatmap.py","phase5_readout.py",
           "phase6c_real_dssp.py","add5_epitope_logic.py","add4_pca_geometry.py",
           "case_thermodynamics.py"]
for name in SCRIPTS:
    p = HERE/name
    if not p.exists(): print("  MISSING",name); continue
    print("running",name,flush=True)
    try:
        runpy.run_path(str(p), run_name="__main__")
    except SystemExit:
        pass
    except Exception as e:
        print(f"   FAILED {type(e).__name__}: {e}")
    import matplotlib.pyplot as plt; plt.close("all")

# remontage the PyMOL figures at 400 dpi (montage of existing renders)
_orig = MPF._ORIGINAL_SAVEFIG
def montage_savefig(self,*a,**k):
    k["dpi"]=400; return _orig(self,*a,**k)
MPF._ORIGINAL_SAVEFIG = montage_savefig
try:
    MPF.remontage_all()
    print("remontaged PyMOL montages at 400 dpi")
except Exception as e:
    print("remontage failed:",e)
print("done")
