"""PyMOL: paint the immunodominance profile onto the lysozyme surface.

The bar-chart version of this result (immunodominance_profile.png) shows WHICH
positions 42 independent antibodies converge on. It cannot show whether those
positions form one contiguous patch on the folded protein or several separate ones --
which is the question that matters for whether the immune response is focused on a
site or merely on a set of exposed residues.

Contact frequency is written into the B-factor column of a high-resolution reference
structure (3D9A, 1.20 A) and rendered as a surface, so colour = fraction of the 42
antibodies contacting that residue. Catalytic residues are marked, because the
lysozyme case study found antibodies converge on the RIM of the active-site cleft
(Asp101, 69%) rather than on the catalytic pair (Glu35 12%, Asp52 19%).

Run:  AA_STUDY_ROOT="$PWD" pymol -cq code/pymol_immunodominance.py
Out:  figures/pymol_immunodominance.png
"""
import json
import os

from pymol import cmd

ROOT = os.environ.get("AA_STUDY_ROOT", os.getcwd())
FIG = os.path.join(ROOT, "figures")
PROF = os.path.join(ROOT, "metadata", "analysis", "immunodominance_profile.json")

CATALYTIC = [35, 52, 101]          # mature numbering: Glu35, Asp52, Asp101


def main():
    with open(PROF) as fh:
        d = json.load(fh)
    prof = {int(k): v for k, v in d["profile"].items()}
    pdb = d["reference_pdb"]
    chain = d["reference_chain"]
    offset = d["offset"]
    n_ab = d["n_antibodies"]

    path = os.path.join(ROOT, "pdb_antibody", pdb + ".pdb")
    cmd.delete("all")
    cmd.load(path, pdb)
    cmd.remove("solvent or hydrogens")
    cmd.remove("not alt ''+A")
    cmd.alter("all", "alt=''")

    # keep only the antigen chain: the antibody would occlude the very surface we
    # are trying to show
    cmd.create("ag", f"{pdb} and chain {chain} and polymer")
    cmd.delete(pdb)

    # contact frequency -> B-factor
    cmd.alter("ag", "b=0.0")
    n_set = 0
    for pos, frac in prof.items():
        resi = pos + offset - 1
        n = cmd.alter(f"ag and resi {resi}", f"b={frac:.4f}")
        if n:
            n_set += 1
    print("positions written to B-factor: %d" % n_set)

    cmd.hide("everything")
    cmd.show("surface", "ag")
    cmd.set("surface_quality", 1)
    cmd.spectrum("b", "white_yellow_orange_red", "ag", minimum=0.0, maximum=1.0)
    cmd.set("ray_opaque_background", 0)
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)

    # mark the catalytic residues
    for p in CATALYTIC:
        r = p + offset - 1
        sel = f"cat{p}"
        cmd.select(sel, f"ag and resi {r}")
        cmd.show("sticks", f"{sel} and sidechain")
        cmd.color("marine", f"{sel} and elem C")

    cmd.orient("ag")
    views = {"front": (1, 0, 0, 0), "back": None}
    out = os.path.join(FIG, "pymol_immunodominance.png")
    cmd.png(out, width=1500, height=1300, dpi=200, ray=1)
    print("wrote %s (n=%d antibodies)" % (out, n_ab))

    cmd.turn("y", 180)
    out2 = os.path.join(FIG, "pymol_immunodominance_back.png")
    cmd.png(out2, width=1500, height=1300, dpi=200, ray=1)
    print("wrote %s" % out2)

    # two-panel montage with a colour key
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from matplotlib import cm, colors

        fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.0))
        for ax, p, lab in ((axes[0], out, "front"), (axes[1], out2, "back (180 deg)")):
            ax.imshow(mpimg.imread(p))
            ax.axis("off")
            ax.set_title(lab, fontsize=10)
        sm = cm.ScalarMappable(norm=colors.Normalize(0, 1),
                               cmap=colors.LinearSegmentedColormap.from_list(
                                   "wyor", ["white", "yellow", "orange", "red"]))
        cb = fig.colorbar(sm, ax=axes, fraction=0.03, pad=0.02)
        cb.set_label(f"fraction of {n_ab} independent antibodies contacting",
                     fontsize=9)
        fig.suptitle("Immunodominance on lysozyme C (3D9A, 1.20 $\\AA$)\n"
                     "blue sticks = catalytic Glu35 / Asp52 / Asp101", fontsize=12)
        fig.savefig(os.path.join(FIG, "pymol_immunodominance_montage.png"),
                    dpi=170, bbox_inches="tight")
        print("wrote montage")
    except Exception as e:
        print("montage failed: %s: %s" % (type(e).__name__, e))


main()
