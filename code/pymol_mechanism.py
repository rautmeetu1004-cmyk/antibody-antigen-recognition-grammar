"""PyMOL: show the mechanism the statistics predicted.

Phase 6b found that the single most antibody-enriched functional-group contact is
the tyrosine phenol -OH hydrogen bonding to the ANTIGEN BACKBONE (+0.56 to N-H,
+0.58 to C=O), and Phase 5 found that backbone readout is the channel where antibody
and general-PPI grammars diverge most (cross-arm r = 0.135). Those are statistical
claims about 36,810 contacts. This renders the actual atoms.

Selected examples are the highest-resolution instances in the corpus, so the
geometry shown is real rather than modelled: every bond drawn is 2.6-3.1 A between
a Tyr OH and a main-chain N or O.

Run:  pymol -cq code/pymol_mechanism.py
Out:  figures/pymol_mechanism_<pdb>.png  and a combined montage
"""
import os
import sys

from pymol import cmd

# PyMOL rewrites __file__ to its own site-packages when running -cq scripts, so the
# study root is taken from the working directory (or AA_STUDY_ROOT) instead.
ROOT = os.environ.get("AA_STUDY_ROOT", os.getcwd())
PDB_DIR = os.path.join(ROOT, "pdb_antibody")
FIG = os.path.join(ROOT, "figures")

# (pdb, resolution, [(ab_chain, ab_resi, ag_chain, ag_resi, ag_atom, dist)])
EXAMPLES = [
    ("3D9A", 1.20, "Lysozyme C", [
        ("H", "350", "C", "700", "O", 2.64),
        ("H", "333", "C", "697", "O", 2.72)]),
    ("2P45", 1.10, "Ribonuclease pancreatic", [
        ("B", "27", "A", "60", "O", 2.68),
        ("B", "27", "A", "62", "N", 2.89)]),
    ("7AQZ", 1.30, "Nanobody complex", [
        ("D", "113", "A", "285", "N", 2.94)]),
]


def setup():
    cmd.set("ray_opaque_background", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("cartoon_transparency", 0.55)
    cmd.set("stick_radius", 0.13)
    cmd.set("dash_gap", 0.25)
    cmd.set("dash_width", 2.5)
    cmd.set("dash_color", "black")
    cmd.set("label_size", 16)
    cmd.set("label_color", "black")
    cmd.set("valence", 0)
    cmd.bg_color("white")


def render(pdb, resolution, antigen_name, bonds):
    path = os.path.join(PDB_DIR, pdb + ".pdb")
    if not os.path.exists(path):
        print("MISSING " + path)
        return None
    cmd.delete("all")
    cmd.load(path, pdb)
    cmd.remove("solvent or hydrogens")
    cmd.remove("not alt ''+A")            # keep the highest-occupancy conformer
    cmd.alter("all", "alt=''")

    ab_chains = sorted({b[0] for b in bonds})
    ag_chains = sorted({b[2] for b in bonds})
    ab_sel = " or ".join("chain " + c for c in ab_chains)
    ag_sel = " or ".join("chain " + c for c in ag_chains)

    cmd.hide("everything")
    cmd.show("cartoon", ab_sel)
    cmd.show("cartoon", ag_sel)
    cmd.color("skyblue", ab_sel)
    cmd.color("wheat", ag_sel)

    for i, (abc, abr, agc, agr, agatom, dist) in enumerate(bonds):
        tyr = "%s_tyr%d" % (pdb, i)
        tgt = "%s_tgt%d" % (pdb, i)
        cmd.select(tyr, "chain %s and resi %s and sidechain" % (abc, abr))
        # the ANTIGEN BACKBONE atom -- the point of the figure
        cmd.select(tgt, "chain %s and resi %s and name N+CA+C+O" % (agc, agr))
        cmd.show("sticks", tyr)
        cmd.show("sticks", tgt)
        cmd.color("marine", tyr + " and elem C")
        cmd.color("orange", tgt + " and elem C")
        cmd.util.cnc(tyr)
        cmd.util.cnc(tgt)
        cmd.distance("%s_hb%d" % (pdb, i),
                     "chain %s and resi %s and name OH" % (abc, abr),
                     "chain %s and resi %s and name %s" % (agc, agr, agatom))
        # label the antigen residue too, so the figure shows WHAT is being read
        cmd.label("chain %s and resi %s and name OH" % (abc, abr),
                  '"Tyr%s"' % abr)
        cmd.label("chain %s and resi %s and name CA" % (agc, agr),
                  '"%s%s (backbone)"' % (agatom, agr))

    # show side chains around the site for context, then frame ON THE BOND --
    # zooming to "visible" framed the whole complex and hid the actual chemistry
    site = " or ".join("%s_tyr%d or %s_tgt%d" % (pdb, i, pdb, i)
                       for i in range(len(bonds)))
    cmd.select("nearby", "byres (polymer within 6 of (%s))" % site)
    cmd.show("sticks", "nearby and sidechain and not (elem H)")
    cmd.set("stick_transparency", 0.6, "nearby and not (%s)" % site)
    cmd.orient(site)
    cmd.zoom(site, 4.5)
    out = os.path.join(FIG, "pymol_mechanism_%s.png" % pdb)
    cmd.set("ray_trace_mode", 0)
    cmd.png(out, width=1500, height=1150, dpi=200, ray=1)
    print("wrote %s  (%s, %.2f A, %d bonds)" % (out, antigen_name, resolution,
                                                len(bonds)))
    return out


if __name__ == "__main__" or True:
    setup()
    made = []
    for pdb, resolution, name, bonds in EXAMPLES:
        p = render(pdb, resolution, name, bonds)
        if p:
            made.append((p, pdb, resolution, name, bonds))

    # montage with matplotlib so the panel labels carry the statistics
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg

        n = len(made)
        if n:
            fig, axes = plt.subplots(1, n, figsize=(5.4 * n, 5.8))
            if n == 1:
                axes = [axes]
            for ax, (p, pdb, resolution, name, bonds) in zip(axes, made):
                ax.imshow(mpimg.imread(p))
                ax.axis("off")
                ds = ", ".join("%.2f $\\AA$" % b[5] for b in bonds)
                ax.set_title("%s  (%.2f $\\AA$)\n%s\nTyr OH $\\rightarrow$ "
                             "backbone: %s" % (pdb, resolution, name, ds),
                             fontsize=10)
            fig.suptitle("Tyrosine phenol -OH hydrogen bonding to the antigen "
                         "BACKBONE\nthe most antibody-enriched functional-group "
                         "contact (+0.56 to N-H, +0.58 to C=O vs general PPI)",
                         fontsize=12.5)
            fig.tight_layout(rect=[0, 0, 1, 0.88])
            out = os.path.join(FIG, "pymol_mechanism_montage.png")
            fig.savefig(out, dpi=170, bbox_inches="tight")
            print("wrote " + out)
    except Exception as e:
        print("montage failed: %s: %s" % (type(e).__name__, e))
