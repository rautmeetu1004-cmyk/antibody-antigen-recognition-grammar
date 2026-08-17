"""PyMOL: the same protein read by an antibody and by its natural partner.

The matched-target analysis reports a median Jaccard overlap of 0.03 between the
epitope and the natural-partner site across 16 proteins solved in both complexes.
That number is abstract. This renders the comparison directly, painting both
footprints onto one copy of the shared protein so that the spatial relationship
between them can be seen rather than inferred.

Two pairs are shown, chosen from opposite ends of the overlap distribution:
angiogenin, where the two footprints partly coincide, and a low-overlap pair where
they do not. Numbering is reconciled by pairwise alignment between the two
structures' copies of the protein, exactly as in the quantitative analysis.

Run:  AA_STUDY_ROOT="$PWD" pymol -cq code/pymol_matched_pair.py
Out:  figures/pymol_matched_pair.png
"""
import json
import os
import sys

from pymol import cmd

ROOT = os.environ.get("AA_STUDY_ROOT", os.getcwd())
FIG = os.path.join(ROOT, "figures")
sys.path.insert(0, os.path.join(ROOT, "code"))

# (uniprot, antibody pdb, ppi pdb, label)
PAIRS = [
    ("P03950", "1H0D", "1A4Y", "Angiogenin"),
    ("P14190", "4WEM", "2J6R", "Low-overlap pair"),
]


def load_positions():
    """Epitope and partner-site positions per matched pair, from the analysis."""
    import case_matched_targets as MT
    contacts = {"antibody": {}, "ppi": {}}
    for arm in ("antibody", "ppi"):
        acc = {}
        for r in MT.C.dedup_rows(arm, MT.C.load_contacts(arm)):
            acc.setdefault((r["pdb_id"], r["ag_chain"]), set()).add(
                r["ag_resi"].strip())
        contacts[arm] = acc
    return contacts, MT


def render(uniprot, ab_pdb, pp_pdb, label, contacts, MT, out_png):
    """Paint the epitope and the natural-partner site onto the PPI copy."""
    idx = MT.antigen_uniprot()
    entry = idx.get(uniprot)
    if not entry:
        print("no entry for " + uniprot)
        return None
    ab = next((x for x in entry["antibody"] if x[0] == ab_pdb), None)
    pp = next((x for x in entry["ppi"] if x[0] == pp_pdb), None)
    if not ab or not pp:
        print("missing chain record for %s" % uniprot)
        return None

    sa, ka = MT.chain_seq_and_positions(ab[0], ab[1], "antibody")
    sb, kb = MT.chain_seq_and_positions(pp[0], pp[1], "ppi")
    m = MT.align_index(sa, sb)
    if not m:
        print("alignment failed for " + uniprot)
        return None

    inv_kb = {v: k for k, v in kb.items()}
    epi_keys, site_keys = set(), set()
    for x in contacts["antibody"].get(ab, set()):
        if x in ka and ka[x] in m:
            r = inv_kb.get(m[ka[x]])
            if r:
                epi_keys.add(r)
    for x in contacts["ppi"].get(pp, set()):
        if x in kb:
            site_keys.add(x)

    both = epi_keys & site_keys
    epi_only = epi_keys - both
    site_only = site_keys - both

    path = os.path.join(ROOT, "pdb_ppi", pp_pdb + ".pdb")
    cmd.delete("all")
    cmd.load(path, pp_pdb)
    cmd.remove("solvent or hydrogens")
    cmd.remove("not alt ''+A")
    cmd.alter("all", "alt=''")
    cmd.create("shared", "%s and chain %s and polymer" % (pp_pdb, pp[1]))
    cmd.delete(pp_pdb)

    cmd.hide("everything")
    cmd.show("surface", "shared")
    cmd.set("surface_quality", 1)
    cmd.color("grey80", "shared")

    def paint(keys, colour, name):
        if not keys:
            return
        sel = " or ".join("resi %s" % k for k in sorted(keys))
        cmd.select(name, "shared and (%s)" % sel)
        cmd.color(colour, name)

    paint(site_only, "skyblue", "site_only")      # natural partner only
    paint(epi_only, "firebrick", "epi_only")      # antibody only
    paint(both, "violetpurple", "both")           # contacted by both

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.orient("shared")
    cmd.png(out_png, width=1400, height=1250, dpi=200, ray=1)
    print("wrote %s  (%s: epitope-only %d, partner-only %d, both %d)"
          % (out_png, label, len(epi_only), len(site_only), len(both)))
    return {"label": label, "uniprot": uniprot, "antibody": ab_pdb, "ppi": pp_pdb,
            "epitope_only": len(epi_only), "partner_only": len(site_only),
            "both": len(both)}


def main():
    contacts, MT = load_positions()
    made = []
    for uniprot, ab_pdb, pp_pdb, label in PAIRS:
        out = os.path.join(FIG, "pymol_matched_%s.png" % uniprot)
        r = render(uniprot, ab_pdb, pp_pdb, label, contacts, MT, out)
        if r:
            made.append((out, r))

    if not made:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from matplotlib.patches import Patch

        fig, axes = plt.subplots(1, len(made), figsize=(6.0 * len(made), 6.0))
        if len(made) == 1:
            axes = [axes]
        for k, (ax, (png, r)) in enumerate(zip(axes, made)):
            ax.imshow(mpimg.imread(png))
            ax.axis("off")
            ax.text(0.02, 0.97, "(%s)" % "abcdefgh"[k], transform=ax.transAxes,
                    fontsize=15, fontweight="bold", va="top")
        handles = [Patch(facecolor="firebrick", label="contacted by antibody only"),
                   Patch(facecolor="skyblue", label="contacted by natural partner only"),
                   Patch(facecolor="mediumpurple", label="contacted by both")]
        fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9,
                   frameon=False)
        fig.tight_layout(rect=[0, 0.06, 1, 1])
        out = os.path.join(FIG, "pymol_matched_pair.png")
        fig.savefig(out, dpi=170, bbox_inches="tight")
        print("wrote " + out)
    except Exception as e:
        print("montage failed: %s: %s" % (type(e).__name__, e))


main()
