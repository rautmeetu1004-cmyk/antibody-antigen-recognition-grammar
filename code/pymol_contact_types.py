"""Render a gallery of interface contact types for the feature-extraction / forces section.

Six real examples from high-resolution antibody-antigen complexes, one per interaction type
the pipeline annotates. Antibody-side carbons in blue, antigen-side in wheat; the contact is a
yellow dashed line with its distance. Montaged 2 x 3.
"""
import os
from pymol import cmd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg

PDBDIR="pdb_antibody"
FIG="figures"
# (label, pdb, ab_chain, ab_resi, ab_atom, ag_chain, ag_resi, ag_atom, mode)
EX=[
 ("Salt bridge","26AP","B","101","NH2","A","215","OE2","atom"),
 ("Hydrogen bond","1A2Y","A","53","OG1","C","19","ND2","atom"),
 ("Van der Waals","10ZO","A","31","CG2","B","665","CD2","atom"),
 ("Aromatic stacking","1XFP","A","32",None,"L","62",None,"ring"),
 ("Cation-pi","1G7I","B","99","NH1","C","23",None,"ring"),
 ("Tyrosine OH to backbone","1A2Y","B","101","OH","C","121","N","atom"),
]

def render(label,pdb,abc,abr,aba,agc,agr,aga,mode,out):
    cmd.reinitialize()
    cmd.load(os.path.join(PDBDIR,f"{pdb}.pdb"),pdb)
    cmd.hide("everything")
    ab=f"{pdb} and chain {abc} and resi {abr}"
    ag=f"{pdb} and chain {agc} and resi {agr}"
    cmd.show("sticks",ab); cmd.show("sticks",ag)
    cmd.color("marine",f"({ab}) and elem C")
    cmd.color("wheat",f"({ag}) and elem C")
    cmd.util.cnc(f"({ab}) or ({ag})")
    cmd.set("stick_radius",0.16)
    # contact indicator
    if mode=="atom":
        a1=f"{ab} and name {aba}"; a2=f"{ag} and name {aga}"
        cmd.distance("dd",a1,a2)
    else:
        # ring centroid pseudo-atoms
        ring_ab="(%s and name CG+CD1+CD2+CE1+CE2+CZ+CE3+CZ2+CZ3+CH2+NE1)"%ab
        ring_ag="(%s and name CG+CD1+CD2+CE1+CE2+CZ+CE3+CZ2+CZ3+CH2+NE1)"%ag
        cmd.pseudoatom("c1",selection=ring_ab); cmd.pseudoatom("c2",selection=ring_ag)
        cmd.distance("dd","c1","c2")
        cmd.hide("everything","c1 or c2")
    cmd.set("dash_color","yellow"); cmd.set("dash_width",3.5)
    cmd.set("label_size",0); cmd.hide("labels")
    cmd.set("cartoon_transparency",1)
    cmd.orient(f"({ab}) or ({ag})")
    cmd.zoom(f"({ab}) or ({ag})",2.5)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background",1); cmd.set("ray_shadows",0)
    cmd.set("antialias",2)
    cmd.png(out,width=1200,height=1000,dpi=400,ray=1)

def main():
    panels=[]
    for label,pdb,abc,abr,aba,agc,agr,aga,mode in EX:
        out=os.path.join(FIG,f"contact_{label.split()[0].lower()}.png")
        try:
            render(label,pdb,abc,abr,aba,agc,agr,aga,mode,out)
            panels.append((out,label,pdb)); print("rendered",label,pdb,flush=True)
        except Exception as e:
            print("FAILED",label,e)
    # montage 2x3
    fig,axes=plt.subplots(2,3,figsize=(13.5,8))
    letters="abcdef"
    for k,(ax,(png,label,pdb)) in enumerate(zip(axes.ravel(),panels)):
        ax.imshow(mpimg.imread(png)); ax.axis("off")
        ax.set_title(f"({letters[k]}) {label}",fontsize=11,loc="left")
        ax.text(0.98,0.02,pdb,transform=ax.transAxes,ha="right",va="bottom",fontsize=8,color="0.4")
    for ax in axes.ravel()[len(panels):]: ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG,"contact_types.png"),dpi=400,bbox_inches="tight")
    print("montage -> figures/contact_types.png")

main()
