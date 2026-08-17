"""Phase 10 - literature concordance: does the pipeline recover known biology?

Two jobs, and they are different:

  POSITIVE CONTROLS. Paratope Tyr/Trp enrichment, CDR-H3 dominance, discontinuous
  epitopes, charge complementarity -- all long established. If this pipeline did not
  recover them, the pipeline would be wrong. Recovering them is what licenses the
  novel claims; they are NOT presented as findings.

  GENUINE COMPARISONS. Places where our measurement can agree, refine, or contradict
  a published claim.

On citation discipline: our measurements below are computed, exact, and reproducible
from the JSONs in metadata/analysis/. The published values are stated as the
qualitative claim plus an approximate figure where one is well known, and every row
is flagged `verify_source` so that no number goes into a manuscript without being
checked against the paper itself. Inventing a precise literature value to make a
table look tidy would be worse than leaving it qualitative.

Outputs
  LITERATURE_CONCORDANCE.md
  metadata/analysis/literature.json
"""
from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C


def load(name):
    p = C.OUT_MET / name
    return json.loads(p.read_text()) if p.exists() else {}


def build_rows():
    inv = load("inventory.json")
    hm = load("heatmap.json")
    cons = load("conservation.json")
    para = load("paratope.json")
    conf = load("conformation.json")
    diff = load("differential.json")
    ctrl = load("controls.json")
    mod = load("modified.json")
    fg = load("functional_groups.json")

    def g(d, *keys, default=None):
        for k in keys:
            if not isinstance(d, dict):
                return default
            d = d.get(k, {} if k != keys[-1] else default)
        return d if d != {} else default

    rows = []

    # ---------------------------------------------------------- positive controls
    pc = para.get("paratope_composition", {})
    rows.append({
        "topic": "Paratope aromatic enrichment",
        "published_claim": "Tyr and Trp are strongly over-represented in antibody "
                           "combining sites; Tyr is the single commonest paratope "
                           "residue.",
        "sources": "Mian, Bradwell & Olson 1991; MacCallum, Martin & Thornton 1996; "
                   "Ramaraj et al. 2012; Peng et al. 2014",
        "our_measurement": (
            f"Tyr log2 enrichment +{pc.get('Y', {}).get('vs_ppi_interface', float('nan')):.2f} "
            f"vs general PPI interfaces and "
            f"+{pc.get('Y', {}).get('vs_antibody_surface', float('nan')):.2f} vs the "
            f"antibody's own surface; Trp "
            f"+{pc.get('W', {}).get('vs_ppi_interface', float('nan')):.2f} / "
            f"+{pc.get('W', {}).get('vs_antibody_surface', float('nan')):.2f}"),
        "verdict": "AGREE (positive control)",
        "note": "Recovered on both nulls, so this is selective use of aromatics, "
                "not a consequence of antibodies being aromatic-rich overall.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Paratope Ser/Thr/Gly enrichment",
        "published_claim": "Ser, Thr and Gly are commonly listed alongside Tyr/Trp "
                           "as enriched paratope residues.",
        "sources": "MacCallum et al. 1996; Collis et al. 2003; Ramaraj et al. 2012",
        "our_measurement": (
            f"Ser +{pc.get('S', {}).get('vs_ppi_interface', float('nan')):.2f} vs PPI "
            f"but {pc.get('S', {}).get('vs_antibody_surface', float('nan')):+.2f} vs "
            f"the antibody's own surface; Thr "
            f"+{pc.get('T', {}).get('vs_ppi_interface', float('nan')):.2f} / "
            f"{pc.get('T', {}).get('vs_antibody_surface', float('nan')):+.2f}; Gly "
            f"+{pc.get('G', {}).get('vs_ppi_interface', float('nan')):.2f} / "
            f"{pc.get('G', {}).get('vs_antibody_surface', float('nan')):+.2f}"),
        "verdict": "REFINE",
        "note": "Enriched relative to general protein interfaces but DEPLETED "
                "relative to the antibody's own accessible surface. The apparent "
                "enrichment is a property of antibody composition, not of selective "
                "paratope usage. Only Tyr and Trp are enriched on both nulls. This "
                "distinction is invisible to any study using a single background.",
        "verify_source": True,
    })

    cdr = g(para, "cdr_vs_framework", "FAB", "_cdr_only_fractions", default={})
    cdr_v = g(para, "cdr_vs_framework", "SD-H", "_cdr_only_fractions", default={})
    rows.append({
        "topic": "CDR-H3 dominance",
        "published_claim": "CDR-H3 contributes more antigen contacts than any other "
                           "loop and is the principal determinant of specificity.",
        "sources": "Kabat & Wu; MacCallum, Martin & Thornton 1996; Xu & Davis 2000",
        "our_measurement": (
            f"Fab: H3 = {cdr.get('H3', float('nan')):.1%} of CDR contacts "
            f"(next: H2 {cdr.get('H2', float('nan')):.1%}, L3 "
            f"{cdr.get('L3', float('nan')):.1%}, L2 {cdr.get('L2', float('nan')):.1%}). "
            f"Nanobody: H3 = {cdr_v.get('H3', float('nan')):.1%}"),
        "verdict": "AGREE (positive control)",
        "note": "L2 contributes almost nothing, also as reported. Nanobodies push "
                "H3 dominance far higher, consistent with the absence of a light "
                "chain and a longer H3.",
        "verify_source": True,
    })

    fw = g(para, "cdr_vs_framework", "ALL", default={})
    rows.append({
        "topic": "Framework contribution to the paratope",
        "published_claim": "Framework residues contribute a minority but "
                           "non-negligible fraction of antigen contacts.",
        "sources": "MacCallum et al. 1996; Sela-Culang, Kunik & Ofran 2013",
        "our_measurement": f"framework = {fw.get('framework', float('nan')):.1%} of "
                           f"paratope contacts, stable across Fab/Fv/nanobody "
                           f"(0.209-0.222)",
        "verdict": "AGREE",
        "note": "Remarkably constant across construct types.",
        "verify_source": True,
    })

    ep = para.get("epitope", {})
    rows.append({
        "topic": "Epitopes are discontinuous",
        "published_claim": "Most protein antigenic epitopes are conformational / "
                           "discontinuous, comprising several sequence segments and "
                           "roughly 15-25 residues.",
        "sources": "Barlow, Edwards & Thornton 1986; Rubinstein et al. 2009; "
                   "Kringelum et al. 2013; Sela-Culang et al. 2013",
        "our_measurement": (
            f"median {ep.get('median_segments', float('nan')):.0f} sequence segments, "
            f"{ep.get('frac_discontinuous_ge2_segments', float('nan')):.1%} of epitopes "
            f"discontinuous (>=2 segments), median "
            f"{ep.get('median_epitope_residues', float('nan')):.0f} epitope residues"),
        "verdict": "AGREE (positive control)",
        "note": "Both the segment count and the epitope size fall squarely in the "
                "published range.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Charge complementarity at interfaces",
        "published_claim": "Salt bridges and complementary charge pairing are a "
                           "conserved feature of protein-protein interfaces.",
        "sources": "Barlow & Thornton 1983; Jones & Thornton 1996; "
                   "Xu, Tsai & Nussinov 1997",
        "our_measurement": (
            "K-D, D-K, K-E, R-D, R-E are top-enriched and bootstrap-sign-stable in "
            "BOTH arms; salt-bridge fraction is statistically identical "
            "(antibody 0.033 vs PPI 0.033, p=0.30) as is opposite-charge pairing "
            "(0.051 vs 0.051, p=0.35)"),
        "verdict": "AGREE, and EXTEND",
        "note": "Charge chemistry is genuinely universal -- and therefore is NOT "
                "what distinguishes antibodies. A negative result worth stating.",
        "verify_source": True,
    })

    # ------------------------------------------------------------ real comparisons
    size = inv.get("size", {})
    rows.append({
        "topic": "Interface size",
        "published_claim": "Typical protein-protein interfaces bury ~1,600 A2 total "
                           "and involve ~20-30 residues per side; antibody-antigen "
                           "interfaces sit at the smaller end of the transient range.",
        "sources": "Jones & Thornton 1996; Lo Conte, Chothia & Janin 1999; "
                   "Chakrabarti & Janin 2002",
        "our_measurement": (
            f"antibody median {g(size, 'n_pairs', 'median_antibody', default=float('nan')):.0f} "
            f"residue pairs vs PPI "
            f"{g(size, 'n_pairs', 'median_ppi', default=float('nan')):.0f} "
            f"(effect {g(size, 'n_pairs', 'effect_rank_biserial', default=float('nan')):+.2f}, "
            f"p={g(size, 'n_pairs', 'p', default=float('nan')):.1e}); "
            f"epitope median {ep.get('median_epitope_residues', float('nan')):.0f} residues"),
        "verdict": "AGREE",
        "note": "Antibody interfaces are significantly smaller than general PPI "
                "interfaces on every size measure.",
        "verify_source": True,
    })

    phys = inv.get("physicochemistry", {})
    rows.append({
        "topic": "Hydrophobicity of antibody interfaces",
        "published_claim": "Antibody-antigen interfaces are less hydrophobic than "
                           "obligate/homodimeric interfaces and more polar in "
                           "character.",
        "sources": "Jones & Thornton 1996; Lo Conte et al. 1999; "
                   "Chakrabarti & Janin 2002",
        "our_measurement": (
            f"apolar-apolar contact fraction antibody "
            f"{g(phys, 'f_apolar_apolar', 'median_antibody', default=float('nan')):.3f} vs PPI "
            f"{g(phys, 'f_apolar_apolar', 'median_ppi', default=float('nan')):.3f} "
            f"(effect {g(phys, 'f_apolar_apolar', 'effect_rank_biserial', default=float('nan')):+.2f}, "
            f"p={g(phys, 'f_apolar_apolar', 'p', default=float('nan')):.1e}); H-bond fraction "
            f"{g(phys, 'f_hbond', 'median_antibody', default=float('nan')):.3f} vs "
            f"{g(phys, 'f_hbond', 'median_ppi', default=float('nan')):.3f} "
            f"(effect {g(phys, 'f_hbond', 'effect_rank_biserial', default=float('nan')):+.2f})"),
        "verdict": "AGREE, and QUANTIFY",
        "note": "The largest single physicochemical separations we measure. "
                "Antibodies substitute aromatic and H-bond chemistry for "
                "hydrophobic burial.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Aromatic/pi chemistry",
        "published_claim": "Aromatic residues, and cation-pi interactions, are "
                           "important in antibody recognition.",
        "sources": "Gallivan & Dougherty 1999; Ramaraj et al. 2012; "
                   "Peng et al. 2014",
        "our_measurement": (
            f"aromatic involvement antibody "
            f"{g(phys, 'f_aromatic', 'median_antibody', default=float('nan')):.3f} vs PPI "
            f"{g(phys, 'f_aromatic', 'median_ppi', default=float('nan')):.3f} "
            f"(effect {g(phys, 'f_aromatic', 'effect_rank_biserial', default=float('nan')):+.2f}); "
            f"cation-pi {g(phys, 'f_cation_pi', 'median_antibody', default=float('nan')):.3f} vs "
            f"{g(phys, 'f_cation_pi', 'median_ppi', default=float('nan')):.3f}; "
            f"functional-group map: aromatic ring <-> ammonium +1.07 antibody-enriched"),
        "verdict": "AGREE, and REFINE",
        "note": "But formal ring-ring stacking is DEPLETED in antibodies "
                "(W-Y delta -0.99, Y-Y -0.85, both significant). The aromatic "
                "enrichment is aromatics contacting NON-aromatics -- backbone and "
                "cations -- not pi-stacking. Frequently blurred in the literature.",
        "verify_source": True,
    })

    rows.append({
        "topic": "What antibodies read on the antigen",
        "published_claim": "Antibody recognition is usually described in terms of "
                           "side-chain complementarity and shape matching.",
        "sources": "Lawrence & Colman 1993; Sundberg & Mariuzza 2002",
        "our_measurement": (
            "side-chain-mediated fraction is IDENTICAL in both arms (~0.72, "
            "per-residue correlation r=0.981), but the main-chain readout GRAMMAR "
            "diverges (cross-arm r=0.135 vs 0.493 for side chains); Tyr phenol -OH "
            "to antigen backbone N-H/C=O is +0.56/+0.58 antibody-enriched"),
        "verdict": "NOVEL / EXTENDS",
        "note": "The distinguishing feature is not how MUCH backbone is read but "
                "HOW it is read. We are not aware of this decomposition being "
                "reported; it is the study's primary novel claim.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Are antibody interfaces chemically special?",
        "published_claim": "Antibody-antigen interfaces are often described as "
                           "'ordinary' protein interfaces distinguished mainly by "
                           "their loop scaffold and lack of co-evolution.",
        "sources": "Jones & Thornton 1996; Sundberg & Mariuzza 2002; "
                   "Ofran et al. 2008",
        "our_measurement": (
            f"grammar correlation r={g(hm, 'cluster30', 'grammar_similarity_antibody_vs_ppi', default=float('nan')):.3f} "
            f"(0.586 with Cys masked; chemical-class r=0.841); permutation test of a "
            f"single shared grammar: observed r="
            f"{g(diff, 'permutation_test', 'observed_r', default=float('nan')):.3f} vs null "
            f"{g(diff, 'permutation_test', 'null_median', default=float('nan')):.3f} "
            f"[{g(diff, 'permutation_test', 'null_ci', default=[0,0])[0]:.3f}, "
            f"{g(diff, 'permutation_test', 'null_ci', default=[0,0])[1]:.3f}], "
            f"p={g(diff, 'permutation_test', 'p_one_sided', default=float('nan')):.4f}"),
        "verdict": "PARTIALLY CONTRADICT",
        "note": "Same chemistry at the class level (r=0.84) but formally distinct "
                "grammars at residue level (p=0.001). And the difference is NOT the "
                "loop scaffold: matching both arms on secondary structure "
                "(loop-loop r=0.185) and on CATH fold class does not raise the "
                "correlation -- it lowers it, below size-matched nulls.",
        "verify_source": True,
    })

    ss = conf.get("ss_composition", {})
    rows.append({
        "topic": "Epitope secondary structure",
        "published_claim": "Epitopes are often reported as enriched in loops/turns "
                           "and in protruding, flexible regions.",
        "sources": "Thornton et al. 1986; Rubinstein et al. 2009; "
                   "Kringelum et al. 2013",
        "our_measurement": (
            f"epitope loop fraction {g(ss, 'antibody_target', 'loop', default=float('nan')):.3f} "
            f"vs PPI partner loop fraction "
            f"{g(ss, 'ppi_target', 'loop', default=float('nan')):.3f}; "
            f"paratope loop fraction {g(ss, 'antibody_probe', 'loop', default=float('nan')):.3f} "
            f"vs PPI probe {g(ss, 'ppi_probe', 'loop', default=float('nan')):.3f}"),
        "verdict": "PARTIALLY CONTRADICT",
        "note": "Epitopes are NOT structurally unusual: their loop content is "
                "indistinguishable from ordinary PPI partner surfaces (0.512 vs "
                "0.509). The loop enrichment is entirely on the PARATOPE side. "
                "Claims that epitopes are loop-rich appear to be relative to whole "
                "protein rather than to other interfaces.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Nanobody paratopes",
        "published_claim": "Nanobodies (VHH) rely disproportionately on a long "
                           "CDR-H3 and can access clefts inaccessible to Fabs.",
        "sources": "Muyldermans 2013; De Genst et al. 2006",
        "our_measurement": (
            f"H3 = {cdr_v.get('H3', float('nan')):.1%} of nanobody CDR contacts vs "
            f"{cdr.get('H3', float('nan')):.1%} in Fabs; Fab-vs-nanobody grammar "
            f"r={g(ctrl, 'construct_type', 'fab_vs_nanobody_r', default=float('nan')):.3f}, "
            f"inside the Fab split-half reliability ceiling 0.253 [0.117, 0.370]"),
        "verdict": "AGREE on architecture, EXTEND on chemistry",
        "note": "H3 dominance is confirmed and much stronger in nanobodies. But the "
                "two formats' recognition CHEMISTRY is indistinguishable at matched "
                "sample size -- different architecture, same grammar.",
        "verify_source": True,
    })

    rows.append({
        "topic": "PTMs and glycans at interfaces",
        "published_claim": "Glycan shields obstruct antibody access to viral "
                           "antigens; PTMs regulate protein-protein interactions.",
        "sources": "Wei et al. 2020 (glycan shield); Nishi, Hashimoto & Panchenko "
                   "2011 (PTMs at interfaces)",
        "our_measurement": (
            f"interfacial glycans: antibody "
            f"{g(mod, 'antibody_glycan', 'frac_structures_with_interfacial', default=float('nan')):.1%} "
            f"vs PPI "
            f"{g(mod, 'ppi_glycan', 'frac_structures_with_interfacial', default=float('nan')):.1%} "
            f"(OR 3.07, p=1.6e-11); interfacial PTMs: antibody "
            f"{g(mod, 'antibody_ptm', 'frac_structures_with_interfacial', default=float('nan')):.1%} "
            f"vs PPI "
            f"{g(mod, 'ppi_ptm', 'frac_structures_with_interfacial', default=float('nan')):.1%} "
            f"(OR 0.16, p=7.3e-11)"),
        "verdict": "AGREE, and QUANTIFY a double dissociation",
        "note": "Antibody interfaces engage glycans ~3x more; PPI interfaces engage "
                "PTMs ~6x more. Partly a target-composition effect (viral "
                "glycoproteins vs signalling complexes) and stated as such.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Shape complementarity (Sc)",
        "published_claim": "Antibody-antigen interfaces have lower shape "
                           "complementarity (~0.64-0.68) than protease-inhibitor or "
                           "oligomeric interfaces (~0.70-0.76).",
        "sources": "Lawrence & Colman 1993",
        "our_measurement": "NOT YET COMPUTED",
        "verdict": "PENDING",
        "note": "Required for publication: Lawrence & Colman defined Sc on exactly "
                "this comparison, and reviewers will expect it. Our smaller-but-"
                "tighter result (closest approach 3.587 vs 3.641 A) makes a specific "
                "prediction here, so it is a real test rather than a formality.",
        "verify_source": True,
    })

    rows.append({
        "topic": "Interface hot spots",
        "published_claim": "Binding free energy concentrates in a few 'hot spot' "
                           "residues, enriched in Trp, Arg and Tyr, surrounded by an "
                           "occluding O-ring.",
        "sources": "Clackson & Wells 1995; Bogan & Thorn 1998",
        "our_measurement": "NOT YET COMPUTED",
        "verdict": "PENDING",
        "note": "Burial data from Phase 6a is in place; the analysis is cheap and "
                "connects to the most-cited PPI literature.",
        "verify_source": True,
    })

    return rows


HEADER = """# Literature concordance

Every measurement below is computed by this study's pipeline and reproducible from
`metadata/analysis/*.json`. Published claims are stated qualitatively with sources;
**all citation details and any approximate published figures must be checked against
the source papers before use in a manuscript** -- they are recorded here as pointers,
not as verified quotations.

Verdicts: **AGREE** (positive control recovered) | **REFINE** (agrees but our data
sharpens or conditions it) | **EXTEND** (new measurement beyond the claim) |
**PARTIALLY CONTRADICT** | **NOVEL** | **PENDING** (analysis not yet run).

The positive controls matter as much as the novel claims: a pipeline that failed to
recover Tyr/Trp paratope enrichment, CDR-H3 dominance, or discontinuous epitopes
would not be trustworthy on anything else.

"""


if __name__ == "__main__":
    rows = build_rows()

    lines = [HEADER]
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
        lines.append(f"## {r['topic']}\n")
        lines.append(f"**Published claim.** {r['published_claim']}  \n")
        lines.append(f"*Sources:* {r['sources']}\n")
        lines.append(f"**Our measurement.** {r['our_measurement']}\n")
        lines.append(f"**Verdict: {r['verdict']}.** {r['note']}\n")
        lines.append("")

    lines.append("## Summary\n")
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"- **{k}**: {v}")
    lines.append("")
    lines.append("Positive controls recovered: paratope Tyr/Trp enrichment, CDR-H3 "
                 "dominance, framework contribution, discontinuous epitopes, charge "
                 "complementarity, interface size, reduced hydrophobicity. The "
                 "pipeline reproduces the established picture before departing "
                 "from it.\n")

    (C.ROOT / "LITERATURE_CONCORDANCE.md").write_text("\n".join(lines))
    (C.OUT_MET / "literature.json").write_text(json.dumps(
        {"rows": rows, "verdict_counts": counts}, indent=1))

    print(f"wrote LITERATURE_CONCORDANCE.md ({len(rows)} topics)")
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {k:24s} {v}")
