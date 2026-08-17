# Literature concordance

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


## Paratope aromatic enrichment

**Published claim.** Tyr and Trp are strongly over-represented in antibody combining sites; Tyr is the single commonest paratope residue.  

*Sources:* Mian, Bradwell & Olson 1991; MacCallum, Martin & Thornton 1996; Ramaraj et al. 2012; Peng et al. 2014

**Our measurement.** Tyr log2 enrichment +1.71 vs general PPI interfaces and +2.70 vs the antibody's own surface; Trp +1.14 / +3.04

**Verdict: AGREE (positive control).** Recovered on both nulls, so this is selective use of aromatics, not a consequence of antibodies being aromatic-rich overall.


## Paratope Ser/Thr/Gly enrichment

**Published claim.** Ser, Thr and Gly are commonly listed alongside Tyr/Trp as enriched paratope residues.  

*Sources:* MacCallum et al. 1996; Collis et al. 2003; Ramaraj et al. 2012

**Our measurement.** Ser +1.14 vs PPI but -0.69 vs the antibody's own surface; Thr +0.21 / -0.94; Gly +0.44 / -0.45

**Verdict: REFINE.** Enriched relative to general protein interfaces but DEPLETED relative to the antibody's own accessible surface. The apparent enrichment is a property of antibody composition, not of selective paratope usage. Only Tyr and Trp are enriched on both nulls. This distinction is invisible to any study using a single background.


## CDR-H3 dominance

**Published claim.** CDR-H3 contributes more antigen contacts than any other loop and is the principal determinant of specificity.  

*Sources:* Kabat & Wu; MacCallum, Martin & Thornton 1996; Xu & Davis 2000

**Our measurement.** Fab: H3 = 37.2% of CDR contacts (next: H2 17.7%, L3 17.3%, L2 2.4%). Nanobody: H3 = 65.7%

**Verdict: AGREE (positive control).** L2 contributes almost nothing, also as reported. Nanobodies push H3 dominance far higher, consistent with the absence of a light chain and a longer H3.


## Framework contribution to the paratope

**Published claim.** Framework residues contribute a minority but non-negligible fraction of antigen contacts.  

*Sources:* MacCallum et al. 1996; Sela-Culang, Kunik & Ofran 2013

**Our measurement.** framework = 21.6% of paratope contacts, stable across Fab/Fv/nanobody (0.209-0.222)

**Verdict: AGREE.** Remarkably constant across construct types.


## Epitopes are discontinuous

**Published claim.** Most protein antigenic epitopes are conformational / discontinuous, comprising several sequence segments and roughly 15-25 residues.  

*Sources:* Barlow, Edwards & Thornton 1986; Rubinstein et al. 2009; Kringelum et al. 2013; Sela-Culang et al. 2013

**Our measurement.** median 4 sequence segments, 97.5% of epitopes discontinuous (>=2 segments), median 22 epitope residues

**Verdict: AGREE (positive control).** Both the segment count and the epitope size fall squarely in the published range.


## Charge complementarity at interfaces

**Published claim.** Salt bridges and complementary charge pairing are a conserved feature of protein-protein interfaces.  

*Sources:* Barlow & Thornton 1983; Jones & Thornton 1996; Xu, Tsai & Nussinov 1997

**Our measurement.** K-D, D-K, K-E, R-D, R-E are top-enriched and bootstrap-sign-stable in BOTH arms; salt-bridge fraction is statistically identical (antibody 0.033 vs PPI 0.033, p=0.30) as is opposite-charge pairing (0.051 vs 0.051, p=0.35)

**Verdict: AGREE, and EXTEND.** Charge chemistry is genuinely universal -- and therefore is NOT what distinguishes antibodies. A negative result worth stating.


## Interface size

**Published claim.** Typical protein-protein interfaces bury ~1,600 A2 total and involve ~20-30 residues per side; antibody-antigen interfaces sit at the smaller end of the transient range.  

*Sources:* Jones & Thornton 1996; Lo Conte, Chothia & Janin 1999; Chakrabarti & Janin 2002

**Our measurement.** antibody median 59 residue pairs vs PPI 86 (effect -0.31, p=2.8e-13); epitope median 22 residues

**Verdict: AGREE.** Antibody interfaces are significantly smaller than general PPI interfaces on every size measure.


## Hydrophobicity of antibody interfaces

**Published claim.** Antibody-antigen interfaces are less hydrophobic than obligate/homodimeric interfaces and more polar in character.  

*Sources:* Jones & Thornton 1996; Lo Conte et al. 1999; Chakrabarti & Janin 2002

**Our measurement.** apolar-apolar contact fraction antibody 0.145 vs PPI 0.261 (effect -0.52, p=1.4e-33); H-bond fraction 0.253 vs 0.200 (effect +0.44)

**Verdict: AGREE, and QUANTIFY.** The largest single physicochemical separations we measure. Antibodies substitute aromatic and H-bond chemistry for hydrophobic burial.


## Aromatic/pi chemistry

**Published claim.** Aromatic residues, and cation-pi interactions, are important in antibody recognition.  

*Sources:* Gallivan & Dougherty 1999; Ramaraj et al. 2012; Peng et al. 2014

**Our measurement.** aromatic involvement antibody 0.431 vs PPI 0.306 (effect +0.53); cation-pi 0.056 vs 0.032; functional-group map: aromatic ring <-> ammonium +1.07 antibody-enriched

**Verdict: AGREE, and REFINE.** But formal ring-ring stacking is DEPLETED in antibodies (W-Y delta -0.99, Y-Y -0.85, both significant). The aromatic enrichment is aromatics contacting NON-aromatics -- backbone and cations -- not pi-stacking. Frequently blurred in the literature.


## What antibodies read on the antigen

**Published claim.** Antibody recognition is usually described in terms of side-chain complementarity and shape matching.  

*Sources:* Lawrence & Colman 1993; Sundberg & Mariuzza 2002

**Our measurement.** side-chain-mediated fraction is IDENTICAL in both arms (~0.72, per-residue correlation r=0.981), but the main-chain readout GRAMMAR diverges (cross-arm r=0.135 vs 0.493 for side chains); Tyr phenol -OH to antigen backbone N-H/C=O is +0.56/+0.58 antibody-enriched

**Verdict: NOVEL / EXTENDS.** The distinguishing feature is not how MUCH backbone is read but HOW it is read. We are not aware of this decomposition being reported; it is the study's primary novel claim.


## Are antibody interfaces chemically special?

**Published claim.** Antibody-antigen interfaces are often described as 'ordinary' protein interfaces distinguished mainly by their loop scaffold and lack of co-evolution.  

*Sources:* Jones & Thornton 1996; Sundberg & Mariuzza 2002; Ofran et al. 2008

**Our measurement.** grammar correlation r=0.418 (0.586 with Cys masked; chemical-class r=0.841); permutation test of a single shared grammar: observed r=0.417 vs null 0.798 [0.728, 0.850], p=0.0010

**Verdict: PARTIALLY CONTRADICT.** Same chemistry at the class level (r=0.84) but formally distinct grammars at residue level (p=0.001). And the difference is NOT the loop scaffold: matching both arms on secondary structure (loop-loop r=0.185) and on CATH fold class does not raise the correlation -- it lowers it, below size-matched nulls.


## Epitope secondary structure

**Published claim.** Epitopes are often reported as enriched in loops/turns and in protruding, flexible regions.  

*Sources:* Thornton et al. 1986; Rubinstein et al. 2009; Kringelum et al. 2013

**Our measurement.** epitope loop fraction 0.512 vs PPI partner loop fraction 0.509; paratope loop fraction 0.670 vs PPI probe 0.462

**Verdict: PARTIALLY CONTRADICT.** Epitopes are NOT structurally unusual: their loop content is indistinguishable from ordinary PPI partner surfaces (0.512 vs 0.509). The loop enrichment is entirely on the PARATOPE side. Claims that epitopes are loop-rich appear to be relative to whole protein rather than to other interfaces.


## Nanobody paratopes

**Published claim.** Nanobodies (VHH) rely disproportionately on a long CDR-H3 and can access clefts inaccessible to Fabs.  

*Sources:* Muyldermans 2013; De Genst et al. 2006

**Our measurement.** H3 = 65.7% of nanobody CDR contacts vs 37.2% in Fabs; Fab-vs-nanobody grammar r=0.179, inside the Fab split-half reliability ceiling 0.253 [0.117, 0.370]

**Verdict: AGREE on architecture, EXTEND on chemistry.** H3 dominance is confirmed and much stronger in nanobodies. But the two formats' recognition CHEMISTRY is indistinguishable at matched sample size -- different architecture, same grammar.


## PTMs and glycans at interfaces

**Published claim.** Glycan shields obstruct antibody access to viral antigens; PTMs regulate protein-protein interactions.  

*Sources:* Wei et al. 2020 (glycan shield); Nishi, Hashimoto & Panchenko 2011 (PTMs at interfaces)

**Our measurement.** interfacial glycans: antibody 14.1% vs PPI 5.1% (OR 3.07, p=1.6e-11); interfacial PTMs: antibody 1.6% vs PPI 9.0% (OR 0.16, p=7.3e-11)

**Verdict: AGREE, and QUANTIFY a double dissociation.** Antibody interfaces engage glycans ~3x more; PPI interfaces engage PTMs ~6x more. Partly a target-composition effect (viral glycoproteins vs signalling complexes) and stated as such.


## Shape complementarity (Sc)

**Published claim.** Antibody-antigen interfaces have lower shape complementarity (~0.64-0.68) than protease-inhibitor or oligomeric interfaces (~0.70-0.76).  

*Sources:* Lawrence & Colman 1993

**Our measurement.** NOT YET COMPUTED

**Verdict: PENDING.** Required for publication: Lawrence & Colman defined Sc on exactly this comparison, and reviewers will expect it. Our smaller-but-tighter result (closest approach 3.587 vs 3.641 A) makes a specific prediction here, so it is a real test rather than a formality.


## Interface hot spots

**Published claim.** Binding free energy concentrates in a few 'hot spot' residues, enriched in Trp, Arg and Tyr, surrounded by an occluding O-ring.  

*Sources:* Clackson & Wells 1995; Bogan & Thorn 1998

**Our measurement.** NOT YET COMPUTED

**Verdict: PENDING.** Burial data from Phase 6a is in place; the analysis is cheap and connects to the most-cited PPI literature.


## Summary

- **AGREE (positive control)**: 3
- **AGREE**: 2
- **PARTIALLY CONTRADICT**: 2
- **PENDING**: 2
- **REFINE**: 1
- **AGREE, and EXTEND**: 1
- **AGREE, and QUANTIFY**: 1
- **AGREE, and REFINE**: 1
- **NOVEL / EXTENDS**: 1
- **AGREE on architecture, EXTEND on chemistry**: 1
- **AGREE, and QUANTIFY a double dissociation**: 1

Positive controls recovered: paratope Tyr/Trp enrichment, CDR-H3 dominance, framework contribution, discontinuous epitopes, charge complementarity, interface size, reduced hydrophobicity. The pipeline reproduces the established picture before departing from it.
