# Findings against the literature — what we confirm, refine, contradict, and add

Every number below is computed by this study and reproducible from
`metadata/analysis/*.json`. Two groups throughout: antibody–antigen (512 structures,
233 cluster30 units) and general protein–protein (2,292 structures, 1,326 units).

> **Novelty caveat, stated up front.** The "NEW" classifications below reflect my
> assessment against known literature, not a systematic search — no literature
> database was queried in producing this study. **Every novelty claim must be
> verified against a proper search before submission.** The confirmations and
> contradictions are safer, because they are anchored to specific well-known claims;
> the novelty column is the one that can be wrong.

---

## A. Confirms the literature — the positive controls

These are not findings. They are the evidence that the pipeline works, and they must
be recovered before anything else in this study is trustworthy.

| # | Established claim | Our evidence | Source |
|---|---|---|---|
| A1 | Paratopes are Tyr/Trp-rich | Tyr **+1.71** vs PPI interfaces and **+2.70** vs the antibody's own surface; Trp **+1.14 / +3.04** — enriched on *both* nulls | Mian 1991; MacCallum 1996; Ramaraj 2012 |
| A2 | CDR-H3 dominates antigen contact | Fab H3 = **37.2%** of CDR contacts (next: H2 17.7%, L3 17.3%); L2 contributes 2.4% | Kabat & Wu; MacCallum 1996; Xu & Davis 2000 |
| A3 | Framework contributes a minority of contacts | **21–22%**, and strikingly constant across Fab / Fv / nanobody | MacCallum 1996; Sela-Culang 2013 |
| A4 | Epitopes are conformational/discontinuous | **97.5%** discontinuous, median **4 segments**, median **22 residues** | Barlow 1986; Rubinstein 2009; Kringelum 2013 |
| A5 | Charge complementarity is a conserved interface feature | K–D, D–K, K–E, R–D, R–E top-enriched and sign-stable in **both** arms | Barlow & Thornton 1983; Jones & Thornton 1996 |
| A6 | Antibody interfaces are smaller than typical PPI | 59 vs 87 residue pairs (effect −0.31, p=3e−13); 350 vs 468 atom contacts | Jones & Thornton 1996; Lo Conte 1999 |
| A7 | Antibody interfaces are less hydrophobic | apolar–apolar contact 0.145 vs 0.261 (effect **−0.52**, p=1e−33) | Jones & Thornton 1996; Chakrabarti & Janin 2002 |
| A8 | Nanobodies rely disproportionately on CDR-H3 | H3 = **65.7%** of nanobody CDR contacts vs 37.2% in Fabs | Muyldermans 2013; De Genst 2006 |
| A9 | Cation–π matters in antibody recognition | cation–π contact fraction 0.056 vs 0.032 (effect +0.34); aromatic ring↔ammonium **+1.07** enriched at functional-group level | Gallivan & Dougherty 1999 |
| A10 | Glycan shield obstructs antibody access | interfacial glycans in **14.1%** of antibody vs 5.1% of PPI structures (OR 3.07, p=1.6e−11) | Wei 2020 and the viral-glycoprotein literature |

## B. Refines the literature — agrees, but our data conditions it

| # | Claim as usually stated | What we measure | The refinement |
|---|---|---|---|
| B1 | Ser, Thr and Gly are enriched paratope residues | Ser **+1.14 vs PPI** but **−0.69 vs the antibody's own surface**; Thr +0.21 / −0.94; Gly +0.44 / −0.45 | The enrichment **inverts** against the right background. It reflects antibody *composition*, not selective paratope usage. Only Tyr, Trp (and weakly Asp, Asn) are enriched on both nulls. Invisible to any single-background study |
| B2 | Antibody interfaces are aromatic-rich | aromatic involvement 0.431 vs 0.306 (effect **+0.53**), but ring↔ring stacking **depleted** (W–Y Δ −0.99, Y–Y −0.85, both significant) and formal stacking is 0.3% of contacts | The aromatics contact **non-aromatics** — backbone and cations. Frequently blurred into "π-stacking" in the literature; it is not stacking |
| B3 | Antibody–antigen interfaces have lower shape complementarity | Sc 0.464 vs 0.467, **p = 0.082 (n.s.)** against general PPI | Lawrence & Colman compared against *obligate and protease–inhibitor* complexes specifically. Against a mixed PPI population the difference disappears. The original claim is not contradicted; its scope is narrower than often quoted |

## C. Contradicts the literature

| # | Established claim | Our evidence | Why the discrepancy |
|---|---|---|---|
| C1 | **Epitopes lie on protruding regions** | protrusion effect **−0.036** (epitope vs the antigen's own surface). The **paratope** protrudes strongly (**+0.293**) | The classic claim appears to be relative to the whole protein, not to other surface patches. Measured against the correct control it vanishes. It is the antibody that protrudes, not the antigen |
| C2 | **Epitopes are loop/turn-enriched** | epitope loop fraction **0.512** vs PPI partner surfaces **0.509** — identical | Same issue: loop enrichment is real versus whole-protein composition, absent versus other interfaces. The loop enrichment in antibody complexes is **entirely paratope-side** (0.670 vs 0.462) |
| C3 | **Antibody interfaces are ordinary interfaces, distinguished mainly by their loop scaffold and lack of co-evolution** | Chemistry is shared at class level (r = **0.841**) but the residue-level grammars are **formally distinct** (permutation p = **0.001**, observed r 0.417 vs null 0.798). And it is **not** the loop scaffold: matching both arms on secondary structure (loop–loop r = **0.185**) or CATH fold **lowers** the correlation, below size-matched nulls | The "ordinary interface" framing survives at the level of chemical classes but fails at residue-pair resolution, and the standard explanation for the difference is excluded |

## D. New findings — not established in the literature to my knowledge

**Requires verification against a systematic search before submission.**

| # | Finding | Evidence | Why it matters |
|---|---|---|---|
| **D1** | **The antibody/PPI difference lives in BACKBONE readout, not side-chain readout** | side-chain-mediated fraction identical (0.716 vs 0.716, per-residue r = **0.981**), but cross-arm grammar correlation is **0.49 for side-chain** vs **0.135 for main-chain** channels | The channel decomposition is standard in protein–nucleic-acid work (base vs backbone) and appears never to have been applied comparatively here. It relocates the entire difference |
| **D2** | **The mechanism: tyrosine phenol –OH hydrogen bonding to the antigen main chain** | most antibody-enriched functional-group contacts: Tyr –OH ↔ backbone N–H **+0.56**, ↔ C=O **+0.58**; verified structurally in **144 bonds across 83 sub-1.8 Å structures**, rendered at 1.10–1.30 Å with 2.64–2.94 Å geometry | Names a specific chemical mechanism for D1, with atomic-resolution confirmation |
| **D3** | **Antibody readout channels are decoupled; PPI channels are redundant** | within-arm side-vs-main map correlation **0.28** (antibody) vs **0.71** (PPI) | Antibodies extract independent information from sequence and shape; co-evolved interfaces read both similarly |
| **D4** | **The shared grammar is largely a grammar of AVOIDANCE** | of 55 universal cells (significant + sign-stable in both arms), the majority are **depletions** — like-charge and charge–apolar mismatch | What protein interfaces agree on is what they refuse to do, more than what they seek |
| **D5** | **Epitopes are chemically indistinguishable from their own antigen's surface** | mean \|log₂ enrichment\| vs own surface: paratope **0.969**, PPI interface 0.492, **epitope 0.258** | The clean structural statement of imposed vs co-evolved binding |
| **D6** | **Antibodies bind a different face of the same protein than its natural partner** | 16 proteins solved with both: median Jaccard **0.033**; **38%** share no contacted position; median 10.8% of epitope inside the partner site | Controls the target-identity confound by construction rather than statistically |
| **D7** | **Antibodies avoid functional sites; natural partners overlap them** | fraction of interface inside a ligand-binding site: **1.6%** (antibody) vs **10.7%** (PPI), effect −0.45, p=8e−09 | Corpus-wide confirmation of D5/D6 |
| **D8** | **Antibody binding is measurably less hydrophobically driven** | ΔCp per 1000 Å² **−93.6** vs **−151.3** cal/mol/K (effect +0.52, p=6e−29); apolar fraction buried 0.50 vs 0.58 | Converts the qualitative "less hydrophobic" claim into burial thermodynamics |
| **D9** | **Antibody interfaces are compact and flat; PPI interfaces large and 3-D** | paratope Rg **10.7 vs 16.3 Å** (effect −0.32); planarity **0.85 vs 0.79** (+0.25); same for epitope side | A geometric signature distinct from the size difference alone |
| **D10** | **PTM/glycan double dissociation** | interfacial PTMs: PPI **9.0%** vs antibody 1.6% (OR 0.16); glycans: antibody **14.1%** vs PPI 5.1% (OR 3.07); both p ≈ 1e−11 | PPI interfaces are chemically regulated; antibody interfaces contend with glycan |
| **D11** | **The antibody grammar is genuinely sharper at equal sample size** | SD of log₂ odds **0.572** (antibody) vs PPI subsampled to the same unit count **0.397 [0.364, 0.437]** | General PPI is a heterogeneous mixture of binding modes whose average washes out preferences; antibodies are one coherent mode |
| **D12** | **Fab and nanobody share one recognition grammar** | observed cross-format r **0.179**, inside the Fab split-half reliability ceiling **0.253 [0.117, 0.370]** | Different architecture (H3 66% vs 37%), same chemistry — so "antibody grammar" is a coherent object |
| **D13** | **A general-PPI potential transfers successfully to antibody interfaces** | AUC **0.833** (PPI potential) vs **0.779** (antibody's own, leave-one-unit-out); 0.810 at matched units | Places a quantitative **upper bound** on how large the grammar difference is. Constrains claims that antibody design needs bespoke matrices |
| **D14** | **Immunodominance is strongly focused and forms one contiguous face** | 42 independent antibodies on lysozyme: **Gly102 contacted by 98%**, Gini **0.49**, 66% of sequence ever contacted; surface render shows a single contiguous immunodominant face | Turns corpus redundancy into a result rather than a nuisance |
| **D15** | **Antibodies converge on the rim of the active site, not the catalytic residues** | Glu35 **12%**, Asp52 **19%**, but substrate-binding Asp101 **69%**; immunodominant positions less conserved (0.47 vs 0.62, p = 0.078, **trend only**) | Relevant to neutralisation mechanism; the conservation result is a trend, not significant |

## E. Methodological contribution

| # | Contribution | Why it is worth reporting |
|---|---|---|
| E1 | **Size-matched nulls for stratified comparisons** | Restricting data to a matched stratum shrinks it, and small samples attenuate correlations by themselves. **Three results here would have been misread without this** — the loop confound, the fold-matched control, and Fab-vs-nanobody |
| E2 | **Reliability ceilings for cross-group correlation** | Comparing correlations computed at different sample sizes compares precisions, not biology. Fab-vs-nanobody looked like a divergence (0.179 vs 0.405) until benchmarked against the Fab corpus's own split-half reliability |
| E3 | **Two-background enrichment** | Interface composition scored against both a foreign interface set *and* the same protein's own surface. This is what overturns B1 and produces D5 |
| E4 | **Kabsch–Sander implemented and validated rather than trusted** | No DSSP binary available; implemented from the 1983 definition and validated at **88.9%** against PDB HELIX/SHEET records over 622,662 residues |
| E5 | **Audit-driven analysis** | 39 checks across two audit scripts, each re-deriving a quantity independently. **Five bugs found and fixed**: PPI double-counting, a dihedral convention error, ASN-as-PTM, a ΔCp sign error, and a UniProt direction mis-mapping. Four produced plausible-looking wrong output |

---

## Summary

| Category | Count |
|---|--:|
| Confirms literature (positive controls) | 10 |
| Refines / conditions literature | 3 |
| Contradicts literature | 3 |
| New findings (pending verification) | 15 |
| Methodological contributions | 5 |

The ten confirmations are what license the fifteen new claims: a pipeline that failed
to recover Tyr/Trp paratope enrichment, H3 dominance, or discontinuous epitopes would
not be trustworthy on anything else.

The strongest new results are **D1/D2/D3** (the difference is in backbone readout,
mechanised by tyrosine hydroxyl) and the **D5–D7 cluster** (epitopes are arbitrary
accessible patches, not evolved or functional surfaces), which together give a
coherent account: *the antibody is the specialised object; the epitope is ordinary
surface.*

**D13 bounds the claim** and should be reported alongside it rather than omitted —
the difference is real, mechanistically specific, and not large enough to break a
general scoring function.
