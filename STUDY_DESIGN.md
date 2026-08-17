# Study design — Antibody–antigen recognition vs general protein–protein interfaces

**Status: Phases 0–5 and 6b complete and audited; Phase 6a running.**
See §6 (Progress log) for what has actually been established, including three
corrections the audit forced.

**Directory:** `AI_v1/apta_tope/antibody_antigen_recognition_study/` (separate, self-contained)

---

## 0. The question

The nucleic-acid study asked whether aptamers read proteins by the same rules as natural
DNA/RNA, and answered: same chemistry, broader protein space, aptamer-specific matrix
needed. This study asks the protein-side counterpart:

> **Do antibodies read protein antigens by the same chemical rules that govern protein–protein
> interfaces generally — and if not, what is the antibody doing differently that makes it a
> universal binding tool?**

Two groups, deliberately parallel to the four-class design:

| Group | Definition | Structures | cluster30 units | Contacts |
|---|---|--:|--:|--:|
| **Group 1 — antibody** | Antibody–antigen, folded protein antigens (≥50 aa floor) | 512 | 233 | 36,810 |
| **Group 2 — PPI** | General non-antibody protein–protein complexes | 2,292 | 1,326 | 812,794 |

Group 2 is the control that tells us which antibody observations are "protein interfaces
being protein interfaces" and which are genuinely antibody-specific. **Every method applied
to Group 1 is applied identically to Group 2.**

### The asymmetry that structures the whole study

This is the single most important design decision, and it differs from the NA arm:

- **Antibody–antigen is directional and non-exchangeable.** One side is a germline-encoded,
  somatically matured, loop-dominated binder; the other is an arbitrary protein surface that
  never evolved to be bound. Rows (paratope) and columns (epitope) mean different things.
  The 20×20 matrix is **asymmetric** and must not be symmetrised.
- **General PPI is exchangeable.** Two co-evolved partners have no probe/target roles. The
  extractor emitted both directions; these are **symmetrised** explicitly.

Comparing an asymmetric map to a symmetric one is only legitimate if stated and handled —
Phase 8 does this with an explicit asymmetry test rather than by assumption.

---

## 1. What already exists (verified, not assumed)

Audited from `~/Documents/heatmap_validation/` — the working area the NA study was
extracted from. It was reported deleted in the NA `HANDOFF.md`; a copy survived. Data is
symlinked into this study, never copied or modified.

**Already done and defensible:**
- Corpus curation with a documented rejection ledger — 3,586 candidates → 774 usable
  (rejected: 1,761 no-antibody, 729 TCR/MHC, 322 no-antigen); detector concordance 0.934
  between an independent sequence rule and a description rule.
- A ≥50-residue antigen floor with sensitivity reported at 0/30/50/80 (774/530/512/481) —
  short peptide epitopes are excluded because they are not protein surfaces.
- Contact extraction, both arms, one row per residue pair: 4.5 Å contact, 4.0 Å salt bridge
  (Barlow–Thornton), 3.5 Å H-bond, 5.5 Å/35° stacking. Atom-level tables also present.
- cluster30/50/95 + UniProt independence units, resolution, antigen length/organism.

**Two integrity flags carried forward, not silently inherited:**
1. This project contains `_FABRICATED_removed/` and `_superseded_contaminated/` directories.
   Phase 0 re-verifies the antibody/PPI tables against the raw PDBs before anything is built
   on them.
2. `salt` counts are **not comparable to the NA arm** (basic N→carboxylate O here,
   basic N→phosphate O there). Same geometry, different acceptor chemistry.

**What does not exist:** any analysis of these two corpora. No 20×20, no heatmap, no DSSP,
no CATH/SCOP, no epitope/paratope logic. That is this study.

---

## 2. Tooling constraints (checked, not assumed)

| Need | Status | Design response |
|---|---|---|
| DSSP secondary structure | **`mkdssp`/`dssp`/`stride` all absent**; `Bio.PDB.DSSP` is only a wrapper around the missing binary | Phase 6a: try `pip install pydssp`; if unavailable, implement Kabsch–Sander directly (backbone H-bond energy, −0.5 kcal/mol cutoff, H/E/T/S assignment). It is ~150 lines of numpy and is validated against a held-out set of PDB `HELIX`/`SHEET` records — which is a *stronger* provenance statement than trusting a binary. |
| CDR assignment (paratope logic) | **ANARCI absent** | Use SAbDab's `POST /download/structures` — SAbDab serves **Chothia-numbered** coordinates, so CDRs are direct residue-number ranges. No numbering software needed. |
| CATH / SCOP fold assignment | `cathdb.info` REST 404s | **RCSB REST works** (`data.rcsb.org`, verified 200) and serves CATH + SCOP2 per polymer entity. Single source, cached locally. |
| SAbDab metadata | **Solved.** Old flat-file URL is dead (site is now a SPA); the live API is `opig.stats.ox.ac.uk/webapps/newsabdab/api` with an OpenAPI spec. Full summary already pulled: 21,914 rows, 16,535 protein-antigen. | Used in Phase 1 to enrich the existing 512 with isotype, species, CDR canonical class, and to define the expansion set. |

Available: numpy 2.4, scipy 1.18, sklearn 1.9, Biopython 1.87 (incl. ShrakeRupley SASA),
matplotlib, PyMOL CLI, python-docx/pptx. Absent: xgboost, AlphaFold/Boltz, any docking.

---

## 3. Chronological plan

Each phase is gated: it produces a named output, and you approve before the next starts.
Phases 1–5 are the core deliverable; 6–10 are the layers you asked to add; 11–12 close it out.

### Phase 0 — Provenance audit *(before anything is built on the data)*
Re-derive contacts for a random 25-structure sample straight from the raw PDBs with an
independently written extractor; require ≥99% pair agreement with the inherited CSVs. Confirm
the `_FABRICATED_removed`/`_superseded_contaminated` directories share no lineage with these
tables. Check the antibody/PPI corpora for overlap (an antibody complex sitting in Group 2
would blunt the comparison) and for the lysozyme/spike redundancy the unit weighting exists
to absorb.
→ `metadata/analysis/provenance_audit.json`. **If this fails, everything downstream stops.**

### Phase 1 — Corpus definition and enrichment
Freeze the two groups. Join SAbDab metadata onto Group 1: isotype, Fab/Fv/VHH/scFv type,
heavy/light subclass, species, antigen family. Decide — with the numbers in hand — whether to
expand Group 1 beyond 512 using the 8,462 protein-antigen PDB entries the API now exposes.
Record antigen-composition statistics (the NA study's target-composition reframing was its
most honest result; the same trap exists here — if antibody antigens are viral surface
proteins and PPI partners are enzymes, fold composition, not antibody chemistry, may drive
every difference).
→ `metadata/analysis/corpus.json`, corpus composition figure.

### Phase 2 — Interaction quantification
Per-arm, unit-weighted contact and bond inventories: contacts per interface, buried area
(ShrakeRupley), bond-type composition (salt/H-bond/vdW/stacking), interface size distributions.
Establishes the denominators every later normalisation uses.
→ `figures/interface_inventory.png`, `metadata/analysis/inventory.json`.

### Phase 3 — **The 20×20 recognition matrices** *(core deliverable)*
For each arm: log2(observed/expected) against the marginal-product null, unit-weighted;
2,000-resample bootstrap **over cluster30 units** (not contacts — within-structure contacts
are massively correlated); Fisher exact + Benjamini–Hochberg over all 400 cells; per-unit sign
replication. cluster30 primary, UniProt as sensitivity. Effect size and evidence shown as
separate visual channels, never merged into one asterisk.
Also the 6×6 chemical-class collapse, which is where the shared-chemistry question is actually
powered — 400 cells over 233 units is thin, and saying so is part of the result.
*(Script written: `code/aa_heatmap.py`, `code/aa_common.py`. Not yet run.)*
→ `figures/heatmap_20x20_{cluster30,uniprot}.png`, `heatmap_chemical_classes.png`,
`heatmap_bond_decomposition.png`, per-cell CSVs.

**Answers:** *what is the 20×20?*, *what are the basis of these interactions?*

### Phase 4 — Bond chemistry: shared or distinct?
Decompose each matrix by bond type and test the shared-chemistry claim directly: correlate
the antibody and PPI maps overall, within bond type, and within chemical class. Identify
"universal" pairs (same sign, significant, both arms) versus arm-specific ones. This mirrors
the NA arm's bonding-conservation study (which found 18 universal bonds, r = 0.87–0.99), so
the two studies' answers are directly comparable.
→ `figures/bond_conservation.png`.

**Answers:** *do they follow a shared chemistry?*, *is there a shared observation across
classes of amino acids?*

### Phase 5 — Does the side partner matter?
The contact tables already carry `hbond_ag_side/main`, `vdw_ag_side/main`,
`ag_side_contacts/ag_main_contacts`. Build **side-chain-mediated vs main-chain-mediated**
matrices separately. Side-chain reading is identity-specific (sequence readout); main-chain
reading is identity-agnostic (shape readout) — the direct protein analogue of the NA arm's
base-vs-backbone readout decomposition. Prediction to be tested, not assumed: antibodies
should skew to side-chain/identity reading of the epitope, since they must discriminate
between similar surfaces.
→ `figures/readout_side_vs_main.png`.

**Answers:** *does side partner matter?*

### Phase 6 — Structural context: DSSP + CATH/SCOP
- **6a DSSP.** Per-residue secondary structure for every interface residue, both arms
  (pydssp or the Kabsch–Sander implementation above, validated against PDB HELIX/SHEET).
  Recompute the 20×20 stratified by secondary structure. Central expected contrast:
  antibody paratopes are **loop**-dominated (CDRs), while general PPI interfaces are far more
  helix/strand — and if so, the "antibody grammar" may partly be **loop** grammar. That
  confound must be tested, not narrated.
- **6b CATH/SCOP** via RCSB REST. Fold-class composition of both groups' antigens/partners;
  re-run the key comparisons **fold-matched** to separate antibody chemistry from
  fold-composition differences. This is the direct analogue of the NA study's
  target-composition control, which is the result that reframed that paper.
→ `figures/dssp_context.png`, `figures/fold_composition.png`.

**Answers:** *how do they differ across [groups]?* — with the confounds actually controlled.

### Phase 7 — Epitope/paratope logic *(antibody-specific)*
Using Chothia-numbered SAbDab coordinates: CDR vs framework contribution per CDR loop
(H1/H2/H3, L1/L2/L3); the H3 dominance test; paratope composition (the classic Tyr/Trp/Ser
enrichment) tested against **both** a general-PPI null and an antibody-surface null; epitope
character (contiguous vs discontinuous, SASA, protrusion, secondary-structure composition);
epitope↔paratope chemical complementarity. Then the payoff question: **what about this logic
makes antibodies a general-purpose binding tool?** — framed as a hypothesis with a
quantitative test (loop-scaffold plasticity + Tyr/Ser-dominated low-cost chemistry allowing
one framework to read many surfaces), not as a claim.
→ `figures/paratope_logic.png`, `figures/epitope_character.png`.

### Phase 8 — Group 1 vs Group 2: differential and inferential statistics
The formal comparison. Cell-wise differential map (antibody − PPI) with bootstrap CIs on the
*difference*; permutation test for global map divergence; Cramér's V effect sizes; the
asymmetry test (is the antibody matrix significantly non-symmetric while PPI is not?);
fold-matched and secondary-structure-matched re-runs from Phase 6. Study-wide FDR accounting
and power analysis — with 400 cells and 233 units, honest power reporting is mandatory.
→ `figures/differential_antibody_vs_ppi.png`, `figures/statistics.png`.

### Phase 9 — Modified amino acids
**Finding already established during design:** the inherited contact tables contain **no**
modified residues — only the 20 standard, plus MSE (selenomethionine, a crystallographic
substitution) in the PPI arm. The extractor filtered them out. So this question **cannot** be
answered from the existing tables and requires going back to the raw PDBs to recover PTMs at
interfaces (SEP/TPO/PTR phospho-, MLY/M3L methyl-lysine, ALY acetyl-lysine, CSO/CSD oxidised
Cys, HYP hydroxyproline, glycosylation sites). Their interface frequency will likely be low;
the honest deliverable may be a well-powered null with a stated detection limit rather than a
grammar. Design commits to reporting whichever it turns out to be.
→ `figures/modified_residues.png`, `metadata/analysis/modified.json`.

### Phase 10 — Literature concordance
Every headline feature linked to prior work, with our statistic beside the published one and
disagreements stated as disagreements:
- Paratope Tyr/Trp/Ser enrichment and H3 dominance (Mian 1991; MacCallum/Martin/Thornton 1996;
  Ramaraj 2012; Peng 2014).
- Antigenicity/epitope propensity and discontinuous epitopes (Kringelum 2013; Rubinstein 2009).
- General PPI interface composition, hot spots, and the O-ring hypothesis
  (Jones & Thornton 1996; Bogan & Thorn 1998; Chakrabarti & Janin 2002; Levy 2010).
- Interface amino-acid propensity/statistical potentials (Glaser 2001; Moal & Fernández-Recio).
- Barlow & Thornton 1983 for the salt-bridge criterion already in use.
Each becomes a row in a concordance table: published claim → our measurement → agree/differ →
if differ, why (corpus, weighting, criterion). Same structure as the NA arm's external
validation, which recovered the Luscombe recognition code 8/10.
→ `LITERATURE_CONCORDANCE.md`, `figures/external_validation.png`.

### Phase 11 — Controls and robustness
Biological-assembly control (crystal packing); resolution deconfound (Group 2 is heavily
cryo-EM, Group 1 mixed — this **will** matter and is planned for, not discovered late);
H-bond criterion robustness; drop-top-N antigen / leave-one-system-out (lysozyme, spike, HIV
Env dominate the antibody corpus); antibody-type stratification (Fab vs VHH — nanobodies have
no light chain and a longer H3, so they are a natural internal contrast).
→ `figures/controls.png`.

### Phase 12 — Deliverables
`FULL_REPORT.md`, `MANUSCRIPT.md`, a `HANDOFF.md` matching the NA study's standard, PyMOL
interface galleries, and — if wanted — `.docx`/`.pptx` via the NA study's builders.

---

## 4. Deliberate design choices, and what would falsify us

- **Unit-weighted everything.** Raw contact counts would let 45 lysozyme complexes shout down
  233 independent clusters. Reported both ways where it changes an answer.
- **Fisher on raw counts, effect size on weighted counts.** Feeding fractional weights to
  Fisher invents significance. Evidence and effect are separate channels throughout.
- **The null is the marginal product**, so the map measures *preference*, not abundance —
  Tyr is common at interfaces; the question is whether Tyr–Tyr is commoner than Tyr's own
  commonness predicts.
- **The comparison can fail honestly.** If the antibody and PPI maps correlate at r ≈ 0.9,
  the answer is "antibodies use ordinary interface chemistry, and their power is in loop
  scaffold plasticity, not novel chemistry" — a real, publishable, and quite likely result.
  If the fold-matched and SS-matched controls in Phases 6/8 erase the differences, we report
  the reframing the way the NA study reported its target-composition reframing.

---

## 4b. Additional protein–protein analyses available to us

Everything below is feasible with the verified toolchain (numpy/scipy/sklearn/Biopython/
PyMOL) and the data already on disk — **no docking, no AlphaFold, no MD**, which are all
absent. Ordered by scientific payoff per unit of work. These extend Phases 3–8; each is
applied to **both groups**, so every one is a comparison, not a description.

**Tier A — high payoff, data already in hand**

- **A1 Hot-spot / O-ring analysis.** Per-residue contact burial and centrality; test the
  Bogan–Thorn hot-spot enrichment (Trp/Arg/Tyr) and the Levy core–rim–support architecture.
  Do antibody paratopes have hot spots in the same places, and are CDR-H3 residues
  preferentially core? Directly links to the most-cited PPI literature.
- **A2 Interface geometry.** Buried SASA, planarity, circularity, gap volume/gap index,
  packing density (Jones & Thornton's descriptor set). The classic discriminator between
  obligate, transient, and antibody interfaces — and a direct literature anchor.
- **A3 Interface conservation.** Per-position conservation from the cluster30/50/95 families
  already computed. Prediction to test: PPI interfaces are conserved (co-evolved), antibody
  **epitopes** are not conserved in the antigen (nothing evolved them to be bound), while
  antibody **paratopes** are hypervariable by design. A three-way conservation contrast is
  a strong, cheap result.
- **A4 Contact-order / sequence-separation structure.** Discontinuous vs contiguous epitopes
  quantified properly (how many sequence segments, how far apart) — the standard epitope
  descriptor, and a real difference from typical PPI interfaces.
- **A5 Pairwise-preference statistical potential + native-vs-decoy scoring.** Turn each 20×20
  into a scoring function; score native interfaces against rigid-body decoys, leave-one-unit-out.
  The exact payoff structure the NA study used (aptamer AUC 0.64; generic map scored aptamers
  *below chance*). The analogous headline here: does the general-PPI potential score antibody
  interfaces below chance? If yes, that is the quantitative case for antibody-specific design
  matrices.
- **A6 Charge/polarity complementarity and electrostatic pairing.** Net charge of each side,
  charge complementarity index, salt-bridge network topology. Antibody paratopes are unusually
  charge-neutral and aromatic-rich — testable against the PPI baseline.

**Tier B — strong, moderate work**

- **B1 Interface residue-propensity vs surface background.** Propensities computed against
  each protein's *own* solvent-accessible surface rather than the interface marginals — the
  Jones & Thornton normalisation. Different null, different question; worth having both.
- **B2 Water-mediated contacts.** Crystal waters are present in the X-ray subset. Bridging
  waters are a well-documented and under-analysed part of antibody recognition. Restricted to
  high-resolution X-ray only (cryo-EM has no ordered waters) — a clean, honest sub-analysis.
- **B3 Rotamer / side-chain entropy cost.** Side-chain conformational restriction at the
  interface, as a proxy for entropic binding cost. Speaks to why Tyr/Ser dominate paratopes:
  low entropic cost, high interaction density.
- **B4 Aromatic architecture.** π–π stacking, cation–π, T-shaped vs parallel geometry (the
  `stack`, `ring_dist`, `ring_angle` columns already exist and are currently unused). Aromatic
  chemistry is the most distinctive feature of paratopes; this is nearly free.
- **B5 Interface clustering / binding-mode typology.** Unsupervised clustering of interfaces by
  their composition vectors — do antibody interfaces form their own cluster, or fall inside the
  PPI distribution? A single figure that answers "are antibodies special?" visually.
- **B6 ML interface predictor with cross-group transfer.** sklearn gradient boosting on
  per-residue features, cluster-aware CV, ablation, and — the informative part — train on PPI,
  test on antibody and vice versa. Transfer failure is the strongest evidence for distinct
  grammars. Mirrors `strata_ml_predictor.py`.
- **B7 Interface-size/affinity scaling.** Where SAbDab or SKEMPI affinities can be joined,
  test buried-area↔affinity scaling per group. Requires an external affinity join; flagged as
  conditional.

**Tier C — worthwhile if the study grows**

- **C1 Hydrophobic patch / desolvation analysis** (obligate PPI is hydrophobic-core-like;
  antibody interfaces are not).
- **C2 Alanine-scanning-in-silico via the potential** — predicted ΔΔG per residue, validated
  against SKEMPI where entries overlap. Conditional on that join.
- **C3 Antibody-germline analysis** — do V-gene families have characteristic recognition
  signatures? Uses SAbDab subclass fields already downloaded.
- **C4 Epitope-overlap / immunodominance mapping** on antigens solved with many different
  antibodies (spike, HIV Env, lysozyme). Turns the corpus redundancy that unit-weighting
  suppresses into its own result.
- **C5 Cross-arm bridge to the nucleic-acid study** — compare the protein-side of NA-binding
  interfaces to the protein-side here, so all four NA classes plus these two groups sit on one
  axis. This is what would unify the two studies into a single paper.

**Recommendation.** Fold **A1–A6** into the core run (they are cheap and they are what
reviewers ask for), take **B4** immediately since the columns already exist, and hold the rest
until Phases 3–8 report. I would not start Tier C before the core comparison is on paper.

## 5. Open decisions for you

1. **Expand Group 1?** 512 → up to ~8,462 protein-antigen entries are now reachable via the
   API. More power, but a fresh download + extraction (hours) and re-curation. My
   recommendation: run Phases 0–5 on the existing 512 first, then decide with results in hand.
2. **Add a Group 3?** TCR–pMHC was explicitly rejected during curation (729 entries). It is
   the natural third class — another immune receptor, loop-mediated, but germline-constrained.
   It would sharpen "what makes antibodies special." Out of scope unless you want it.
3. **Nanobodies as their own group** rather than a stratification within Group 1.
4. **Which extension tier** from §4b enters the core run. My recommendation: A1–A6 plus B4
   in the core; B1–B3, B5–B6 held until Phases 3–8 report; Tier C deferred.

---

---

## 6. Progress log — what is actually established

### Completed

**Phase 0 — provenance audit: PASS.** Contacts re-derived independently from raw
coordinates for a 40-structure sample: antibody F1 = 0.9990 (2,048 pairs), PPI
F1 = 1.0000 (2,563 pairs), median distance deviation 0.0025 Å (CSV rounding). Zero
PDBs shared between arms. Tables resolve outside every quarantine directory.
**Conclusion: nothing fabricated; no deletion warranted.**

**Phase 1 — manifest.** `metadata/corpus_manifest.csv`, 3,084 structures. Antibody
arm is **100% X-ray**, median 1.83 Å vs PPI 1.80 Å — the resolution confound
predicted in §2 is *absent*, but the corpus contains **none** of SAbDab's 4,786
cryo-EM protein-antigen entries, so large viral glycoprotein antigens are
systematically missing. 203/512 are nanobodies (SD-H) — promoted to their own
stratum rather than pooled.

**Phase 2 — interface physics.** Antibody interfaces are *smaller* (59 vs 87 residue
pairs), *tighter* (3.587 vs 3.641 Å), far more **aromatic** (+0.53 effect), far less
**apolar–apolar** (−0.52), more **H-bonded** (+0.44), more **cation–π** (+0.34).
Salt bridges and opposite-charge pairing are **statistically identical** (p = 0.30,
0.35). Antibodies substitute aromatic/H-bond chemistry for hydrophobic burial —
an enthalpy-leaning strategy where general PPI is entropy-leaning.

**Phase 3 — the 20×20.** Grammar r(antibody, PPI) = **0.418** (cluster30), 0.395
(UniProt sensitivity — agrees, so not a clustering artefact). Chemical-class
r = **0.841**. Shared chemistry at class level, divergent grammar at residue level.

**Phase 4 — conservation.** Masking cysteine *raises* r to **0.586** (disulfides were
suppressing agreement). 55 universal cells, and the majority are **depletions** —
what interfaces agree on is what they refuse to do. Salt-bridge map correlation is
**not estimable** (12 cells, CI spans zero).

**Phase 5 — readout.** Side-chain-mediated fraction is ~0.72 in *both* arms, and
per-residue readout correlates at r = 0.981 — a clean null. But the **main-chain
grammar diverges** (cross-arm r = 0.135 vs 0.493 for side chains), and the antibody's
two channels are decoupled (r = 0.28) where PPI's are redundant (r = 0.71).

**Phase 6b — functional groups** (added at your request). 14-group chemistry from the
atom tables; map similarity r = 0.618. Antibody-enriched: **cation–π** (aromatic ring
↔ ammonium +1.07) and **Tyr phenol –OH ↔ antigen backbone** (+0.56 N–H, +0.58 C=O).
This names the mechanism behind Phase 5's divergent backbone grammar. Ring↔ring
stacking is *depleted* in antibodies, so aromatic enrichment is aromatics contacting
non-aromatics.

### Audit (`code/audit_phases.py`, 13/13 machinery checks pass)

Independent re-derivation confirmed BH matches textbook exactly, log2-odds is exactly
zero on a rank-1 table, the vectorised bootstrap matches explicit unit resampling, and
glycine reads exactly 0.000 side-chain. Three defects were found and **fixed**:

1. **PPI contacts were double-counted** (406,397 mirrored pairs, 100.0% correspondence).
   Correct for unit-weighted matrices, wrong for Fisher. Fixed via `dedup_rows()`;
   PPI significant cells 318 → **306**/400.
2. **"Contacts per probe residue" was an artefact** — retracted; now a clean null
   (+0.029, p = 0.50).
3. **Side-chain readout significance is estimator-dependent** (p = 0.037 vs 0.955).
   Effect is negligible either way; claim softened to "no meaningful difference".

### Running / next

- **Phase 6a** (running): Kabsch–Sander secondary structure implemented from scratch
  and **validated at 91.0% agreement** against PDB HELIX/SHEET records (1,659
  residues); plus φ/ψ, χ₁/χ₂ rotamers, and Shrake–Rupley burial. Antibody arm in
  full (512); PPI a seeded sample of 1,000/2,292 — reported, not implied.
- **Phase 6c** — torsion/dihedral analysis (your request): φ/ψ of read backbone,
  χ₁ rotamer restriction as an entropic-cost proxy.
- **Phase 7** — epitope/paratope logic, using SAbDab's **IMGT-numbered** coordinates
  (verified: served structures carry IMGT insertion codes 111A–112G, not Chothia
  27B/100x), so CDRs are fixed ranges with no ANARCI needed. Nanobodies must be
  analysed within construct type — 40% of the corpus has no light chain.
- **Case studies** (2, your request): **lysozyme C** (38 antibody vs 12 natural PPI
  complexes) and **GFP** (7 vs 6, includes nanobodies, no natural immune history).
  31 proteins are read by both an antibody and a natural partner — the direct
  analogue of the nucleic-acid study's DUX4 case.
- Phases 8–12 as designed, ending in a **Springer-format manuscript**.

### On SAbDab's census (established during Phase 1)

Inclusion is structural, not textual: chains are scanned for numerable Ig variable
domains with OPIG's ANARCII engine, then classified (FV/FAB/SD-H/VNAR; constructs
SCFV/DIABODY/…), with antigens auto-typed per entity. The unit is the **antibody
instance**, not the PDB entry — 21,914 instances over 11,458 entries; 8,462 entries
carry a protein antigen. Non-redundancy is a *service* (default 70% identity), not a
baked-in filter. Note the axis difference: SAbDab clusters by **antibody** sequence,
this study by **antigen** cluster30.
