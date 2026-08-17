# Manuscript plan

Springer format. Written before drafting so the argument is fixed before the prose.

---

## 1. The objective, in one sentence

**Determine whether antibodies read protein surfaces by the same chemical rules that
govern co-evolved protein–protein interfaces — measured with a single pipeline on
both, and with the obvious confounds (what is bound, rather than how) explicitly
tested rather than assumed away.**

---

## 2. The literature gap

Antibody–antigen interfaces are among the best-characterised objects in structural
biology, and general protein–protein interfaces likewise. The gap is not that either
is unstudied; it is that **they have almost never been measured against each other on
one footing**, and that the comparisons which do exist leave the central confound
untested.

Four specific gaps this study addresses:

**(a) Cross-study comparison masquerading as evidence.** Statements like "antibody
interfaces are less hydrophobic" or "smaller" or "less shape-complementary" are
generally assembled from separate papers using different contact criteria, different
corpora, different H-bond definitions, and — most importantly — different or absent
redundancy control. Forty-five lysozyme complexes weigh as much as forty-five
independent systems in an unweighted count. We measure both classes with one
extractor, one criterion set, and cluster30 independence units throughout.

**(b) "What is bound" versus "how it is bound" is not separated.** Antibody paratopes
are CDR loops; antibody antigens are disproportionately β-rich viral surface
proteins. Any observed chemical difference could therefore be loop chemistry, or fold
composition, rather than antibody chemistry. This confound is stated in reviews but
rarely tested. We test it directly, and — critically — against **size-matched nulls**,
because restricting data to a matched stratum shrinks it and small samples attenuate
correlations on their own.

**(c) No residue-pair recognition grammar for protein–protein comparison.** The
pairwise "recognition code" framing is standard for protein–nucleic acid interfaces
(Luscombe-style amino-acid × base matrices) but has not been applied comparatively to
antibody versus general protein interfaces with proper independence weighting and
per-cell inference.

**(d) Readout-channel decomposition has never crossed over.** Protein–nucleic acid
work routinely separates *base* readout (identity-specific) from *backbone* readout
(shape-specific). The exact analogue exists for protein–protein interfaces —
side-chain versus main-chain mediated contact — and to our knowledge has not been
used to compare antibody and general interfaces. It turns out to be where the
difference lives.

**What we plan to achieve.** A matched, redundancy-controlled, confound-tested
answer to "are antibodies chemically special?", with the mechanism named at
functional-group resolution and validated on atomic-resolution structures — plus an
explicit account of what is *not* special, which is a large part of the answer.

---

## 3. The central claim, and its honest limits

**Claim.** Antibody and general protein–protein interfaces share their chemistry at
the level of chemical classes (r = 0.84) but use **formally distinct residue-level
grammars** (permutation test against a single-grammar null: observed r = 0.417 vs
null 0.798 [0.728, 0.850], p = 0.001). The difference is **not** explained by loop
context, fold class, resolution, redundancy, H-bond criterion, or antibody format.
It is concentrated in **backbone (shape) readout** — antibody and PPI side-chain
grammars correlate at r = 0.49 but their main-chain grammars at only r = 0.14 — and
the mechanism is **tyrosine phenol –OH hydrogen bonding to the antigen main chain**,
together with genuine cation–π, and notably *not* π-stacking.

**Limits we state up front.** The difference does not break transferability: a
general-PPI potential discriminates native from decoy antibody interfaces as well as
an antibody-specific one. Electrostatics are indistinguishable between the classes.
Side-chain readout *fraction* is indistinguishable. Side-chain rotamer restriction is
indistinguishable. The corpus is 100% X-ray and contains none of SAbDab's ~4,800
cryo-EM protein-antigen entries.

---

## 4. Figures, and the question each answers

### Main text

**Fig. 1 — Corpus, validation and reproducibility**
*Question: can these numbers be trusted?*
Panels: corpus composition and resolution matching (1.83 vs 1.80 Å median);
redundancy structure (largest cluster 8.8%, 149/233 singleton units); independent
re-extraction agreement (F1 0.9990 / 1.0000); Kabsch–Sander validation against PDB
HELIX/SHEET (0.889 over 622,662 residues).
*Sources:* `provenance_audit.json`, `corpus.json`, `structural_pass.json`.
Unusual to lead with, and deliberate: the comparison is only as good as the matching.

**Fig. 2 — Interface physics and thermodynamics**
*Question: what forces hold these interfaces together, and how do they differ?*
Panels: bond-type composition; the physicochemical panel (aromatic +0.53, apolar
−0.52, H-bond +0.44, cation–π +0.34, salt bridge n.s., charge pairing n.s.);
interface size and closest approach; **ΔCp and apolar/polar burial** (Addition 2).
*Answers:* antibodies substitute aromatic/H-bond chemistry for hydrophobic burial,
over a smaller and tighter interface — with the enthalpy/entropy inference converted
into a burial-thermodynamic measurement rather than left as composition-talk.
*Figures:* `phase2_physicochemistry.png`, `phase2_interface_inventory.png`,
`add2_thermodynamics.png`.

**Fig. 3 — The 20×20 recognition grammar**
*Question: do the two classes prefer the same residue pairings?*
Both arms, cluster30 primary with UniProt sensitivity inset, effect size and evidence
as separate visual channels; 6×6 chemical-class collapse.
*Answers:* shared at class level (r = 0.841), divergent at residue level (r = 0.418;
0.586 with Cys masked).
*Figures:* `heatmap_20x20_cluster30.png`, `heatmap_chemical_classes.png`.

**Fig. 4 — Conservation: what is universal, what is specific**
*Question: which pairings are shared, and is the sharing about attraction or
avoidance?*
Cell-by-cell agreement scatter with CI; conservation by bond type; universal /
arm-specific / opposed classification.
*Answers:* 55 universal cells, **mostly depletions** — the shared grammar is largely
a grammar of avoidance (like-charge, charge–apolar mismatch). Salt-bridge map
correlation is reported as **not estimable** (12 cells, CI spans zero).
*Figure:* `phase4_conservation.png`.

**Fig. 5 — Readout channels and functional-group mechanism** ← *the novel claim*
*Question: is the difference in what they read, or how they read it?*
Panels: side-chain vs main-chain matrices per arm; per-residue side-chain fraction
(r = 0.981 between arms — a clean null); cross-arm correlation by channel (side 0.49
vs main 0.14); functional-group 14×14 difference map; **PyMOL mechanism montage**
showing Tyr OH → antigen backbone at 1.10–1.30 Å (2.64–2.94 Å bonds).
*Answers:* identical channel *proportions*, divergent backbone *grammar*, with the
responsible chemistry named and structurally verified. Ring↔ring stacking is
*depleted* — aromatic enrichment is aromatics touching non-aromatics.
*Figures:* `phase5_readout.png`, `phase6b_functional_groups.png`,
`pymol_mechanism_montage.png`.

**Fig. 6 — Confound battery**
*Question: is any of this just loops, folds, resolution, or redundancy?*
Panels: secondary-structure composition per side; loop-confound correlations with
size-matched nulls; CATH fold composition and fold-matched control; resolution
strata; drop-top-N and leave-one-cluster-out; H-bond cutoff; Fab vs nanobody against
the split-half reliability ceiling.
*Answers:* every confound tested; none raises the correlation. Matching on context
*lowers* it (loop–loop r = 0.185, below its size-matched null). Epitopes are
structurally ordinary (loop fraction 0.512 vs 0.509) — the loop enrichment is
entirely paratope-side.
*Figures:* `phase6c_loop_confound.png`, `phase6d_folds.png`, `phase11_controls.png`.

**Fig. 7 — Paratope architecture and epitope character**
*Question: how does the antibody divide the labour, and what does it read?*
Panels: CDR vs framework by construct type (Fab H3 = 37%, nanobody H3 = 66%,
framework 21–22% throughout); paratope composition **against two nulls**; epitope
contiguity (median 4 segments, 97.5% discontinuous, 22 residues); **shape
complementarity** (Addition 3); **immunodominance map** — 42 independent antibodies
on lysozyme, Gly102 hit by 98%, Gini 0.49.
*Answers:* recovers the known architecture as a positive control, and shows that
Ser/Thr/Gly "paratope enrichment" **inverts sign** against the antibody's own
surface — a compositional artefact invisible to single-background studies.
*Figures:* `phase7_paratope.png`, `add3_shape_complementarity.png`,
`immunodominance_profile.png`.

**Fig. 8 — Formal comparison and the limits of the difference**
*Question: is the difference statistically real, and does it matter practically?*
Panels: differential map with CI-on-difference and BH; permutation test against the
single-grammar null; matrix asymmetry (antibody 0.757 vs PPI 0.028); dynamic-range
test; **transfer test** (Addition 1).
*Answers:* the grammars differ formally (p = 0.001) and the antibody grammar is
genuinely sharper at equal sample size — **but** a general-PPI potential transfers
successfully to antibody interfaces (AUC 0.833 vs the antibody potential's 0.779).
The difference is real and mechanistically specific yet not large enough to break a
general scoring function. Reported as a limit on the claim, not buried.
*Figures:* `phase8_differential.png`, `add1_transfer.png`.

**Table 1** — corpus. **Table 2** — literature concordance (16 topics, verdicts).

### Supplementary
Bond-decomposition panels; UniProt sensitivity maps; Ramachandran/χ₁/burial
distributions; PTM and glycan double dissociation (`phase9_modified.png`); per-cell
CSVs; the full audit report (39 checks across two audit scripts); every sensitivity
re-run.

---

## 5. Section-by-section: what each part must do

**Introduction.** Establish that both interface classes are well studied *separately*,
then make the gap concrete with the four points in §2. End on the question, not on a
claim.

**Methods.** Corpus and curation ledger; contact criteria; **independence units and
why weighting is not optional**; the statistical conventions (effect size from
weighted counts, evidence from raw counts, never mixed); the bootstrap and
permutation designs; **an explicit subsection on size-matched nulls and reliability
ceilings**, because three separate results in this study would have been misread
without them and that is a transferable methodological point.

**Results.** In figure order. Positive controls stated as positive controls before
any novel claim is made.

**Discussion.** What the grammar difference means mechanistically (loop scaffold +
disposable binder → enthalpic, aromatic, backbone-reading strategy); why transfer
nonetheless survives; what this implies for antibody design tools; and the negative
results as results.

**Limitations.** 100% X-ray, no cryo-EM half of the field; PPI structural analyses on
a 1,000/2,292 sample; Sc is a contact-surface approximation, not Connolly; CATH
coverage asymmetric (43% antibody vs 80% PPI); no PB electrostatics or MD; the
transfer test's potentials differ in contact count even when unit-matched.

**Data and code availability.** Every figure regenerable from `code/`; all
intermediate JSON/CSV included.

---

## 6. What must still be finished before drafting

1. **ΔCp** (running) — Fig. 2 depends on it.
2. **Sc** (running) — Fig. 7 depends on it.
3. **Two case studies** (lysozyme, GFP) — matched antigen read by antibody and by a
   natural partner; strengthens Fig. 5 or Fig. 7.
4. **PyMOL immunodominance surface render** — Fig. 7 panel.
5. Optional but valuable: hot spots (Bogan–Thorn) for Fig. 2; bound-vs-unbound CDR
   change from SAbDab's 1,150 apo structures for Fig. 7.

## 7. Journal positioning

With the above: a rigorous comparative-anatomy paper. Realistic targets *Proteins*,
*mAbs*, *BMC Bioinformatics*, *Bioinformatics Advances*. The transfer test coming
back negative removed the strongest "so what", so a higher tier would need either
experimental affinity validation (SKEMPI/AB-Bind) or the apo/bound conformational
analysis to produce something striking. That assessment is deliberately conservative.
