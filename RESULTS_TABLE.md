# Study register — every analysis, result, observation and figure

Two groups throughout: **antibody–antigen** (512 structures, 233 cluster30 units,
36,810 residue contacts) and **general protein–protein** (2,292 structures, 1,326
units, 406,397 contacts after deduplication). Every statistic is unit-weighted;
effect sizes come from weighted counts and evidence from raw counts, never mixed.

`effect` = rank-biserial (0 = distributions overlap completely, ±1 = fully separated).

---

## Foundation

| # | Analysis | Result | Observation | Figure / file |
|---|---|---|---|---|
| 0 | **Provenance audit** — re-extract contacts independently from raw coordinates | antibody F1 **0.9990** (2,048 pairs), PPI F1 **1.0000** (2,563); median distance deviation 0.0025 Å; 0 shared PDBs between arms | Inherited tables are re-derivable and uncontaminated; no deletion warranted | `provenance_audit.json` |
| 1 | **Corpus manifest** — 3,084 structures catalogued | antibody 100% X-ray, median **1.83 Å** vs PPI **1.80 Å**; 203/512 nanobodies; largest cluster 8.8%, 149/233 singleton units | Resolution confound absent by construction; but corpus contains **none** of SAbDab's 4,786 cryo-EM entries | Fig. 1 · `corpus_manifest.csv` |
| — | **Machinery audits** (2 scripts, 39 checks) | 39/39 pass; BH matches textbook exactly; log2-odds exactly 0 on rank-1 table; LOO subtraction exact to 3e−15 | Found and fixed **4 bugs**: PPI double-counting, dihedral convention, ASN-as-PTM, ΔCp sign | `audit_phases.json`, `audit_full.json` |

## Chemistry and physics

| # | Analysis | Result | Observation | Figure |
|---|---|---|---|---|
| 2 | **Interface physics** | aromatic **+0.53**, apolar–apolar **−0.52**, H-bond **+0.44**, cation–π **+0.34**, closest approach −0.41; salt bridge **n.s. (p=0.30)**, charge pairing **n.s. (p=0.35)** | Antibodies substitute aromatic/H-bond chemistry for hydrophobic burial. **Electrostatics do not distinguish them** | Fig. 2 · `phase2_physicochemistry.png` |
| A2 | **Burial thermodynamics (ΔCp)** | apolar fraction buried 0.50 vs 0.58 (eff −0.52); ΔCp/1000 Å² **−93.6 vs −151.3** cal/mol/K (eff +0.52, p=6e−29) | PPI binding is **1.6× more hydrophobically driven per unit area**. Confirms the enthalpy-leaning inference as measurement | Fig. 2 · `add2_thermodynamics.png` |
| 6b | **Functional groups** (14×14, 1.95M atom contacts) | map r = **0.618**; antibody-enriched: aromatic ring↔ammonium **+1.07**, Tyr phenol–OH↔backbone N–H **+0.56** / C=O **+0.58**; ring↔ring **depleted** | Names the mechanism: **cation–π and tyrosine-reads-backbone**, explicitly *not* π-stacking | Fig. 5 · `phase6b_functional_groups.png` |
| 9 | **PTMs and glycans** | glycans at interface antibody **14.1%** vs PPI 5.1% (OR 3.07, p=1.6e−11); PTMs antibody 1.6% vs PPI **9.0%** (OR 0.16, p=7.3e−11) | Clean **double dissociation**: PPI interfaces are PTM-regulated, antibody interfaces contend with glycan shield | Fig. 2 panel · `phase9_modified.png` |

## Recognition grammar

| # | Analysis | Result | Observation | Figure |
|---|---|---|---|---|
| 3 | **The 20×20 matrices** | grammar r = **0.418** (cluster30), 0.395 (UniProt sensitivity); chemical-class r = **0.841** | **Shared chemistry at class level, divergent grammar at residue level** | Fig. 3 · `heatmap_20x20_cluster30.png` |
| 4 | **Conservation** | Cys-masked r rises to **0.586**; 55 universal cells, **majority are depletions**; salt-bridge map **not estimable** (12 cells, CI spans 0) | The shared grammar is largely a **grammar of avoidance** — like-charge and charge–apolar mismatch | Fig. 4 · `phase4_conservation.png` |
| 5 | **Readout channels** | side-chain fraction **0.716 vs 0.716** (n.s.), per-residue r = **0.981**; cross-arm side-chain r = 0.49 vs **main-chain r = 0.135**; within-arm side/main r 0.28 (Ab) vs 0.71 (PPI) | **The divergence is in backbone (shape) readout, not sequence readout.** Antibody channels are decoupled; PPI channels redundant | Fig. 5 · `phase5_readout.png` |
| 8 | **Differential + inference** | permutation test: observed r **0.417** vs null **0.798** [0.728, 0.850], **p = 0.001**; 33/400 cells differ; asymmetry 0.757 vs 0.028; dynamic range genuinely sharper at equal n | Grammars are **formally distinct**, and the antibody matrix is genuinely directional | Fig. 8 · `phase8_differential.png` |
| A1 | **Transfer test** | antibody interfaces: own potential AUC 0.779, **PPI potential 0.833** (0.810 at matched units) | Shared chemistry is sufficient for transfer — **a quantitative bound on how large the grammar difference is** | Fig. 8 · `add1_transfer.png` |

## Confound battery — every one tested, none explains the difference

| # | Confound | Test | Result |
|---|---|---|---|
| 6c | **Loop context** | grammar r within matched secondary structure, vs size-matched null | loop–loop r **0.185** vs null [0.247, 0.390] → **below null; matching context lowers r** |
| 6d | **Fold class** | CATH-matched strata vs size-matched null | mainly-α 0.167, β 0.260, α/β 0.172 — **none rises**. (CATH coverage 43% Ab vs 80% PPI) |
| 11 | **Resolution** | high (<2 Å) vs low strata | 0.356 / 0.230, neither raises r |
| 11 | **Redundancy** | drop top 1/2/3/5 clusters; leave-one-cluster-out ×20 | r 0.402–0.457; LOO range 0.382–0.457 — **no single antigen drives it** |
| 11 | **H-bond criterion** | re-derive at 3.0 / 3.2 / 3.5 Å | cross-arm r 0.300 / 0.317 / 0.276 — stable |
| 11 | **Construct type** | Fab vs nanobody, against split-half reliability ceiling | observed 0.179 inside ceiling **0.253 [0.117, 0.370]** → formats share one grammar |

## Architecture: paratope and epitope

| # | Analysis | Result | Observation | Figure |
|---|---|---|---|---|
| 6a | **Structural pass** | Kabsch–Sander implemented from scratch, **88.9%** agreement with PDB HELIX/SHEET over 622,662 residues; 1,512 structures, 0 failures | No DSSP binary existed; validated rather than trusted | `phase6c_conformation.png` |
| 6c | **Secondary structure** | paratope **67% loop** vs PPI probe 46%; epitope **51.2% loop** vs PPI target 50.9% | **Epitopes are structurally ordinary.** The loop enrichment is entirely paratope-side | Fig. 6 · `phase6c_loop_confound.png` |
| 6c | **Torsions / rotamers** | interface vs background rotamer shift ~1–2 points in **both** arms | No difference in side-chain conformational restriction — entropy cost is not the discriminator | Supp · `phase6c_conformation.png` |
| 7 | **Paratope architecture** | IMGT mapping coverage **1.000**; Fab H3 = **37.2%** of CDR contacts, nanobody H3 = **65.7%**; framework 21–22% across all formats | Recovers H3 dominance as a positive control; framework contribution strikingly constant | Fig. 7 · `phase7_paratope.png` |
| 7 | **Paratope composition, two nulls** | Tyr +1.71/+2.70, Trp +1.14/+3.04; **Ser +1.14 / −0.69**, Thr +0.21/−0.94, Gly +0.44/−0.45 | **Ser/Thr/Gly enrichment inverts** against the antibody's own surface — a compositional artefact invisible to single-background studies | Fig. 7 |
| 7 | **Epitope character** | median **4 segments**, **97.5% discontinuous**, 22 residues | Textbook conformational epitopes — positive control | Fig. 7 |
| A5 | **Interface vs its own surface** | mean \|log₂ enrichment\|: paratope **0.969**, PPI interface 0.492, **epitope 0.258** | **Epitopes are barely distinguishable from their own antigen's surface** — antibodies impose chemistry on an arbitrary patch | Fig. 7 · `add5_epitope_logic.png` |
| A5 | **Protrusion** | paratope **+0.293**, PPI −0.163, **epitope −0.036** | The classic "epitopes protrude" claim **does not survive** against a same-surface control. It is the *paratope* that protrudes | Fig. 7 |
| A3 | **Shape complementarity** | antibody 0.464 vs PPI 0.467, **p = 0.082 (n.s.)** | Indistinguishable from the mixed PPI population (L&C compared against obligate complexes specifically) | Fig. 7 · `add3_shape_complementarity.png` |
| A4 | **Geometry** | paratope Rg **10.7 vs 16.3 Å** (eff −0.32); planarity **0.85 vs 0.79** (+0.25); epitope Rg 10.5 vs 16.0 | Antibody interfaces are **compact and flat**; PPI interfaces large and three-dimensional | Fig. 7 · `add4_geometry.png` |
| A4 | **PCA / separability** | grouped-CV classifier **AUC 0.875**; **93.8%** of antibody interfaces inside the PPI 95% ellipse (PC1+PC2 = 39% variance) | **Visually overlapping, statistically separable** — the signature is distributed across descriptors, not on one axis | Fig. 3b · `add4_pca.png` |

## Case study and validation

| # | Analysis | Result | Observation | Figure |
|---|---|---|---|---|
| — | **Immunodominance** (42 antibodies, one antigen) | **Gly102 contacted by 98%**, Asn103 83%, Asn19 83%, Trp63 76%; 66% of sequence ever contacted; **Gini 0.49** | Corpus redundancy turned into a result: strongly uneven targeting | Fig. 7 · `immunodominance_profile.png` |
| — | **Lysozyme case study** | catalytic Glu35 **12%**, Asp52 **19%** contacted, but substrate-binding Asp101 **69%**; immunodominant positions less conserved (0.47 vs 0.62, **p = 0.078, trend**) | Antibodies target the **rim** of the active-site cleft, not the catalytic residues; and prefer variable positions | `case_lysozyme.png` |
| — | **PyMOL mechanism** | 144 Tyr-OH→backbone bonds in 83 sub-1.8 Å structures; rendered at 1.10–1.30 Å, bonds 2.64–2.94 Å | Statistical claim **verified as atoms** | Fig. 5 · `pymol_mechanism_montage.png` |
| 10 | **Literature concordance** | 16 topics: 8 AGREE (positive controls), 2 PARTIALLY CONTRADICT, 1 NOVEL, 1 REFINE | Pipeline recovers known biology before departing from it | Table 2 · `LITERATURE_CONCORDANCE.md` |

---

## The negative results (reported, not buried)

| Finding | Numbers |
|---|---|
| Salt bridges / charge pairing identical | 0.033 vs 0.033, p = 0.30; 0.051 vs 0.051, p = 0.35 |
| Side-chain readout fraction identical | 0.716 vs 0.716, per-residue r = 0.981 |
| Rotamer restriction identical | ~1–2 point shift in both arms |
| Shape complementarity indistinguishable | 0.464 vs 0.467, p = 0.082 |
| General-PPI potential transfers successfully | AUC 0.833 vs 0.779 |
| Epitopes structurally ordinary | 51.2% vs 50.9% loop; protrusion −0.036 |
| Fab and nanobody grammars indistinguishable | 0.179, inside reliability ceiling |

## Methodological point worth its own section

**Three results would have been misread without size-matched nulls or reliability
ceilings** — the loop confound, the fold-matched control, and Fab-vs-nanobody.
Restricting data shrinks it, and small samples attenuate correlations on their own,
so a correlation that *drops* under stratification can mean either a real effect or
merely less data. Every stratified comparison here is therefore tested against a
null drawn at the same sample size. This is transferable and belongs in Methods.
