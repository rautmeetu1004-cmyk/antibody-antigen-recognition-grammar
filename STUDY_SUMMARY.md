# Study summary — what each phase asks, gives, and found

**Scope:** antibody–antigen (512 structures, 233 cluster30 units, 36,810 contacts)
vs general protein–protein (2,292 structures, 1,326 units, 406,397 contacts).
**Outputs:** 30 scripts · 33 figures · 28 result JSONs · 39 audit checks · 5 bugs
found and fixed.

---

## Foundation — is the data trustworthy?

| Phase | Asks | Gives | Headline result |
|---|---|---|---|
| **0** Provenance | Were these inherited tables fabricated? | Independent re-extraction from raw coordinates | antibody F1 **0.9990**, PPI **1.0000**, 0 shared PDBs → **clean, no deletion needed** |
| **1** Manifest | What exactly is in the corpus? | 3,084-structure catalogue with exclusion reasons | 100% X-ray, **1.83 vs 1.80 Å matched**; 203/512 nanobodies; **no cryo-EM at all** |
| **Audits** | Does the machinery compute what it claims? | 39 independent re-derivations | 39/39 pass; **5 bugs caught**: PPI double-count, dihedral convention, ASN-as-PTM, ΔCp sign, UniProt direction |

## Chemistry — what forces are at work?

| Phase | Asks | Gives | Headline result |
|---|---|---|---|
| **2** Physics | Which forces hold each interface together? | Bond/chemistry composition per interface | aromatic **+0.53**, apolar **−0.52**, H-bond **+0.44**; **salt bridges identical (p=0.30)** |
| **A2** Thermodynamics | Enthalpy- or entropy-driven? | ΔCp from per-atom buried surface | ΔCp/1000 Å² **−93.6 vs −151.3** → PPI is **1.6× more hydrophobically driven** |
| **6b** Functional groups | Which chemical groups actually touch? | 14×14 map from 1.95M atom contacts | **Tyr –OH → antigen backbone +0.56/+0.58**; cation–π +1.07; **ring↔ring depleted** |
| **9** Modified residues | Do PTMs and glycans participate? | Corpus-wide PTM/glycan scan | **Double dissociation**: glycans 14.1% vs 5.1%; PTMs 1.6% vs 9.0% (both p≈1e−11) |

## Grammar — do they use the same rules?

| Phase | Asks | Gives | Headline result |
|---|---|---|---|
| **3** The 20×20 | Same residue-pair preferences? | Two recognition matrices + class collapse | **class r = 0.841, residue r = 0.418** → shared chemistry, divergent grammar |
| **4** Conservation | What is universal? | Universal / arm-specific classification | 55 universal cells, **majority are depletions** → a **grammar of avoidance** |
| **5** Readout | In *what* they read, or *how*? | Side-chain vs main-chain decomposition | fraction identical (0.716/0.716) but **main-chain grammar r = 0.135 vs side-chain 0.49** |
| **8** Inference | Is the difference real? | Permutation test + differential map | observed **0.417** vs null **0.798** → **p = 0.001**; asymmetry 0.757 vs 0.028 |
| **A1** Transfer | Does it matter practically? | Native-vs-decoy scoring, cross-potential | PPI potential **0.833** > antibody's own 0.779 → **difference does not break transfer** |

## Confounds — is it an artefact? *(all six tested against size-matched nulls)*

| Phase | Confound | Verdict |
|---|---|---|
| **6c** | Loop context | loop–loop r **0.185**, *below* null → **matching context lowers r** |
| **6d** | CATH fold class | α 0.167 / β 0.260 / α-β 0.172 → **none rises** |
| **11** | Resolution | 0.356 (high) / 0.230 (low) → no rescue |
| **11** | Redundancy | drop-top-5 and leave-one-out: r 0.382–0.457 → **no single antigen drives it** |
| **11** | H-bond criterion | 3.0/3.2/3.5 Å → r 0.300/0.317/0.276 → stable |
| **11** | Construct type | Fab-vs-nanobody 0.179, inside ceiling [0.117, 0.370] → **one grammar** |

## Architecture — how is the binding built?

| Phase | Asks | Gives | Headline result |
|---|---|---|---|
| **6a** Structural pass | Local structure of every interface residue | SS + φ/ψ + χ₁/χ₂ + burial, 1,019,571 residues | Kabsch–Sander written from scratch, **validated 88.9%** vs PDB records |
| **6c** SS composition | Are paratopes and epitopes unusual? | Per-side secondary structure | paratope **67% loop** vs PPI 46%; **epitope 51.2% vs PPI 50.9% — identical** |
| **7** Paratope | How is labour divided? | IMGT CDR assignment (coverage **1.000**) | Fab H3 **37.2%**, nanobody H3 **65.7%**, framework 21–22% throughout |
| **7** Two-null composition | Is Tyr/Ser enrichment real? | Enrichment vs PPI *and* vs own surface | **Ser inverts** (+1.14 → −0.69); only Tyr/Trp survive both |
| **A5** Epitope logic | Is an epitope chemically special? | Interface vs its *own* protein's surface | paratope **0.969**, PPI 0.492, **epitope 0.258**; protrusion +0.293 vs **−0.036** |
| **A3** Shape complementarity | Are antibodies less complementary? | Sc on contact surface | 0.464 vs 0.467, **p = 0.082 (n.s.)** |
| **A4** Geometry + PCA | Are they a separable class? | 20-descriptor space + patch geometry | Rg **10.7 vs 16.3 Å**, planarity 0.85 vs 0.79; **AUC 0.875 but 93.8% overlap** |
| **A6** Functional sites | Do antibodies target functional sites? | Ligand-proximity site overlap | **1.6% vs 10.7%** (effect −0.45, p=8e−09) → antibodies avoid them |

## Case studies — the same question at single-system resolution

| Study | Asks | Headline result |
|---|---|---|
| **Immunodominance** | Where do 42 independent antibodies hit one antigen? | **Gly102 in 98%**, Gini **0.49**, one contiguous face on the surface render |
| **Lysozyme** | Active site or its rim? Conserved or variable? | catalytic Glu35 **12%** / Asp52 **19%** vs substrate-binding Asp101 **69%**; immunodominant positions less conserved (p = 0.078, **trend**) |
| **Matched targets** | Same protein, antibody vs natural partner — same face? | 16 pairs: median Jaccard **0.033**, **38% share no position** |
| **PyMOL mechanism** | Is the Tyr-OH claim real at atomic resolution? | **144 bonds in 83 sub-1.8 Å structures**, rendered at 1.10–1.30 Å (2.64–2.94 Å) |

## Synthesis

| Phase | Gives |
|---|---|
| **10** Literature concordance | 16 topics: 10 confirm, 3 refine, 3 contradict → `LITERATURE_CONCORDANCE.md` |
| **—** Findings register | 10 confirm / 3 refine / 3 contradict / 15 new / 5 methodological → `FINDINGS_VS_LITERATURE.md` |
| **—** Manuscript plan | 8 figures, each mapped to the question it answers → `MANUSCRIPT_PLAN.md` |

---

## The one-paragraph answer

Antibodies and co-evolved protein partners share their chemistry at the level of
chemical classes (r = 0.84) but use **formally distinct residue-level grammars**
(p = 0.001) that survive every confound tested. The difference is **not** in what
fraction of the surface they read by side chain (identical, 0.72) but in **how they
read backbone** (cross-arm r = 0.14 vs 0.49), mechanised by **tyrosine hydroxyl
hydrogen-bonding to the antigen main chain** plus cation–π — not π-stacking.
Antibodies achieve this over a **smaller, flatter, less hydrophobically driven**
interface (ΔCp −94 vs −151 cal/mol/K per 1000 Å²). Meanwhile the **epitope is
ordinary**: chemically indistinguishable from the rest of its own antigen's surface
(0.258 vs the paratope's 0.969), non-protruding, on a different face from the natural
partner (Jaccard 0.033), and avoiding functional sites (1.6% vs 10.7%). The
specialised object is the antibody; the epitope is whatever surface was available.
The difference is nonetheless bounded — a general-PPI potential still scores antibody
interfaces well (AUC 0.833).
