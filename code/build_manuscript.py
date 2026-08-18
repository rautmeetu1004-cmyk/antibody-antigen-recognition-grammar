"""Build the Springer-format manuscript and its supplementary file as .docx.

Numbers are pulled from metadata/analysis/*.json wherever practical so the text
cannot drift from the results it describes.
"""
from __future__ import annotations

import json
import pathlib
import sys

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

M = C.OUT_MET
FIG = C.OUT_FIG


def J(name):
    p = M / name
    return json.loads(p.read_text()) if p.exists() else {}


# ----------------------------------------------------------------- doc scaffolding
def new_doc():
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)
    for s in doc.sections:
        s.top_margin = Inches(1)
        s.bottom_margin = Inches(1)
        s.left_margin = Inches(1)
        s.right_margin = Inches(1)
    return doc


def para(doc, text, size=11, italic=False, bold=False, align=None, space=6):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.size = Pt(size)
    r.italic = italic
    r.bold = bold
    p.paragraph_format.space_after = Pt(space)
    if align:
        p.alignment = align
    return p


def heading(doc, text, level=1):
    h = doc.add_heading(text, level=level)
    for r in h.runs:
        r.font.color.rgb = RGBColor(0, 0, 0)
        r.font.name = "Times New Roman"
    return h


def table(doc, rows, header=True, widths=None):
    t = doc.add_table(rows=0, cols=len(rows[0]))
    t.style = "Table Grid"
    for i, row in enumerate(rows):
        cells = t.add_row().cells
        for j, val in enumerate(row):
            cells[j].text = ""
            p = cells[j].paragraphs[0]
            r = p.add_run(str(val))
            r.font.size = Pt(9)
            r.bold = (i == 0 and header)
    return t


def figure(doc, fname, caption, width=6.2):
    path = FIG / fname
    if path.exists():
        doc.add_picture(str(path), width=Inches(width))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    p = doc.add_paragraph()
    r = p.add_run(caption)
    r.font.size = Pt(9)
    p.paragraph_format.space_after = Pt(12)


# ------------------------------------------------------------------------- content
def build_manuscript():
    hm = J("heatmap.json")
    inv = J("inventory.json")
    cons = J("conservation.json")
    rd = J("readout.json")
    conf = J("conformation.json")
    para_j = J("paratope.json")
    diff = J("differential.json")
    ctrl = J("controls.json")
    mod = J("modified.json")
    thermo = J("thermodynamics.json")
    epi = J("epitope_logic.json")
    fs = J("functional_sites.json")
    mt = J("matched_targets.json")
    pca = J("pca_geometry.json")
    sc = J("shape_complementarity.json")
    tr = J("transfer.json")
    sp = J("structural_pass.json")
    imm = J("immunodominance.json")

    doc = new_doc()

    # ---------------------------------------------------------------- title block
    p = doc.add_paragraph()
    r = p.add_run("Antibodies read protein sequence conventionally and protein "
                  "shape unconventionally: a matched comparison of antibody and "
                  "co-evolved protein interfaces")
    r.bold = True
    r.font.size = Pt(15)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    para(doc, "", size=6)
    para(doc, "Meetrayu Raut, Amrita Sinharay, Anup Kale",
         size=11, align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, "Computational Biotechnology Laboratory",
         size=10, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, "Correspondence: rautk2024@gmail.com",
         size=10, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)

    # ------------------------------------------------------------------- abstract
    heading(doc, "Abstract", 1)
    para(doc,
         "Antibodies bind protein surfaces with which they have no evolutionary "
         "history, and are thought to exploit physicochemical features common to all "
         "protein surfaces rather than features peculiar to their targets. Testing "
         "that proposition requires measuring antibody and co-evolved interfaces by "
         "one procedure, which is rarely done. The absence of a matched comparator "
         "carries a specific cost: several important observations about epitopes are "
         "null results, and a null cannot be interpreted without knowing what a "
         "positive one looks like. We measured 512 antibody-antigen and 2,292 "
         "protein-protein complexes with a single pipeline, weighting every statistic "
         "over 30 percent sequence-identity clusters. The two classes agree closely "
         "when chemistry is described by chemical class (r = 0.84) and differ at "
         "residue-pair resolution, where a permutation test against the null of one "
         "shared grammar gives 0.42 against a null of 0.80 (p = 0.001). The "
         "difference survives matching on secondary structure, fold class, "
         "resolution, redundancy, hydrogen-bond criterion and antibody format. "
         "Decomposing contacts by the moiety being read localises it. Both classes "
         "read side chains in the same proportion, 0.72, and their side-chain "
         "grammars correlate at 0.49, but their main-chain grammars correlate at "
         "0.14. Functional-group analysis attributes this to the tyrosine phenol "
         "hydroxyl hydrogen bonding to the antigen main chain, while aromatic "
         "stacking is depleted. Antibody interfaces are smaller, flatter and less "
         "hydrophobically driven, with a heat-capacity change of -94 against -151 cal "
         "per mol per K per 1000 square angstroms buried. Measuring three interface "
         "types by one procedure also places existing null results on a scale: "
         "compositional distinguishability from a protein's own accessible surface is "
         "0.97 for the paratope and 0.48 for a co-evolved interface, but 0.25 for the "
         "epitope. The specialised component of an antibody-antigen interface is the "
         "antibody.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    para(doc, "Keywords: antibody-antigen recognition; protein-protein interfaces; "
              "structural bioinformatics; interface chemistry; epitope; paratope",
         size=10)

    # --------------------------------------------------------------- introduction
    heading(doc, "1. Introduction", 1)
    para(doc,
         "The antibody-antigen interface occupies an unusual position in structural "
         "biology. It is among the most intensively solved classes of protein "
         "complex, but it is also the least evolutionarily conserved [22, 25]. "
         "The production of antibodies revolves around a germline repertoire, "
         "diversified by somatic hypermutation, which produces binders against "
         "surfaces that have no evolutionary history of being bound [26]. This is the "
         "opposite of the situation in an obligate oligomer or an enzyme-inhibitor "
         "pair, where both partners have been shaped by selection acting on the "
         "interface itself [9, 14]. The question that follows is simple to state but "
         "surprisingly hard to answer. When an antibody reads a protein surface, does "
         "it use the same chemical rules that co-evolved partners use?",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Both interface classes have been described in detail, but rarely against "
         "each other on a common footing. Statements that antibody interfaces are "
         "smaller, more polar or less shape-complementary are typically assembled "
         "from separate studies with different contact definitions, different "
         "structural corpora and, most importantly, different redundancy treatment "
         "[5, 12, 15, 20]. The consequence is not merely imprecision. For instance, a "
         "corpus in which forty-five lysozyme complexes are counted as forty-five "
         "independent observations will report the chemistry of lysozyme rather than "
         "the chemistry of antibodies.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The absence of a matched comparator has a specific and under-appreciated "
         "cost, which is that null results become uninterpretable. Analyses of B-cell "
         "epitopes have repeatedly found that epitope composition does not differ "
         "significantly from equally exposed antigen surface, and that epitope "
         "residues cannot be related to accessibility or to a protrusion index "
         "[11, 21]. These are informative observations, and they are among the "
         "reasons that structure-based epitope prediction has proved difficult. What "
         "they cannot establish on their own is whether the absence of a difference "
         "makes epitopes unusual, because an interface known to be shaped by "
         "selection was not measured by the same procedure for comparison. Similar "
         "reasoning applies to reports that paratope aromatic residues contact "
         "antigen main chain atoms [15, 19]. The observation is well founded, but the "
         "quantity it should be compared against has been missing.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Additionally, antibody paratopes are formed by hypervariable loops [13, 16], "
         "and antibody antigens are drawn disproportionately from beta-rich viral "
         "surface proteins. Any chemical difference between antibody and general "
         "interfaces could therefore reflect what is being bound rather than how it "
         "is being bound. This confound is widely acknowledged and seldom tested. "
         "Testing it requires comparing the two classes within matched structural "
         "context, and it requires a null model for that comparison, because "
         "restricting data to a matched stratum reduces sample size and attenuates "
         "correlations for purely statistical reasons.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "In this paper, we addressed both difficulties directly, the missing "
         "comparator and the structural confound. A single extraction pipeline was "
         "applied to antibody-antigen and general protein-protein complexes, and all "
         "statistics were weighted over independence units defined by 30 percent "
         "sequence identity clustering of the bound protein. Every stratified "
         "comparison was benchmarked against a null drawn at the same sample size. We "
         "then asked three questions in sequence. First, what physical forces hold "
         "each interface class together. Second, whether the residue-pair "
         "preferences, which we refer to as the recognition grammar, are the same in "
         "both classes. Third, if they differ, whether the difference can be "
         "attributed to structural context, and where it cannot, which chemistry is "
         "responsible.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # -------------------------------------------------------------------- methods
    heading(doc, "2. Materials and Methods", 1)

    heading(doc, "2.1 Corpus construction", 2)
    para(doc,
         "Two groups were assembled from experimentally determined structures. "
         "Group 1 comprised antibody-antigen complexes with a folded protein "
         "antigen. Antibody chains were identified by an independent sequence rule "
         "and cross-checked against entity description text, with the sequence rule "
         "taking precedence and the description serving as an audit of it. The two "
         "detectors agreed on 93.4 percent of scored chains. Complexes containing "
         "T-cell receptor or major histocompatibility components were excluded, as "
         "were entries without an identifiable antigen. A floor of 50 residues was "
         "applied to the antigen, because short peptide epitopes are not protein "
         "surfaces and would not support the comparison this study makes. "
         "Sensitivity of the corpus to that floor was recorded at 0, 30, 50 and 80 "
         "residues, giving 774, 530, 512 and 481 structures respectively. The 50 "
         "residue floor was adopted, yielding 512 structures.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Group 2 comprised general protein-protein complexes with no antibody "
         "component, giving 2,292 structures that entered the analysis. Because two "
         "co-evolved partners have no natural probe and target roles, each complex "
         "was extracted in both directions and flagged accordingly. This doubling is "
         "correct for unit-weighted matrices, where each half-interface is clustered "
         "by the identity of its own protein and each deserves a vote, but it is "
         "incorrect for any test that requires a true sample size. Deduplication to "
         "a single direction was therefore applied before all Fisher exact tests and "
         "all per-interface distributions.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Antibody metadata were enriched from the Structural Antibody Database, "
         "which classifies constructs by scanning chains for numerable "
         "immunoglobulin variable domains rather than by text matching. Of the 512 "
         "structures in Group 1, 506 carried a construct annotation: 254 Fab, 203 "
         "single-domain heavy chain, 48 Fv and 1 single-domain light chain. The two "
         "groups are closely matched in data quality, with median resolutions of "
         "1.83 and 1.80 angstroms.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "2.2 Feature extraction", 2)
    para(doc,
         "Feature extraction is the foundation of every result reported here, and it "
         "was performed identically for both groups. Three layers were produced.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The first layer is the residue-pair contact table. A contact was recorded "
         "between two residues on opposite sides of an interface when any pair of "
         "their heavy atoms fell within 4.5 angstroms, a threshold in standard "
         "use for interface definition and consistent with the contact analyses "
         "of Jones and Thornton [9] and of Chakrabarti and Janin [5]. Each contact was annotated "
         "with the interaction types it supports: a salt bridge when a basic "
         "side-chain nitrogen lay within 4.0 angstroms of a carboxylate oxygen, "
         "following the ion-pair criterion of Barlow and Thornton [1]; a hydrogen bond "
         "when a nitrogen or oxygen donor-acceptor pair fell within 3.5 angstroms; a "
         "van der Waals contact within 4.5 angstroms; and aromatic stacking when "
         "ring centroids lay within 5.5 angstroms at an interplanar angle below 35 "
         "degrees. Each contact was further annotated by the moiety of the target "
         "residue involved, separating side-chain from main-chain mediation. This "
         "separation is what later allows sequence readout and shape readout to be "
         "analysed as independent channels. Extraction produced 36,810 residue pairs "
         "for Group 1 and 406,397 for Group 2 after deduplication.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The second layer is the atom-level contact table, retaining the identity "
         "of both atoms in every contacting pair. This layer supports analysis at "
         "the resolution of chemical functional groups rather than residue types, "
         "which matters because a residue label conflates moieties that behave "
         "differently. Tyrosine, for example, presents both an aromatic ring and a "
         "hydrogen-bonding hydroxyl, and these participate in different chemistry. "
         "Every heavy atom was assigned to one of fourteen functional groups: "
         "guanidinium, ammonium, carboxylate, amide, hydroxyl, phenol hydroxyl, "
         "aromatic ring, indole nitrogen, indole ring, imidazole, thiol or "
         "thioether, aliphatic carbon, backbone amide and backbone carbonyl. "
         "Backbone was deliberately split into its donor and acceptor because those "
         "groups behave differently. This layer comprises 182,442 mapped atom "
         "contacts for Group 1 and 1,767,916 for Group 2.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The third layer is per-residue structural context. Because this pass "
         "requires repeated solvent-accessibility calculation it was run over "
         "1,512 complexes rather than the full corpus: all 512 antibody "
         "complexes, and a seeded random sample of 1,000 of the 2,292 "
         "protein-protein complexes, which still exceeds the antibody arm in both "
         "structures and independence units. Every analysis that depends on "
         "secondary structure, torsion angles or buried surface therefore uses "
         "this subsample on the protein-protein side, while the contact and "
         "atom-level analyses use the complete corpus. The pass yielded "
         "1,019,571 annotated residues. "
         "Secondary structure was assigned by the Kabsch and Sander algorithm [10], "
         "implemented directly for this study because no DSSP executable was "
         "available in the analysis environment. The implementation places the amide "
         "hydrogen along the direction opposite the preceding carbonyl, computes the "
         "electrostatic hydrogen-bond energy with the standard coefficient, applies "
         "a threshold of -0.5 kcal per mol, and derives turns, bridges and ladders "
         "from the resulting hydrogen-bond map. Because an implementation of this "
         "kind should not be trusted without evidence, it was validated against the "
         "HELIX and SHEET records deposited with each structure, giving 88.9 percent "
         "three-state agreement across 622,662 residues. The implementation was "
         "additionally validated against two independent references. Against pydssp, "
         "a separate implementation of the same algorithm, three-state agreement "
         "was 98.25 percent over 432 chains and 79,977 residues. Against the DSSP "
         "program itself, using annotations computed on the complete structure and "
         "obtained from PDB-REDO, agreement was 97.67 percent with a Cohen kappa "
         "of 0.965 over 65,336 residues in 97 structures, with residues matched "
         "between the two files by alpha-carbon coordinate rather than by residue "
         "number. The second comparison is the more informative because DSSP "
         "operates on the whole structure while our implementation operates per "
         "chain. Recall against DSSP was 0.998 for helix and 0.950 for strand, and "
         "the disagreement is almost entirely one-directional: 940 residues that "
         "DSSP assigns to strand are assigned to loop here, against 3 in the "
         "opposite direction. That asymmetry is the expected signature of "
         "per-chain assignment, since a bridge partner in another chain cannot be "
         "seen, and it bounds the effect at five percent of strand residues. "
         "Agreement with deposited HELIX and SHEET records is lower still, at 88.9 "
         "percent, but that figure is bounded by the reference rather than by the "
         "code, since those records are author annotations with no loop category. It is lower "
         "than the agreement typically quoted between independent secondary-"
         "structure programs, which exceeds 95 percent, and the discrepancy "
         "requires comment. Deposited HELIX and SHEET records are author "
         "annotations rather than a geometric derivation, they are known to "
         "differ from algorithmic assignment particularly at element termini "
         "and for short or distorted elements, and they provide no "
         "loop category against which agreement can be scored. Our figure is "
         "therefore agreement with a coarser reference, not with DSSP itself, and "
         "should be read as a lower bound. All secondary-structure claims in "
         "this work are three-state and comparative, applied identically to "
         "both groups, so a uniform assignment error would affect both arms "
         "equally and cannot generate a difference between them. Backbone torsions phi and "
         "psi and side-chain torsions chi1 and chi2 were computed from coordinates, "
         "with the dihedral routine asserted against an independent reference "
         "implementation before use. Solvent-accessible surface area was computed by "
         "the Shrake-Rupley algorithm [23] in the complex and in isolation, giving the "
         "surface buried on association for each residue and, in a separate pass "
         "resolved to atoms, separately for apolar and polar atoms.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "2.3 Independence units and weighting", 2)
    para(doc,
         "Structural corpora are heavily redundant, and the redundancy is not "
         "random. The antibody corpus contains 45 complexes of a single lysozyme "
         "cluster and 39 of a spike protein cluster. Every statistic in this study "
         "is therefore weighted over independence units defined by clustering the "
         "bound protein at 30 percent sequence identity. Each contact carries a "
         "weight equal to the reciprocal of the number of structures in its unit, so "
         "that a cluster represented by 45 structures contributes as much as a "
         "cluster represented by one. Group 1 resolves to 233 units and Group 2 to "
         "1,326. A parallel analysis using UniProt accession as the unit was run "
         "throughout as a sensitivity check. Resampling, where used, resamples units "
         "rather than contacts, because contacts within a structure are strongly "
         "correlated and a contact-level resample would report intervals that are "
         "far too narrow.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "2.4 Statistical framework", 2)
    para(doc,
         "Each statistic below is stated with the reason it was chosen, since the "
         "choice of null model determines what a number means.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Recognition grammar. For each group we constructed a 20 by 20 matrix whose "
         "entry is the unit-weighted count of contacts between a given residue on "
         "each side. Counts alone are uninformative because some residues are simply "
         "common at interfaces. Each cell was therefore converted to the base-2 "
         "logarithm of the ratio of observed count to the count expected from the "
         "product of the row and column marginals. This quantity measures "
         "preference: how much more often a pair is in contact than the abundance of "
         "its two residues predicts. The antibody matrix was kept asymmetric, "
         "because paratope and epitope are not exchangeable roles, while the "
         "protein-protein matrix was symmetrised, because two co-evolved partners "
         "are.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Effect size and evidence were kept in separate channels throughout. Effect "
         "sizes were computed from unit-weighted counts, which express the biology "
         "we intend to measure. Significance was computed from raw counts using "
         "Fisher exact tests, because a test requires a true sample size and feeding "
         "fractional weights to it would manufacture significance. Multiple testing "
         "across the 400 cells of each matrix was controlled by the "
         "Benjamini-Hochberg procedure.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Confidence intervals on matrix quantities were obtained by resampling "
         "independence units with replacement 2,000 times and recomputing the "
         "statistic. This was implemented as a multinomial draw over precomputed "
         "per-unit matrices followed by a single matrix product, which is "
         "mathematically identical to explicit resampling and was verified as such. "
         "We note that a percentile interval on a correlation is attenuated by "
         "resampling noise, so such intervals are reported as stability intervals "
         "rather than unbiased confidence bounds.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Test of grammar identity. To ask whether the two matrices could be samples "
         "of one underlying grammar, independence units from both groups were pooled "
         "and randomly reassigned to two groups of the original sizes, and the "
         "cross-group correlation recomputed, 1,000 times. This permutation null is "
         "constructed in the same way as the observation, so unlike the bootstrap it "
         "is not attenuated, and it tests exactly the hypothesis of interest.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Size-matched nulls. Several analyses restrict both groups to a matched "
         "context, such as loop residues only or a single fold class. Restriction "
         "reduces sample size, and small samples attenuate correlations on their "
         "own, so a correlation that falls under stratification cannot be "
         "interpreted without a reference. For every stratified comparison we drew "
         "random subsets of the unrestricted data at the same size and recomputed "
         "the statistic, giving the value expected from sample size alone. Three "
         "results in this study would have been misread without this control.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Reliability ceilings. Comparing correlations computed at different sample "
         "sizes compares precision rather than biology. Where two subgroups of "
         "unequal size were compared, we established the ceiling by splitting the "
         "larger subgroup in half at the smaller subgroup's size and correlating the "
         "halves with each other.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Per-interface comparisons. Distributional comparisons between the two "
         "groups used the Mann-Whitney U test on unit-level means, so that each "
         "independence unit contributes one observation. Effect sizes are reported "
         "as the rank-biserial correlation, which is zero when two distributions "
         "coincide and plus or minus one when they are completely separated. This "
         "quantity is reported alongside every p value, because with corpora of this "
         "size statistical significance is easy to obtain and is not by itself "
         "informative.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Burial thermodynamics. The change in heat capacity on association was "
         "estimated from buried surface using the parameterisation of Murphy and "
         "Freire [17], in which apolar and polar buried areas enter with coefficients of "
         "0.45 and -0.26 cal per mol per K per square angstrom. Because the relation "
         "is written in terms of the change in accessible area, which is negative on "
         "binding, buried areas enter with a negative sign. A large negative heat "
         "capacity change is the classical signature of hydrophobically driven, "
         "entropy-dominated association. Values were also normalised per 1000 square "
         "angstroms buried, since antibody interfaces are smaller and an unnormalised "
         "comparison would partly restate that.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Interface scoring and transfer. Each matrix was used as a threading "
         "potential in which an interface score is the mean log-odds over its "
         "contacts. Decoys preserved the interface geometry, the contact graph and "
         "the residue composition, permuting only which residue occupied which "
         "position, so that discrimination reflects recognition chemistry rather "
         "than interface size. Discrimination is reported as the area under the "
         "curve over 50 decoys per interface, averaged within unit. The potential "
         "used to score a given interface excluded that interface's own unit.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Shape complementarity was computed in the manner of Lawrence and Colman [12], "
         "using the negative dot product of surface normals at nearest-neighbour "
         "points weighted by a Gaussian in separation, with the weighting constant "
         "set to 0.5 per square angstrom. Points were placed on the van der Waals "
         "contact surface with accessibility determined by a 1.4 angstrom probe. "
         "This approximates the molecular surface used in the original definition "
         "but omits re-entrant regions, so absolute values are not comparable with "
         "published figures. The comparison between groups, computed identically for "
         "both, is what the analysis supports.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "2.5 Confound controls", 2)
    para(doc,
         "Six potential artefacts were tested. Secondary structure was controlled by "
         "recomputing the cross-group correlation within matched assignments, "
         "including loop residues on both sides. Fold composition was controlled "
         "using CATH class assignments retrieved for the bound chain. Resolution was "
         "controlled by stratifying both groups above and below 2.0 angstroms. "
         "Redundancy was controlled by dropping the largest antigen clusters and by "
         "leaving out each of the twenty largest clusters in turn. The hydrogen-bond "
         "criterion was controlled by rederiving hydrogen-bond matrices from atom "
         "coordinates at 3.0, 3.2 and 3.5 angstroms. Antibody format was controlled "
         "by comparing Fab and single-domain constructs against the reliability "
         "ceiling described above. Every stratified comparison carries its "
         "size-matched null.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "2.6 Verification", 2)
    para(doc,
         "Contacts were re-derived from raw coordinates by an independently written "
         "extractor for a random sample of structures, giving F1 agreement of 0.9990 "
         "for Group 1 across 2,048 residue pairs and 1.0000 for Group 2 across 2,563 "
         "pairs, with a median deviation in closest-approach distance of 0.0025 "
         "angstroms. Thirty-nine further automated checks re-derive reported "
         "quantities by independent routes, including an independent "
         "reimplementation of the multiple-testing procedure and confirmation that "
         "the vectorised unit bootstrap reproduces explicit resampling. These are "
         "documented in the supplementary information.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "2.7 Software", 2)
    para(doc,
         "Analyses were implemented in Python 3.14.6 with NumPy 2.4.6, SciPy 1.18.0, "
         "scikit-learn 1.9.0, Biopython 1.87 and Matplotlib 3.11.0. Structural "
         "rendering used PyMOL 3.1.0 and document generation python-docx 1.2.0. Fold "
         "assignments and antibody metadata were retrieved from public web services. "
         "All code, intermediate tables and figures are provided.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # -------------------------------------------------------- results & discussion
    heading(doc, "3. Results and Discussion", 1)

    heading(doc, "3.1 Corpus quality and the comparability of the two groups", 2)
    para(doc,
         "The two groups are well matched in the respects that matter for a chemical "
         "comparison. Median resolution is 1.83 angstroms for antibody complexes and "
         "1.80 for protein-protein complexes, so contact assignments are of "
         "comparable reliability. Redundancy is present but tractable: the largest "
         "antigen cluster accounts for 8.8 percent of Group 1, and 149 of its 233 "
         "units are singletons.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "One property of the corpus deserves emphasis because it bounds the "
         "generality of everything that follows. Group 1 consists entirely of X-ray "
         "structures. The Structural Antibody Database contains several thousand "
         "cryo-electron microscopy entries with protein antigens, and none of them "
         "is represented here. Those entries are enriched in large viral "
         "glycoproteins, which are precisely the targets of much contemporary "
         "therapeutic antibody work. The corpus is therefore of unusually high and "
         "uniform quality, which benefits the chemistry, and it is not "
         "representative of the field as a whole.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "3.2 The physical forces at each interface", 2)
    para(doc,
         "Both interface classes are built from the same two dominant interactions. "
         "Van der Waals contacts account for 64 percent of antibody contacts and 69 "
         "percent of protein-protein contacts, and hydrogen bonds for 25 and 20 "
         "percent respectively. Salt bridges are a minor component in both, at "
         "roughly 3 percent, and formal aromatic stacking is negligible at 0.3 "
         "percent.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Where the classes differ, the pattern is a substitution rather than an "
         "addition. Antibody interfaces show markedly higher aromatic involvement "
         "and markedly lower apolar contact, with rank-biserial effects of 0.53 and "
         "-0.52. Hydrogen bonding is more frequent, with an effect of 0.44, as are "
         "cation-pi contacts [8], at 0.34. The interfaces are smaller, with a median of "
         "59 residue pairs against 87, yet closer packed, with a mean "
         "closest-approach distance of 3.59 against 3.64 angstroms. Both figures sit "
         "within the ranges reported for the two classes separately [14, 20].",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Several comparisons are notable for showing no difference, and these require "
         "a different form of evidence. A non-significant difference test is "
         "absence of evidence rather than evidence of absence, so equivalence was "
         "tested directly by two one-sided tests against a margin set at 20 percent "
         "of the largest composition difference this study measures on the same "
         "scale. Salt-bridge frequency, contact frequency between oppositely "
         "charged residues and aromatic stacking frequency are all statistically "
         "equivalent between the two classes on that criterion, with equivalence p "
         "below 0.0001 in each case and bootstrap intervals on the difference "
         "spanning less than 0.013 throughout. Charge complementarity is therefore "
         "general protein interface physics rather than an antibody characteristic, "
         "and an account of what distinguishes antibodies must look elsewhere. This "
         "is a positive finding about sameness, not a failure to detect a "
         "difference.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The composition difference has a thermodynamic counterpart. The fraction "
         "of buried surface that is apolar is 0.50 for antibody interfaces and 0.58 "
         "for protein-protein interfaces. The corresponding heat-capacity change, "
         "normalised per 1000 square angstroms buried, is -94 against -151 cal per "
         "mol per K, with an effect of 0.52. A less negative heat-capacity change "
         "indicates a smaller contribution from the classical hydrophobic effect. "
         "Antibody association is therefore measurably less entropically driven than "
         "association between co-evolved partners, which is consistent with a binder "
         "that cannot co-evolve a complementary hydrophobic core with its target and "
         "must instead purchase affinity through directional interactions [24].",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "phase2_physicochemistry.png",
           "Fig. 1 Physicochemical composition of antibody-antigen and general "
           "protein-protein interfaces. Each panel gives the distribution, across "
           "cluster30 independence units, of the fraction of interface contacts of "
           "one type: (a) hydrogen bonds, (b) salt bridges, (c) van der Waals "
           "contacts, (d) aromatic stacking, (e) opposite-charge pairs, (f) "
           "like-charge pairs, (g) apolar-apolar pairs, (h) contacts involving an "
           "aromatic residue, (i) cation-pi contacts, and (j) contacts read through "
           "the target residue side chain. Red is antibody-antigen and blue is "
           "general protein-protein, with dashed lines at the medians. Effect sizes "
           "quoted in the text are rank-biserial correlations.")

    heading(doc, "3.3 The recognition grammar", 2)
    para(doc,
         "Collapsing the 20 by 20 matrices to six chemical classes gives a "
         "correlation of 0.84 between the two groups. At this level the classes "
         "agree: charged residues pair with oppositely charged residues, apolar with "
         "apolar, and the pattern is shared. At the level of individual residue "
         "pairs the correlation is 0.42, rising to 0.59 when cysteine is excluded. "
         "Cysteine is the single most divergent residue between the groups, because "
         "inter-chain disulfides are frequent in protein-protein complexes and "
         "essentially absent from paratopes, and its exclusion increases rather than "
         "decreases the measured agreement.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The gap between 0.84 and 0.42 is the central observation of this study. "
         "The two classes obey the same chemistry when chemistry is described "
         "coarsely, and they diverge when it is described at the resolution of "
         "residue identity. A permutation test resolves whether that divergence is "
         "real. Pooling independence units from both groups and reassigning them at "
         "random gives cross-group correlations centred on 0.80, with a 95 percent "
         "range of 0.73 to 0.85. The observed value of 0.42 lies far below this "
         "range, with p = 0.001 at 1,000 permutations. The two grammars are not two "
         "samples of one grammar.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The composition of the shared portion is informative. Of the 55 cells that "
         "are significant and sign-stable in both groups, the majority are "
         "depletions rather than enrichments. Pairs such as lysine with alanine, "
         "glutamate with aspartate and lysine with arginine are avoided in both "
         "classes, while the shared enrichments are confined largely to charge "
         "complementarity. What protein interfaces agree on is principally what they "
         "will not do.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "heatmap_20x20_cluster30.png",
           "Fig. 2 Amino-acid recognition matrices. (a) Antibody-antigen, with "
           "paratope residue on the vertical axis and epitope residue on the "
           "horizontal axis. (b) General protein-protein, with the two exchangeable "
           "partners on the two axes. Colour is the base-2 logarithm "
           "of observed over expected contact frequency, where expected is the "
           "product of marginals. Crosses mark cells failing Benjamini-Hochberg "
           "control at q < 0.05 across 400 tests; outlines mark cells whose sign is "
           "stable in at least 95 percent of 2,000 unit resamples.")

    heading(doc, "3.4 The difference lies in backbone readout", 2)
    para(doc,
         "A protein surface can be read in two ways that carry different "
         "information. Contacts to side chains report residue identity, and are "
         "therefore sequence readout. Contacts to the main chain report local "
         "conformation, since the backbone amide and carbonyl are chemically "
         "identical in every residue except proline and glycine, and are therefore "
         "shape readout. Separating these channels was originally motivated as a "
         "check on the interpretability of the matrices. Had most contacts been "
         "main-chain mediated, a residue-identity matrix would have been recording "
         "which residues present their backbone rather than which chemistry is "
         "preferred.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The check returned a near-null. Both groups read close to 0.72 of their "
         "contacts through target side chains, and the per-residue pattern of "
         "side-chain against main-chain reading correlates between the groups at "
         "0.98. The difference in fraction is small but not zero: antibodies read "
         "0.016 more of their contacts through side chains, with a bootstrap "
         "interval on the difference of 0.002 to 0.032 that excludes zero, and "
         "equivalence at our stated margin is not established for this quantity. "
         "The relevant point is one of magnitude. A shift of 1.6 percentage points "
         "in how much backbone is read is an order of magnitude smaller than the "
         "divergence in how it is read, reported next.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The channels nevertheless behave very differently across the two classes. "
         "When the matrices are rebuilt separately for each channel, the "
         "cross-group correlation is 0.49 for side-chain mediated contacts and "
         "0.14 for main-chain mediated contacts. The divergence between antibody "
         "and general interfaces therefore resides almost entirely in how the "
         "backbone is read rather than in how side chains are read.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Within each class the two channels also carry distinguishable "
         "information, and quantifying that required its own null. Permuting the "
         "moiety label of every contact within its interface, which preserves "
         "interface identity, residue composition and the overall side to main "
         "ratio, gives the correlation expected if moiety carried no grammar. Both "
         "classes fall below their permuted null, at 0.28 against 0.53 for "
         "antibodies and 0.55 against 0.94 for protein-protein interfaces, each at "
         "p = 0.005. Both classes therefore read sequence and shape through "
         "channels that encode genuinely different preferences. We note "
         "specifically that the raw difference between 0.28 and 0.55 should not be "
         "read as antibody channels being more independent than co-evolved ones. "
         "The two nulls differ substantially because the protein-protein corpus is "
         "an order of magnitude larger, and relative to their own nulls the two "
         "classes are comparably decoupled.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "phase5_readout.png",
           "Fig. 3 Sequence readout and shape readout as separate channels. Each "
           "matrix is restricted to contacts mediated by one moiety of the target "
           "residue. (a) Antibody-antigen, side-chain-mediated. (b) "
           "Antibody-antigen, main-chain-mediated. (c) General protein-protein, "
           "side-chain-mediated. (d) General protein-protein, main-chain-mediated. "
           "Colour is the base-2 logarithm of observed over expected contact "
           "frequency. (e) Fraction of contacts each target residue makes through "
           "its side chain in the two classes; residues left of the dashed line are "
           "read mainly through the backbone.")

    heading(doc, "3.5 The chemistry responsible", 2)
    para(doc,
         "Descending from residues to functional groups identifies the interaction "
         "that distinguishes the two classes. Across 1.95 million atom-atom "
         "contacts, the functional-group maps of the two groups correlate at 0.62, "
         "higher than the residue-level correlation, which indicates that chemistry "
         "is more conserved than residue identity. The largest antibody-specific "
         "enrichments are contacts between an aromatic ring and an ammonium group, "
         "at 1.07, and contacts between the tyrosine phenol hydroxyl and the antigen "
         "main chain, at 0.56 to the backbone amide and 0.58 to the backbone "
         "carbonyl.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "That paratope aromatic side chains contact antigen main chain atoms has "
         "been reported before, in analyses of the origins of antibody specificity [15, 19]. "
         "What a matched comparison adds is a baseline, a moiety and a negative "
         "control. The baseline shows how much more often this occurs than at a "
         "general protein interface. The moiety is specifically the phenol hydroxyl "
         "rather than the ring, which matters because the two do different chemistry. "
         "And the negative control is aromatic stacking, which is depleted rather "
         "than enriched in antibody interfaces: tryptophan-tyrosine and "
         "tyrosine-tyrosine pairs differ by -0.99 and -0.85 from protein-protein "
         "interfaces, and formal stacking accounts for 0.3 percent of contacts in "
         "both groups. Aromatic enrichment in a paratope therefore does not mean "
         "stacking. The aromatic residues contact backbone and cationic groups, not "
         "one another.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Taken with the channel analysis of the previous section, this identifies "
         "the chemistry behind the backbone-readout divergence. Antibodies read the "
         "antigen main chain using tyrosine hydroxyls as hydrogen-bond donors and "
         "acceptors, an interaction that couples a directional hydrogen bond to the "
         "packing and polarisability of a large rigid ring, and that is available to "
         "a germline-encoded loop without requiring any complementary adaptation in "
         "the target.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Because a statistical enrichment is not the same as a demonstrated "
         "interaction, we examined the geometry directly. The corpus contains 144 "
         "hydrogen bonds between a tyrosine hydroxyl and an antigen main-chain "
         "nitrogen or oxygen in structures at better than 1.8 angstroms, distributed "
         "across 83 structures. Their distances lie between 2.6 and 3.1 angstroms, "
         "which is unambiguous hydrogen-bonding geometry.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "pymol_mechanism_montage.png",
           "Fig. 4 Tyrosine hydroxyl hydrogen bonds to the antigen main chain in "
           "three high-resolution complexes: (a) 3D9A at 1.20 angstroms, a lysozyme "
           "C complex, (b) 2P45 at 1.10 angstroms, a ribonuclease complex, and (c) "
           "7AQZ at 1.30 angstroms, a nanobody complex. Antibody is shown in blue "
           "and antigen "
           "in wheat, with the tyrosine side chain and the contacted backbone atoms "
           "as sticks.")

    heading(doc, "3.6 The difference is not an artefact of structural context", 2)
    para(doc,
         "Antibody paratopes are loop-dominated. Secondary-structure assignment from "
         "DSSP confirms this: 64 percent of paratope contacts are contributed by "
         "loop residues, against 42 percent for the probe side of protein-protein "
         "complexes. The obvious concern is that the recognition grammar reported "
         "above is loop grammar rather than antibody grammar.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Restricting both groups to loop residues on both sides does not reduce the "
         "difference. It increases it. Using secondary structure from the DSSP "
         "program, available for 2,767 of the 2,804 entries and covering 97 "
         "percent of contacts in each group, the cross-group correlation falls "
         "from 0.387 for all contacts to 0.190 for loop-to-loop contacts. Because "
         "restriction also reduces sample size, this was compared with random "
         "subsets of the unrestricted data at the same size, which give 0.348 with "
         "a 95 percent range of 0.252 to 0.428. The observed value lies below that "
         "range. Restricting the target side to loop behaves the same way, at 0.239 "
         "against a null of 0.372. Restricting it to strand or helix gives 0.264 "
         "and 0.235, both inside their nulls and so uninformative rather than "
         "contrary. Matching structural context therefore sharpens the difference "
         "where it discriminates at all, and never explains it. The same test run "
         "with the in-house assignment gives 0.185 for loop-to-loop against a null "
         "of 0.323, so the conclusion does not depend on which assignment is used.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Fold composition differs between the groups as expected, with 50 percent "
         "of antibody antigens in mainly-beta CATH classes against 22 percent of "
         "protein-protein partners. Matching on CATH class does not raise the "
         "correlation in any stratum. Resolution stratification, removal of the "
         "largest antigen clusters, leave-one-cluster-out resampling and variation "
         "of the hydrogen-bond distance criterion all leave the result intact. "
         "Finally, Fab and single-domain antibodies differ substantially in "
         "architecture [7, 18], and their cross-format grammar correlation of 0.18 falls "
         "inside the reliability ceiling of 0.25, with a range of 0.12 to 0.37, "
         "obtained by splitting the Fab corpus in half at the same sample size. "
         "That overlap is consistent with the two formats sharing one grammar, but "
         "it does not establish it. Holding this comparison to the same "
         "equivalence standard applied to the null results in Section 3.2, a two "
         "one-sided test against a margin of 0.085 does not reject non-equivalence, "
         "and the interval on the difference between the cross-format correlation "
         "and the within-format ceiling runs from -0.30 to +0.08. With 96 "
         "independence units on the single-domain side the comparison is "
         "underpowered, and it should be read as failing to detect a difference "
         "rather than as demonstrating equivalence.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "One alternative explanation deserves direct treatment, because it would "
         "account for the residue-level divergence without any difference in "
         "recognition chemistry. Hypervariable loops and general protein surfaces "
         "differ in amino-acid composition, antibody paratopes being notably rich "
         "in tyrosine, serine and glycine, and a matrix of raw contact counts would "
         "reproduce that compositional difference rather than any preference. The "
         "design of the statistic forecloses this. Every cell reports the observed "
         "count divided by the count expected from the product of the two "
         "residues' own marginal abundances at that interface class, so a residue "
         "that is simply common contributes proportionally to both the observed and "
         "the expected term and cancels. The quantity measured is preference "
         "relative to availability, not abundance. The two-background analysis of "
         "paratope composition reported below makes the same point from the "
         "opposite direction, since it shows which apparent enrichments survive "
         "when the antibody's own surface composition is used as the reference "
         "and which do not.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "phase6c_loop_confound.png",
           "Fig. 5 Structural context of the two interface classes. (a) Cross-group "
           "grammar correlation computed within matched secondary-structure strata, "
           "with the number of antibody contacts in each stratum printed on the bar "
           "and the all-contacts value marked by the dashed line. (b) "
           "Secondary-structure composition of each interface side, unit-weighted, "
           "for the antibody paratope, the antibody epitope and the two partners of "
           "a general protein-protein complex.")

    heading(doc, "3.7 The paratope is specialised and the epitope is not", 2)
    para(doc,
         "Correspondence between our coordinates and IMGT numbering was established "
         "by exact alpha-carbon coordinate matching against IMGT-numbered [13] "
         "structures of the same entries rather than by sequence alignment, "
         "which is least reliable at exactly the insertion-heavy loops of "
         "interest. Median mapping coverage was 1.000. Because IMGT "
         "numbering is a gapped scheme, positions absent from a given domain "
         "simply have no residue to map and contribute nothing, while "
         "insertions carry insertion codes and are assigned to the region "
         "containing their parent position. Region boundaries follow Lefranc "
         "and colleagues, with CDR1 at positions 27 to 38, CDR2 at 56 to 65 "
         "and CDR3 at 105 to 117. Loop contribution percentages are "
         "therefore computed over residues actually present, not over "
         "nominal scheme positions, so a short loop is not penalised for the "
         "gaps in its numbering.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Mapping every paratope residue onto IMGT numbering recovers the expected "
         "architecture. In Fab complexes, CDR-H3 contributes 37 percent of "
         "complementarity-determining region contacts, more than any other loop, "
         "while CDR-L2 contributes 2 percent. In single-domain antibodies, which "
         "have no light chain, CDR-H3 contributes 66 percent. Framework residues "
         "contribute 21 to 22 percent of paratope contacts and this proportion is "
         "almost constant across formats. These observations agree with the "
         "established picture [16, 22, 26] and serve here as evidence that the "
         "pipeline behaves correctly.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Paratope composition, however, requires care in the choice of background. "
         "Scored against general protein-protein interfaces, the paratope is "
         "enriched in tyrosine, tryptophan and serine, which reproduces a "
         "long-standing observation. Scored against the antibody's own accessible "
         "surface, tyrosine and tryptophan remain strongly enriched at 2.70 and "
         "3.04, but serine, threonine and glycine reverse sign, at -0.69, -0.94 and "
         "-0.45. Their apparent enrichment reflects the composition of antibody "
         "surfaces rather than selective use at the binding site. A study using a "
         "single background cannot distinguish these cases.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Applying the same two-background logic to the epitope places a known "
         "observation on a scale. Kringelum and colleagues [11] reported that epitope "
         "composition does not deviate significantly from equally exposed antigen "
         "surface, and that epitope residues could not be related to relative "
         "accessibility or to a protrusion index. Our measurements agree. The value "
         "of measuring three interface types with one pipeline is that the absence "
         "of a difference becomes interpretable, because there is something to "
         "compare it with. The mean absolute compositional enrichment of an "
         "interface relative to its own protein's accessible surface is 0.97 for the "
         "paratope and 0.48 for a general protein-protein interface, but 0.25 for "
         "the epitope. An epitope therefore carries roughly half the chemical "
         "distinguishability of an ordinary co-evolved interface and a quarter that "
         "of a paratope. A null result becomes a positioned one.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The geometry behaves the same way. Comparing each interface with its own "
         "protein's surface across independence units, the paratope protrudes "
         "strongly, with a matched-pairs effect of 0.58 at p = 3 x 10^-25, while the "
         "epitope shows no detectable protrusion preference, with an effect of -0.12 "
         "that does not reach significance at p = 0.44. General protein-protein "
         "interfaces are significantly less protruding than their surrounding "
         "surface, at -0.39. Two points follow. The frequently repeated statement "
         "that epitopes occupy protruding regions holds relative to the whole "
         "protein but not relative to other exposed surface, which is consistent "
         "with the earlier failure to relate epitopes to a protrusion index. And the "
         "protrusion in an antibody complex belongs to the antibody. Epitope "
         "secondary structure is the one axis on which the epitope is not quite "
         "ordinary. Under DSSP assignment 48.8 percent of epitope contacts come "
         "from loop residues against 44.9 percent for protein-protein partner "
         "surfaces, a difference of 3.9 points. That is real, and it is an order "
         "of magnitude smaller than the paratope excess over the same baseline, "
         "which is 22 points, so the loop character of an antibody complex remains "
         "overwhelmingly a property of the antibody. An earlier version of this "
         "comparison used a per-chain assignment that misses inter-chain bridges "
         "and reported the two as indistinguishable; the residual difference "
         "emerges only once strands are fully recovered.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Two further analyses support the same conclusion from different "
         "directions. Sixteen proteins in the corpus are solved both in complex with an "
         "antibody and in complex with a natural partner. This is a small sample, "
         "and the smallest on which any statement in this paper rests, so it is "
         "offered as corroboration of the corpus-wide results rather than as "
         "independent evidence. Among those 16 the median Jaccard overlap between "
         "the epitope and the natural partner site is 0.03, and 38 percent of pairs "
         "share no contacted position at all. Antibodies bind a different face of "
         "the same protein. Corpus-wide, the fraction of interface residues lying "
         "within a ligand-binding site is 1.6 percent for antibodies and 10.7 "
         "percent for protein-protein interfaces. Natural partners engage the "
         "functional surface; antibodies do not.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "pymol_matched_pair.png",
           "Fig. 6 The same protein read by an antibody and by its natural "
           "partner, with both footprints painted onto one copy of the shared "
           "protein. (a) Angiogenin, among the higher-overlap pairs. (b) A "
           "low-overlap pair. Red marks residues contacted only by the antibody, "
           "blue only by the natural partner, and purple by both.")
    figure(doc, "add5_epitope_logic.png",
           "Fig. 7 Interface character scored against each protein's own accessible "
           "surface. (a) Residue enrichment of the antibody paratope, the antibody "
           "epitope and the general protein-protein interface, each relative to the "
           "surface of the protein that carries it. (b) Mean absolute enrichment for "
           "the same three interface types, which measures how chemically "
           "distinguishable each is from its own surface. (c) Rank-biserial effect "
           "sizes for relative accessibility and protrusion, comparing interface "
           "residues with non-interface surface residues.")

    heading(doc, "3.8 Interface geometry and separability", 2)
    para(doc,
         "Antibody interfaces are compact and flat. The radius of gyration of the "
         "contacted patch is 10.7 angstroms for paratopes and 16.3 for "
         "protein-protein interfaces, and planarity is 0.85 against 0.79. Epitopes "
         "show the same compactness, at 10.5 angstroms. Part of the difference in "
         "radius restates the difference in interface size, but planarity is "
         "comparatively independent of size and points the same way. Because the "
         "antibody corpus is entirely X-ray, and well-diffracting proteins might be "
         "smaller and more rigid than proteins generally, this result was tested "
         "against the possibility that it reflects crystallisability rather than "
         "antibody chemistry. The premise does not hold: antibody antigens in this "
         "corpus are larger than protein-protein partners, with medians of 169 and "
         "135 residues. The difference also survives stratification. Within every "
         "quartile of bound-protein length the paratope patch remains more compact "
         "and more planar, with planarity effects between 0.12 and 0.39, and the "
         "difference is likewise present when both groups are restricted to "
         "structures better than 2.0 angstroms. Six "
         "hypervariable loops presented at the end of a beta-sandwich form a small "
         "and relatively flat binding face, and the antigen surface they engage is "
         "correspondingly flat.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Considering all measured interface properties together clarifies how "
         "large the overall difference is. In a principal component projection of "
         "twenty descriptors, 93.8 percent of antibody interfaces fall within the "
         "95 percent ellipse of the protein-protein distribution, so the two classes "
         "overlap heavily along any single axis. A classifier trained on the same "
         "descriptors with cross-validation grouped by independence unit "
         "nevertheless separates them with an area under the curve of 0.875. The "
         "antibody signature is real but distributed across many properties rather "
         "than concentrated in one.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "add4_pca.png",
           "Fig. 8 Separability of the two interface classes in descriptor space. "
           "(a) Projection of all interfaces onto the first two principal components "
           "of the descriptor set, with class means marked by crosses. (b) The ten "
           "largest logistic-regression coefficients, where negative values indicate "
           "descriptors that push an interface towards the antibody class and "
           "positive values towards the protein-protein class.")

    heading(doc, "3.9 How large is the difference in practice", 2)
    para(doc,
         "A formally significant difference between two grammars does not "
         "necessarily have practical consequences, and it is worth measuring how far "
         "the difference extends. Each matrix was used as a threading potential to "
         "discriminate native interfaces from decoys in which residue identities "
         "were permuted among interface positions. A potential built from "
         "protein-protein interfaces discriminates antibody interfaces with an area "
         "under the curve of 0.833, which exceeds the 0.779 obtained with a "
         "potential built from antibody interfaces themselves under leave-one-unit-"
         "out validation. The comparison persists when the protein-protein potential "
         "is subsampled to the same number of independence units, giving 0.810.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "This result requires care, because it appears to conflict with the docking "
         "literature. Antibody-antigen complexes are known to be handled poorly by "
         "generic docking potentials [6], and asymmetric potentials trained specifically "
         "on antibody data substantially improve pose prediction [4]. Our measurement "
         "does not contradict that work, because the two tasks are different in an "
         "important way. Rigid-body docking requires a potential to select the "
         "correct relative orientation of two proteins from a large space of "
         "incorrect ones, so its performance depends jointly on pair chemistry and "
         "on geometric complementarity. Our decoy set holds geometry fixed and "
         "permutes only residue identities, so it isolates pair chemistry alone. A "
         "potential can therefore capture pair chemistry well enough to identify the "
         "native sequence assignment while still failing to identify the native "
         "pose. Read together, the two observations locate the deficiency of generic "
         "potentials on antibody complexes in the geometric rather than the "
         "chemical term.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The bound this places on our own claim is worth stating plainly. The "
         "chemistry shared between the two interface classes, which the class-level "
         "correlation of 0.84 quantifies, is sufficient for a general potential to "
         "remain effective at the chemical discrimination task. The residue-level "
         "divergence is real, mechanistically specific and robust to every confound "
         "tested, and it is not large enough to make general pair statistics "
         "inapplicable. We also note that the antibody potential is estimated from "
         "fewer contacts even at matched unit counts, so a contact-matched "
         "comparison would refine the estimate.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The practical reading for computational tooling follows from the same "
         "distinction. Scoring functions used in docking and in complex-structure "
         "prediction combine a pairwise chemical term with terms for shape "
         "complementarity, desolvation and pose plausibility. Our result indicates "
         "that the pairwise chemical term learned from general interfaces is "
         "largely transferable to antibody complexes, while the established "
         "difficulty of generic potentials on antibody targets is better attributed "
         "to the geometric terms. That is consistent with what is known about "
         "antibody interfaces from this study: they are smaller, flatter and more "
         "compact than co-evolved interfaces, and they engage a surface patch "
         "carrying no complementary adaptation, which is a harder geometric "
         "recognition problem rather than a different chemical one. For tool "
         "development the implication is that effort spent on antibody-specific "
         "pairwise potentials may return less than effort spent on antibody-"
         "specific treatment of interface geometry and loop conformation. We offer "
         "this as a directed hypothesis rather than a demonstrated result, since "
         "testing it requires a docking benchmark that this study does not "
         "include.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "3.10 If epitopes are ordinary surface, why are they not random?", 2)
    para(doc,
         "The preceding sections establish that an epitope is not chemically, "
         "structurally or functionally distinguished from the rest of an antigen's "
         "exposed surface. That raises an obvious objection. If nothing marks a "
         "patch as an epitope, antibody binding should be distributed across the "
         "accessible surface, and it plainly is not. Resolving this requires a case "
         "in which many independent antibodies have been raised against one antigen, "
         "so that targeting can be observed rather than inferred. The corpus "
         "provides one: 42 independent complexes of lysozyme C, which the weighting "
         "scheme deliberately compresses to a single independence unit everywhere "
         "else in this study.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Targeting is strongly uneven. Sixty-six percent of the sequence is "
         "contacted by at least one antibody, so most of the surface is accessible "
         "to the repertoire, yet the distribution of contact frequency has a Gini "
         "coefficient of 0.49. Glycine 102 is contacted by 98 percent of the 42 "
         "antibodies, and asparagine 103, asparagine 19, tryptophan 63 and arginine "
         "73 by 74 to 83 percent. Rendered on the surface, the frequently contacted "
         "positions form one contiguous face rather than a scattered set of exposed "
         "residues, while the opposite face is almost untouched. The face is "
         "discontinuous in sequence and contiguous in space, which is the classical "
         "description of a conformational determinant [2], and its concentration of "
         "binding onto a few positions parallels the hot-spot organisation of general "
         "protein interfaces [3].",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The two observations are compatible, and their reconciliation is the point "
         "of this section. It rests on one antigen and 42 structures that the "
         "weighting scheme treats as a single independence unit everywhere else, so "
         "what follows is an interpretation consistent with the corpus-wide results "
         "rather than an independent test of them, and it should not be generalised "
         "to antigens with different surface topography. Convergence does not require that the targeted patch be "
         "chemically special. A contiguous, accessible, conformationally coherent "
         "face is sufficient, because that is what a set of hypervariable loops can "
         "engage. Selection acts on the antibody, which must present a complementary "
         "surface, not on the antigen, which need only be exposed. Immunodominance is "
         "therefore better read as a constraint of the paratope than as a property of "
         "the epitope, which is consistent with every asymmetry reported above.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "Two further details support this reading. The dominant face lies at the rim "
         "of the active-site cleft rather than on the catalytic residues. Glutamate "
         "35 and aspartate 52, the catalytic pair, are contacted by 12 and 19 percent "
         "of antibodies, whereas aspartate 101, which contributes to substrate "
         "binding and sits at the cleft edge, is contacted by 69 percent. Antibodies "
         "converge near the functional site without engaging it, which matches the "
         "corpus-wide finding that only 1.6 percent of antibody interface residues "
         "lie within a ligand-binding site. Second, alignment of 236 UniRef50 "
         "homologues indicates that frequently targeted positions are less conserved "
         "than rarely targeted ones, with mean identities of 0.47 and 0.62. This "
         "comparison does not reach significance and is reported as a trend, but its "
         "direction is consistent with targeting of surface that is under no "
         "selective pressure to be bound.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "pymol_immunodominance_montage.png",
           "Fig. 9 Contact frequency of 42 independently raised antibodies mapped "
           "onto the lysozyme C surface, shown from (a) the dominant face and (b) "
           "the same molecule rotated by 180 degrees. Pale yellow marks rarely "
           "contacted positions and orange frequently contacted ones; catalytic "
           "residues are coloured blue.")

    para(doc,
         "The same series answers a further question that the corpus-wide "
         "thermodynamics cannot. If 42 antibodies converge on one face, do they "
         "also converge on one thermodynamic strategy for engaging it? They do "
         "not. Across the series the heat-capacity change per 1000 square "
         "angstroms buried has a median of -39 cal per mol per K but ranges from "
         "-10 to -178, an eighteen-fold spread, and the coefficient of variation "
         "of 0.70 is only modestly below the 0.90 seen across antibody interfaces "
         "generally. Neither the number of immunodominant positions contacted nor "
         "the total buried surface predicts the value, with Spearman "
         "correlations of 0.15 and -0.02. What does stay narrow is the apolar "
         "fraction of buried surface, between 0.38 and 0.62 with a median of 0.42, "
         "which is lower than the antibody corpus median of 0.50 and indicates "
         "that this particular antigen face is engaged with unusually little "
         "hydrophobic burial. Spatial convergence therefore does not imply "
         "thermodynamic convergence. One surface patch admits many "
         "energetically distinct solutions, which is what one expects if the "
         "patch is selected for accessibility and contiguity rather than for a "
         "particular chemistry, and it is consistent with the finding that the "
         "epitope carries little compositional signature of its own.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    figure(doc, "case_thermodynamics.png",
           "Fig. 10 Thermodynamic signature across 42 independent antibodies against "
           "one antigen. (a) Heat-capacity change per unit buried surface for the "
           "lysozyme series against all antibody interfaces. (b) The same quantity "
           "against the number of immunodominant positions contacted. (c) Apolar "
           "fraction of buried surface against total buried surface.")
    heading(doc, "3.11 Synthesis", 2)
    para(doc,
         "The results assemble into a consistent account. An antibody is a "
         "specialised binding apparatus: it presents a small, flat, protruding, "
         "aromatic-rich surface that is strongly distinguishable in composition from "
         "the rest of the antibody, and it reads the antigen main chain by rules "
         "that co-evolved interfaces do not use, principally through tyrosine "
         "hydroxyl hydrogen bonds and cation-pi contacts. Because it cannot "
         "co-evolve with its target, it purchases affinity through directional "
         "chemistry rather than through hydrophobic burial, and its thermodynamic "
         "signature reflects this.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)
    para(doc,
         "The epitope, by contrast, is ordinary. It is chemically almost "
         "indistinguishable from the rest of the antigen surface, it does not "
         "protrude, its secondary structure matches that of other interfaces, it "
         "lies on a different face from the natural binding partner, and it avoids "
         "functional sites. The asymmetry between the two sides of an "
         "antibody-antigen interface is therefore not merely a matter of which side "
         "was selected. It is measurable in composition, geometry and chemistry, and "
         "it is the reason the interface as a whole differs from a co-evolved one.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------------------------------------------------------------- limitations
    heading(doc, "4. Limitations", 1)
    para(doc,
         "The antibody corpus consists entirely of X-ray structures and therefore "
         "excludes the cryo-electron microscopy literature, which is where most "
         "large viral glycoprotein complexes are found. Structural analyses on the "
         "protein-protein group used a random sample of 1,000 of the 2,292 "
         "complexes, which is reported rather than implied. Shape complementarity "
         "was computed on an approximation to the molecular surface that omits "
         "re-entrant regions, so its absolute values are not comparable with "
         "published figures although the between-group comparison is internally "
         "valid. CATH assignments cover 43 percent of antibody antigens against 80 "
         "percent of protein-protein partners, so the fold-matched control is weaker "
         "than the secondary-structure control. The secondary-structure "
         "implementation operates on one chain at a time, so beta bridges formed "
         "between chains are not detected and their residues are assigned as loop. "
         "This inflates loop content and deflates strand content in both groups. "
         "Because the same procedure is applied to each, the comparative claims are "
         "unaffected, but the absolute loop fractions in Section 3.6 should be read "
         "as upper bounds. Comparison against DSSP annotations computed on the "
         "complete structure quantifies the effect: strand recall is 0.950, so "
         "approximately five percent of strand residues are assigned to loop "
         "here, while helix recall is 0.998 and the reverse error is negligible "
         "at three residues in 65,336. Because the same procedure is applied to "
         "both groups, this shifts absolute loop fractions upward in each without "
         "creating a difference between them. Per-residue evolutionary "
         "conservation could not be computed corpus-wide in this environment and is "
         "reported only for the lysozyme case study. Interface-versus-surface "
         "comparisons in Section 3.7 are computed on unit-level medians, paired within "
         "independence unit. Finally, rigorous Poisson-Boltzmann electrostatics, molecular "
         "dynamics and docking were outside the scope of the analysis environment.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # ---------------------------------------------------------------- conclusions
    heading(doc, "5. Conclusions", 1)
    para(doc,
         "Antibodies and co-evolved protein partners share the chemistry of protein "
         "association at the level of chemical classes, and differ measurably at the "
         "level of residue-pair preference. The difference is not explained by loop "
         "context, fold composition, resolution, redundancy, contact criteria or "
         "antibody format. It is concentrated in the reading of the polypeptide "
         "backbone and is carried principally by tyrosine hydroxyl hydrogen bonds to "
         "the antigen main chain. Antibody interfaces are smaller, flatter and less "
         "hydrophobically driven than co-evolved interfaces, while the epitopes they "
         "engage are ordinary protein surface by every measure we applied. The "
         "specialisation in an antibody-antigen complex belongs to the antibody, not "
         "to the site it binds.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    # -------------------------------------------------------------- declarations
    heading(doc, "Declarations", 1)
    para(doc, "Funding: this research received no external funding.", size=10)
    para(doc, "Conflicts of interest: the authors declare none.", size=10)
    para(doc, "Data availability: all structures analysed are publicly available "
              "from the Protein Data Bank. Derived contact tables, per-cell "
              "matrices, structural feature tables and analysis outputs are provided "
              "as supplementary data.", size=10)
    para(doc, "Code availability: all analysis code is provided.", size=10)
    para(doc, "Author contributions: MR designed the study, implemented the analysis pipeline and drafted the manuscript. AS and AK contributed to study design and interpretation and revised the manuscript. All authors read and approved the final manuscript.", size=10)

    heading(doc, "References", 1)
    refs = [
        "Barlow D, Thornton J. Ion-pairs in proteins. Journal of Molecular Biology. 1983;168(4):867-885. doi:10.1016/s0022-2836(83)80079-5.",
        "Barlow DJ, Edwards MS, Thornton JM. Continuous and discontinuous protein antigenic determinants. Nature. 1986;322(6081):747-748. doi:10.1038/322747a0.",
        "Bogan AA, Thorn KS. Anatomy of hot spots in protein interfaces. Journal of Molecular Biology. 1998;280(1):1-9. doi:10.1006/jmbi.1998.1843.",
        "Brenke R, Hall DR, Chuang G, Comeau SR, Bohnuud T, Beglov D, Schueler-Furman O, Vajda S, Kozakov D. Application of asymmetric statistical potentials to antibody–protein docking. Bioinformatics. 2012;28(20):2608-2614. doi:10.1093/bioinformatics/bts493.",
        "Chakrabarti P, Janin J. Dissecting protein–protein recognition sites. Proteins: Structure, Function, and Bioinformatics. 2002;47(3):334-343. doi:10.1002/prot.10085.",
        "Chuang G, Kozakov D, Brenke R, Comeau SR, Vajda S. DARS (Decoys As the Reference State) Potentials for Protein-Protein Docking. Biophysical Journal. 2008;95(9):4217-4227. doi:10.1529/biophysj.108.135814.",
        "De Genst E, Silence K, Decanniere K, Conrath K, Loris R, Kinne J, Muyldermans S, Wyns L. Molecular basis for the preferential cleft recognition by dromedary heavy-chain antibodies. Proceedings of the National Academy of Sciences. 2006;103(12):4586-4591. doi:10.1073/pnas.0505379103.",
        "Gallivan JP, Dougherty DA. Cation-π interactions in structural biology. Proceedings of the National Academy of Sciences. 1999;96(17):9459-9464. doi:10.1073/pnas.96.17.9459.",
        "Jones S, Thornton JM. Principles of protein-protein interactions. Proceedings of the National Academy of Sciences. 1996;93(1):13-20. doi:10.1073/pnas.93.1.13.",
        "Kabsch W, Sander C. Dictionary of protein secondary structure: Pattern recognition of hydrogen‐bonded and geometrical features. Biopolymers. 1983;22(12):2577-2637. doi:10.1002/bip.360221211.",
        "Kringelum JV, Nielsen M, Padkjær SB, Lund O. Structural analysis of B-cell epitopes in antibody:protein complexes. Molecular Immunology. 2013;53(1-2):24-34. doi:10.1016/j.molimm.2012.06.001.",
        "Lawrence MC, Colman PM. Shape Complementarity at Protein/Protein Interfaces. Journal of Molecular Biology. 1993;234(4):946-950. doi:10.1006/jmbi.1993.1648.",
        "Lefranc M, Pommié C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental &amp; Comparative Immunology. 2003;27(1):55-77. doi:10.1016/s0145-305x(02)00039-3.",
        "Conte LL, Chothia C, Janin J. The atomic structure of protein-protein recognition sites 1 1Edited by A. R. Fersht. Journal of Molecular Biology. 1999;285(5):2177-2198. doi:10.1006/jmbi.1998.2439.",
        "MacCallum RM, Martin AC, Thornton JM. Antibody-antigen Interactions: Contact Analysis and Binding Site Topography. Journal of Molecular Biology. 1996;262(5):732-745. doi:10.1006/jmbi.1996.0548.",
        "Mian I, Bradwell AR, Olson AJ. Structure, function and properties of antibody binding sites. Journal of Molecular Biology. 1991;217(1):133-151. doi:10.1016/0022-2836(91)90617-f.",
        "Murphy KP, Freire E. Thermodynamics of structural stability and cooperative folding behavior in proteins. Advances in Protein Chemistry. 1992;43:313-361. doi:10.1016/s0065-3233(08)60556-2.",
        "Muyldermans S. Nanobodies: Natural Single-Domain Antibodies. Annual Review of Biochemistry. 2013;82(1):775-797. doi:10.1146/annurev-biochem-063011-092449.",
        "Peng H, Lee KH, Jian J, Yang A. Origins of specificity and affinity in antibody–protein interactions. Proceedings of the National Academy of Sciences. 2014;111(26). doi:10.1073/pnas.1401131111.",
        "Ramaraj T, Angel T, Dratz EA, Jesaitis AJ, Mumey B. Antigen–antibody interface properties: Composition, residue interactions, and features of 53 non-redundant structures. Biochimica et Biophysica Acta (BBA) - Proteins and Proteomics. 2012;1824(3):520-532. doi:10.1016/j.bbapap.2011.12.007.",
        "Rubinstein ND, Mayrose I, Halperin D, Yekutieli D, Gershoni JM, Pupko T. Computational characterization of B-cell epitopes. Molecular Immunology. 2008;45(12):3477-3489. doi:10.1016/j.molimm.2007.10.016.",
        "Sela-Culang I, Kunik V, Ofran Y. The Structural Basis of Antibody-Antigen Recognition. Frontiers in Immunology. 2013;4. doi:10.3389/fimmu.2013.00302.",
        "Shrake A, Rupley J. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. Journal of Molecular Biology. 1973;79(2):351-371. doi:10.1016/0022-2836(73)90011-9.",
        "Spolar RS, Record MT. Coupling of Local Folding to Site-Specific Binding of Proteins to DNA. Science. 1994;263(5148):777-784. doi:10.1126/science.8303294.",
        "Sundberg EJ, Mariuzza RA. Molecular recognition in antibody-antigen complexes. Advances in Protein Chemistry. 2002;119-160. doi:10.1016/s0065-3233(02)61004-6.",
        "Xu JL, Davis MM. Diversity in the CDR3 Region of VH Is Sufficient for Most Antibody Specificities. Immunity. 2000;13(1):37-45. doi:10.1016/s1074-7613(00)00006-6.",
    ]
    for i_, r_ in enumerate(refs, 1):
        para(doc, f"{i_}. {r_}", size=9, space=2)

    out = C.ROOT / "MANUSCRIPT.docx"
    doc.save(str(out))
    return out


# ------------------------------------------------------------------ supplementary
def build_supplementary():
    doc = new_doc()
    p = doc.add_paragraph()
    r = p.add_run("Supplementary Information")
    r.bold = True
    r.font.size = Pt(15)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    para(doc, "Antibodies read protein sequence conventionally and protein shape "
              "unconventionally: a matched comparison of antibody and co-evolved "
              "protein interfaces",
         size=11, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)

    heading(doc, "Contents", 1)
    table(doc, [
        ["Section", "Contents", "Items"],
        ["S1", "Corpus composition and redundancy structure", "Table S1"],
        ["S2", "Antigen length floor sensitivity", "Table S2"],
        ["S3", "Interface physicochemistry, full comparison", "Table S3"],
        ["S4", "Confound controls with size-matched nulls", "Table S4"],
        ["S5", "Paratope composition against two backgrounds", "Table S5"],
        ["S6", "Burial thermodynamics", "Table S6"],
        ["S7", "Interface geometry", "Table S7"],
        ["S8", "Verification checks and corrected errors", "Text"],
        ["S9", "Analysis code and the purpose of each script",
         "Tables S8 to S10"],
        ["S10", "Supplementary figures", "Figs. S1 to S20"],
        ["S11", "Supplementary data files", "Table S11"],
    ])
    para(doc, "Table of contents for the supplementary information.", size=9)

    heading(doc, "S1. Corpus composition", 1)
    table(doc, [
        ["Property", "Antibody-antigen", "Protein-protein"],
        ["Structures analysed", "512", "2,292"],
        ["cluster30 independence units", "233", "1,326"],
        ["Residue-pair contacts", "36,810", "406,397"],
        ["Atom-atom contacts mapped", "182,442", "1,767,916"],
        ["Median resolution (A)", "1.83", "1.80"],
        ["Experimental method", "100% X-ray", "predominantly X-ray"],
        ["Largest cluster share", "8.8%", "3.6%"],
        ["Singleton units", "149 of 233", "719 of 1,326"],
    ])
    para(doc, "Table S1 Corpus composition and redundancy structure.", size=9)

    heading(doc, "S2. Antigen length floor sensitivity", 1)
    table(doc, [
        ["Minimum antigen length (residues)", "Structures retained"],
        ["0", "774"], ["30", "530"], ["50 (adopted)", "512"], ["80", "481"],
    ])
    para(doc, "Table S2 Sensitivity of the antibody corpus to the antigen length "
              "floor.", size=9)

    heading(doc, "S3. Interface physicochemistry, full comparison", 1)
    inv = J("inventory.json")
    rows = [["Property", "Antibody median", "PPI median", "Effect", "p"]]
    labels = {
        "n_pairs": "Residue pairs per interface",
        "n_probe_residues": "Probe-side interface residues",
        "n_target_residues": "Target-side interface residues",
        "n_atom_contacts": "Atom contacts per interface",
        "mean_min_dist": "Mean closest approach (A)",
        "f_hbond": "Hydrogen bond fraction",
        "f_salt": "Salt bridge fraction",
        "f_vdw": "van der Waals fraction",
        "f_stack": "Aromatic stacking fraction",
        "f_opp_charge": "Opposite-charge pair fraction",
        "f_like_charge": "Like-charge pair fraction",
        "f_apolar_apolar": "Apolar-apolar fraction",
        "f_aromatic": "Aromatic involvement",
        "f_cation_pi": "Cation-pi fraction",
        "f_target_sidechain_read": "Read via target side chain",
    }
    allstats = {**inv.get("size", {}), **inv.get("physicochemistry", {})}
    for k, lab in labels.items():
        s = allstats.get(k)
        if not s or s.get("p") is None:
            continue
        rows.append([lab, f"{s['median_antibody']:.3f}", f"{s['median_ppi']:.3f}",
                     f"{s['effect_rank_biserial']:+.3f}", f"{s['p']:.2e}"])
    table(doc, rows)
    para(doc, "Table S3 Interface properties compared across cluster30 units by "
              "Mann-Whitney U test. Effect is the rank-biserial correlation.", size=9)

    heading(doc, "S4. Confound controls", 1)
    table(doc, [
        ["Confound", "Test", "Observed r", "Reference", "Verdict"],
        ["Secondary structure", "loop to loop", "0.185", "null 0.323 [0.247, 0.390]",
         "below null"],
        ["Secondary structure", "target strand", "0.134", "null 0.275 [0.177, 0.366]",
         "below null"],
        ["Fold class", "mainly alpha", "0.167", "null 0.242 [0.168, 0.336]",
         "not raised"],
        ["Fold class", "mainly beta", "0.260", "null 0.312 [0.214, 0.394]",
         "not raised"],
        ["Fold class", "alpha/beta", "0.172", "null 0.279 [0.175, 0.353]",
         "not raised"],
        ["Resolution", "better than 2.0 A", "0.356", "null 0.401 [0.373, 0.453]",
         "not raised"],
        ["Resolution", "2.0 A or worse", "0.230", "null 0.264 [0.171, 0.348]",
         "not raised"],
        ["Redundancy", "drop 5 largest clusters", "0.449", "baseline 0.406",
         "unchanged"],
        ["Redundancy", "leave one cluster out", "0.382 to 0.457", "baseline 0.406",
         "unchanged"],
        ["H-bond criterion", "3.0 / 3.2 / 3.5 A", "0.300 / 0.317 / 0.276", "-",
         "stable"],
        ["Antibody format", "Fab vs single domain", "0.179",
         "ceiling 0.253 [0.117, 0.370]", "within ceiling"],
    ])
    para(doc, "Table S4 Confound controls. Each stratified comparison is referenced "
              "to a null drawn at the same sample size, because restriction reduces "
              "sample size and attenuates correlation independently of biology.",
         size=9)

    heading(doc, "S5. Paratope composition against two backgrounds", 1)
    pj = J("paratope.json").get("paratope_composition", {})
    rows = [["Residue", "log2 vs PPI interface", "log2 vs antibody own surface"]]
    for a in C.AAS:
        d = pj.get(a)
        if d:
            rows.append([a, f"{d['vs_ppi_interface']:+.2f}",
                         f"{d['vs_antibody_surface']:+.2f}"])
    table(doc, rows)
    para(doc, "Table S5 Paratope residue enrichment scored against two backgrounds. "
              "Residues whose sign differs between the two columns are enriched "
              "relative to other interfaces but not relative to what the antibody "
              "itself presents.", size=9)

    heading(doc, "S6. Burial thermodynamics", 1)
    th = J("thermodynamics.json").get("stats", {})
    rows = [["Quantity", "Antibody", "Protein-protein", "Effect", "p"]]
    for k, lab in (("dA_total", "Total buried surface (A^2)"),
                   ("frac_apolar_buried", "Apolar fraction of buried surface"),
                   ("dCp_cal_mol_K", "dCp (cal/mol/K)"),
                   ("dCp_per_1000A2", "dCp per 1000 A^2 (cal/mol/K)")):
        s = th.get(k)
        if s:
            rows.append([lab, f"{s['median_antibody']:.2f}",
                         f"{s['median_ppi']:.2f}",
                         f"{s['effect_rank_biserial']:+.3f}", f"{s['p']:.2e}"])
    table(doc, rows)
    para(doc, "Table S6 Buried surface and heat-capacity change computed from "
              "atom-level solvent accessibility.", size=9)

    heading(doc, "S7. Interface geometry", 1)
    g = J("pca_geometry.json").get("geometry_medians", {})
    rows = [["Quantity", "Antibody", "Protein-protein", "Effect", "p"]]
    for k, lab in (("probe_rg", "Probe patch radius of gyration (A)"),
                   ("probe_planarity", "Probe patch planarity"),
                   ("probe_n_spatial_patches", "Probe spatial sub-patches"),
                   ("target_rg", "Target patch radius of gyration (A)"),
                   ("target_planarity", "Target patch planarity"),
                   ("target_n_spatial_patches", "Target spatial sub-patches")):
        s = g.get(k)
        if s:
            rows.append([lab, f"{s['antibody']:.2f}", f"{s['ppi']:.2f}",
                         f"{s['effect_rank_biserial']:+.3f}", f"{s['p']:.2e}"])
    table(doc, rows)
    para(doc, "Table S7 Geometric properties of the contacted patch.", size=9)

    heading(doc, "S8. Verification checks", 1)
    para(doc,
         "Thirty-nine automated checks re-derive reported quantities by independent "
         "routes. Five errors were identified and corrected during development, all "
         "of which produced superficially plausible output: duplicated "
         "protein-protein contacts inflating Fisher sample size; a dihedral sign "
         "convention displacing every side-chain rotamer by 60 degrees, which passed "
         "a Ramachandran sanity check because it preserved the sign of phi; "
         "asparagine N-glycosylation sites recorded as post-translational "
         "modifications through MODRES records; an omitted sign in the heat-capacity "
         "relation that inverted the thermodynamic interpretation; and a direction "
         "mismatch when matching partner identity across the two arms, which "
         "compared each protein against its partner rather than against itself.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "S8b. Secondary-structure cross-validation", 1)
    v = J("dssp_validation.json")
    table(doc, [
        ["Quantity", "Value"],
        ["Reference implementation", "pydssp 0.9.1, numpy backend"],
        ["Chains compared", f"{v.get('n_chains_compared', 0):,}"],
        ["Residues compared", f"{v.get('n_residues_compared', 0):,}"],
        ["Three-state agreement", f"{v.get('three_state_agreement', 0):.4f}"],
        ["Cohen kappa", f"{v.get('cohens_kappa', 0):.4f}"],
        ["Failures", str(v.get("n_failed", 0))],
        ["Helix recall / precision",
         f"{v['per_state']['helix']['recall_vs_pydssp']:.3f} / "
         f"{v['per_state']['helix']['precision_vs_pydssp']:.3f}"],
        ["Strand recall / precision",
         f"{v['per_state']['strand']['recall_vs_pydssp']:.3f} / "
         f"{v['per_state']['strand']['precision_vs_pydssp']:.3f}"],
        ["Loop recall / precision",
         f"{v['per_state']['loop']['recall_vs_pydssp']:.3f} / "
         f"{v['per_state']['loop']['precision_vs_pydssp']:.3f}"],
    ])
    para(doc, "Table S8b Agreement between the in-house Kabsch-Sander implementation "
              "and an independent implementation of the same algorithm. Both assign "
              "from a single ordered backbone array, so this establishes fidelity to "
              "the algorithm and does not address inter-chain beta bridges.", size=9)

    heading(doc, "S8c. Validation against the DSSP program", 1)
    w = J("dssp_real_validation.json")
    table(doc, [
        ["Quantity", "Value"],
        ["Reference", "DSSP, whole structure, via PDB-REDO"],
        ["Structures compared", str(w.get("n_structures_compared", 0))],
        ["Residues compared", f"{w.get('n_residues_compared', 0):,}"],
        ["Three-state agreement", f"{w.get('three_state_agreement', 0):.4f}"],
        ["Cohen kappa", f"{w.get('cohens_kappa', 0):.4f}"],
        ["Helix recall", f"{w.get('helix_recall', 0):.4f}"],
        ["Strand recall", f"{w.get('strand_recall', 0):.4f}"],
        ["DSSP strand assigned loop here",
         f"{w.get('strand_called_loop_by_us', 0):,} of "
         f"{w.get('strand_total_in_dssp', 0):,}"],
    ])
    para(doc, "Table S8c Agreement with the DSSP program itself. Residues were "
              "matched by alpha-carbon coordinate. The one-directional strand "
              "disagreement quantifies the per-chain limitation of the in-house "
              "implementation.", size=9)

    heading(doc, "S9. Analysis code", 1)
    para(doc,
         "All analyses were implemented as standalone scripts sharing a single "
         "statistical engine. Each script writes its results to a JSON file and its "
         "figures to disk, so every number reported in the manuscript can be traced "
         "to the script that produced it. Scripts are listed below in the order they "
         "are run.",
         align=WD_ALIGN_PARAGRAPH.JUSTIFY)

    heading(doc, "S9.1 Shared engine and corpus preparation", 2)
    table(doc, [
        ["Script", "Purpose"],
        ["aa_common.py",
         "Statistical engine shared by every analysis. Loads contact tables, "
         "assigns independence units and weights, builds 20 by 20 matrices, "
         "computes log odds against the marginal-product null, runs the "
         "unit bootstrap and the Fisher plus Benjamini-Hochberg procedure, and "
         "deduplicates the two-direction protein-protein tables where a true "
         "sample size is required."],
        ["phase0_provenance_audit.py",
         "Re-derives contacts from raw coordinates with an independently written "
         "extractor and compares them with the inherited tables. Also tests corpus "
         "overlap between the two groups and quantifies redundancy."],
        ["phase1_manifest.py",
         "Builds the 3,084-structure corpus manifest, joins Structural Antibody "
         "Database metadata, and records exclusion reasons for structures that do "
         "not enter the analysis."],
    ])
    para(doc, "Table S8 Shared engine and corpus preparation scripts.", size=9)

    heading(doc, "S9.2 Feature extraction and core analyses", 2)
    table(doc, [
        ["Script", "Purpose"],
        ["phase2_inventory.py",
         "Collapses contacts to one record per interface and computes bond-type "
         "composition, charge and aromatic content, interface size and packing."],
        ["aa_heatmap.py",
         "Builds the 20 by 20 recognition matrices for both groups, at both "
         "independence-unit definitions, with bootstrap stability and multiple "
         "testing control, plus the chemical-class collapse and the bond-type "
         "decomposition."],
        ["phase4_conservation.py",
         "Correlates the two grammars overall, per bond type and with cysteine "
         "masked, and classifies cells as universal, group-specific or opposed."],
        ["phase5_readout.py",
         "Splits every contact between target side chain and target main chain and "
         "rebuilds the matrices separately for the two readout channels."],
        ["phase6a_structural.py",
         "Single structural pass over 1,512 complexes. Implements Kabsch and Sander "
         "secondary structure, computes backbone and side-chain torsions, and "
         "computes solvent-accessible surface in complex and in isolation."],
        ["phase6a_fix_dihedrals.py",
         "Recomputes torsion angles after correction of the dihedral sign "
         "convention, asserting the routine against a reference implementation "
         "before writing."],
        ["phase6b_functional_groups.py",
         "Maps every atom-atom contact to one of fourteen chemical functional "
         "groups and builds the functional-group recognition maps."],
        ["phase6c_conformation.py",
         "Tests whether the antibody signature is loop context, by recomputing the "
         "cross-group correlation within matched secondary structure. Also produces "
         "Ramachandran, rotamer and burial distributions."],
        ["phase6d_folds.py",
         "Retrieves CATH and SCOP assignments and repeats the comparison within "
         "matched fold class."],
        ["phase7_paratope.py",
         "Transfers IMGT numbering onto the corpus by exact alpha-carbon coordinate "
         "matching, assigns complementarity-determining regions, and scores "
         "paratope composition against two backgrounds."],
        ["phase8_differential.py",
         "Permutation test of grammar identity, differential map with intervals on "
         "the difference, matrix asymmetry test and dynamic-range test."],
        ["phase9_modified.py",
         "Scans raw coordinates for post-translational modifications and glycans at "
         "interfaces, which the contact tables exclude by construction."],
        ["phase10_literature.py",
         "Assembles the concordance between measured quantities and published "
         "claims."],
        ["phase11_controls.py",
         "Resolution, redundancy, hydrogen-bond criterion and antibody-format "
         "controls."],
    ])
    para(doc, "Table S9 Feature extraction and core analysis scripts.", size=9)

    heading(doc, "S9.3 Extended analyses, case studies and verification", 2)
    table(doc, [
        ["Script", "Purpose"],
        ["add1_transfer.py",
         "Uses each matrix as a threading potential and measures native-versus-decoy "
         "discrimination with self and cross potentials under "
         "leave-one-unit-out validation."],
        ["add2_thermodynamics.py",
         "Computes buried apolar and polar surface at atom level and derives the "
         "heat-capacity change on association."],
        ["add3_shape_complementarity.py",
         "Computes shape complementarity from surface points and normals on the van "
         "der Waals contact surface."],
        ["add4_pca_geometry.py",
         "Computes patch radius of gyration, planarity and sub-patch counts, then "
         "projects all interface descriptors and measures class separability with "
         "unit-grouped cross-validation."],
        ["add5_epitope_logic.py",
         "Scores each interface against its own protein's accessible surface, and "
         "measures relative accessibility and protrusion."],
        ["add6_functional_sites.py",
         "Defines functional sites by proximity to bound biological ligands and "
         "measures interface overlap with them."],
        ["immunodominance.py",
         "Normalises antigen residue numbering across depositions and maps contact "
         "frequency across all antibodies raised against one antigen."],
        ["case_lysozyme.py",
         "Case study combining immunodominance, catalytic-site targeting and "
         "sequence conservation from aligned homologues."],
        ["case_matched_targets.py",
         "Identifies proteins solved both with an antibody and with a natural "
         "partner and measures overlap between the two interfaces."],
        ["pymol_mechanism.py",
         "Renders tyrosine hydroxyl hydrogen bonds to the antigen main chain in "
         "high-resolution complexes."],
        ["pymol_immunodominance.py",
         "Paints antibody contact frequency onto the antigen surface."],
        ["audit_phases.py",
         "Thirteen checks on the statistical machinery, including an independent "
         "reimplementation of the multiple-testing procedure and verification that "
         "the vectorised bootstrap matches explicit resampling."],
        ["audit_full.py",
         "Twenty-six checks across all analyses, each re-deriving a reported "
         "quantity by an independent route."],
        ["build_manuscript.py",
         "Generates the manuscript and this supplementary file, drawing numerical "
         "values from the result files so that text cannot drift from results."],
    ])
    para(doc, "Table S10 Extended analyses, case studies and verification scripts.",
         size=9)

    heading(doc, "S10. Supplementary figures", 1)
    figure(doc, "heatmap_chemical_classes.png",
           "Fig. S1 Chemical-class collapse of the recognition matrices. (a) Antibody-antigen. (b) General protein-protein. Cell values are log2(observed/expected).", width=5.9)
    figure(doc, "heatmap_bond_decomposition.png",
           "Fig. S2 Recognition matrices decomposed by bond type. (a) to (d) Antibody-antigen, showing salt bridge, hydrogen bond, van der Waals and aromatic stacking respectively. (e) to (h) The same four decompositions for general protein-protein interfaces.", width=5.9)
    figure(doc, "heatmap_20x20_uniprot.png",
           "Fig. S3 Recognition matrices computed with UniProt accession as the independence unit, as a sensitivity analysis on the cluster30 definition used throughout. (a) Antibody-antigen. (b) General protein-protein.", width=5.9)
    figure(doc, "phase4_conservation.png",
           "Fig. S4 Conservation of recognition chemistry between the two classes. (a) Cell-by-cell agreement, coloured by the chemical class of the probe residue. (b) Correlation by bond type, with and without cysteine. (c) Classification of cells as universal, class-specific or opposed.", width=5.9)
    figure(doc, "phase6b_functional_groups.png",
           "Fig. S5 Functional-group recognition maps over fourteen chemical groups. (a) Antibody-antigen. (b) General protein-protein. (c) Difference between them, with orange indicating antibody enrichment.", width=5.9)
    figure(doc, "phase6d_folds.png",
           "Fig. S6 Fold composition and the fold-matched control. (a) CATH class of the recognised protein in each group. (b) Cross-class grammar correlation within matched fold class, against size-matched nulls shown in grey.", width=5.9)
    figure(doc, "phase11_controls.png",
           "Fig. S7 Robustness controls. (a) Resolution stratification against size-matched nulls. (b) Effect of removing the largest antigen clusters. (c) Sensitivity to the hydrogen-bond distance criterion. (d) Fab against single-domain antibodies, with the within-format reliability ceiling.", width=5.9)
    figure(doc, "phase8_differential.png",
           "Fig. S8 Formal comparison of the two grammars. (a) Differential map, antibody minus protein-protein, with dots marking cells whose interval excludes zero. (b) Permutation test against the null of one shared grammar. (c) Matrix asymmetry for each class.", width=5.9)
    figure(doc, "phase7_paratope.png",
           "Fig. S9 Paratope architecture and epitope character. (a) Share of paratope contacts by CDR loop and framework, split by construct type. (b) Paratope residue enrichment scored against two backgrounds. (c) Number of sequence segments per epitope.", width=5.9)
    figure(doc, "phase9_modified.png",
           "Fig. S10 Modified residues and glycans at interfaces. (a) Post-translational modifications. (b) Glycans. Bars give the fraction of structures with at least one such group contacting the partner.", width=5.9)
    figure(doc, "add1_transfer.png",
           "Fig. S11 Native-versus-decoy discrimination using each matrix as a threading potential. (a) Antibody interfaces scored by the antibody potential under leave-one-unit-out validation and by the protein-protein potential. (b) The equivalent comparison for protein-protein interfaces.", width=5.9)
    figure(doc, "add2_thermodynamics.png",
           "Fig. S12 Burial thermodynamics from atom-level solvent accessibility. (a) Fraction of buried surface that is apolar. (b) Heat-capacity change per 1000 square angstroms buried. (c) Total buried surface. Distributions are over cluster30 units.", width=5.9)
    figure(doc, "add3_shape_complementarity.png",
           "Fig. S13 Shape complementarity of the two interface classes, computed on the van der Waals contact surface. Values are not directly comparable with published figures because re-entrant surface is omitted; the between-class comparison is internally valid.", width=5.9)
    figure(doc, "add4_geometry.png",
           "Fig. S14 Geometric spread of the contacted patch. (a) to (c) Probe side: radius of gyration, planarity, and number of spatially separate sub-patches. (d) to (f) The same three quantities for the target side.", width=5.9)
    figure(doc, "add6_functional_sites.png",
           "Fig. S15 Overlap between interfaces and ligand-binding sites. (a) Fraction of interface residues lying within a ligand-binding site. (b) Ligands most often defining those sites.", width=5.9)
    figure(doc, "phase6c_conformation.png",
           "Fig. S16 Backbone dihedrals, side-chain rotamers and burial. (a) to (c) Antibody-antigen: Ramachandran distribution of interface residues, chi1 for interface against background, and buried surface per residue. (d) to (f) The same three for general protein-protein interfaces.", width=5.9)
    figure(doc, "immunodominance_profile.png",
           "Fig. S17 Per-position contact frequency across 42 independently raised antibodies against lysozyme C, in canonical numbering. Labelled positions are the eight most frequently contacted.", width=5.9)
    figure(doc, "case_lysozyme.png",
           "Fig. S18 Lysozyme C case study. (a) Fraction of the 42 antibodies contacting each position, with catalytic residues marked by dashed lines. (b) Sequence conservation at each position from 236 aligned UniRef50 homologues.", width=5.9)
    figure(doc, "case_matched_targets.png",
           "Fig. S19 Proteins solved both with an antibody and with a natural partner. (a) Distribution of Jaccard overlap between the epitope and the partner site across 16 matched proteins. (b) The highest-overlap pairs individually.", width=5.9)
    figure(doc, "phase2_interface_inventory.png",
           "Fig. S20 Interface size and packing. (a) to (f) Residue pairs per interface, interface residues on each side, atom contacts, contacts per probe residue, and mean closest-approach distance. Distributions are over cluster30 units.", width=5.9)
    figure(doc, "add7_equivalence.png",
           "Fig. S21 Equivalence testing of the reported null results. Points give the difference of means between the classes, bars the bootstrap interval, and the shaded band the equivalence margin. Green marks quantities established as equivalent.", width=5.9)
    figure(doc, "add7_channel_null.png",
           "Fig. S22 Permutation null for readout-channel independence. (a) Antibody-antigen. (b) General protein-protein. The shaded band is the correlation expected if the moiety of a contact carried no grammar information.", width=5.9)
    figure(doc, "add7_crystallisability.png",
           "Fig. S23 Test of whether the interface geometry difference reflects crystallisability. (a) Radius of gyration and (b) planarity effects within quartiles of bound-protein length and within the high-resolution subset.", width=5.9)

    heading(doc, "S11. Supplementary data files", 1)
    table(doc, [
        ["File", "Contents"],
        ["corpus_manifest.csv", "All 3,084 structures with descriptions, chains, "
                               "resolution, units and exclusion reasons"],
        ["matrix_antibody_cluster30.csv", "Per-cell log odds, q values, bootstrap "
                                          "stability and counts"],
        ["matrix_ppi_cluster30.csv", "As above for the protein-protein group"],
        ["fgroup_matrix_antibody.csv", "Functional-group contact matrix"],
        ["fgroup_matrix_ppi.csv", "As above for the protein-protein group"],
        ["interface_inventory.csv", "Per-interface bond and chemistry composition"],
        ["structural_residues.csv", "1,019,571 residues with secondary structure, "
                                    "torsions and buried surface"],
        ["interface_geometry.csv", "Patch radius, planarity and sub-patch counts"],
        ["thermodynamics.csv", "Per-interface buried apolar and polar surface"],
        ["shape_complementarity.csv", "Per-interface shape complementarity"],
        ["paratope_regions.csv", "Every paratope contact with its IMGT region"],
        ["differential_cells.csv", "Per-cell differences with intervals"],
        ["modified_contacts.csv", "Modified residues and glycans at interfaces"],
    ])
    para(doc, "Table S11 Supplementary data files provided with this manuscript.",
         size=9)

    out = C.ROOT / "SUPPLEMENTARY.docx"
    doc.save(str(out))
    return out


if __name__ == "__main__":
    a = build_manuscript()
    b = build_supplementary()
    print("wrote", a)
    print("wrote", b)
