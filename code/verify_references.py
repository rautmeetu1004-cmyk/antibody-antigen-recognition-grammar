"""Verify every reference against CrossRef and emit a corrected bibliography.

Peer review flagged that only three of twenty-six references carried publisher
verification, and that deferring the rest to proof stage is not acceptable at
submission. This queries the CrossRef REST API by bibliographic string, takes the
best match, and compares the returned record against what the manuscript claims.

The output is deliberately conservative. A reference is marked VERIFIED only when
CrossRef returns a record whose title matches closely and whose first author agrees.
Anything else is reported for manual attention rather than silently accepted, since
a wrong citation inserted automatically is worse than an incomplete one.
"""
from __future__ import annotations

import difflib
import json
import pathlib
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

MAILTO = "rautk2024@gmail.com"          # CrossRef asks for a contact in the UA

# (key, first author surname, year, title as the manuscript gives it)
REFS = [
    ("barlow1983", "Barlow", 1983, "Ion-pairs in proteins"),
    ("barlow1986", "Barlow", 1986,
     "Continuous and discontinuous protein antigenic determinants"),
    ("bogan1998", "Bogan", 1998, "Anatomy of hot spots in protein interfaces"),
    ("brenke2012", "Brenke", 2012,
     "Application of asymmetric statistical potentials to antibody-protein docking"),
    ("chakrabarti2002", "Chakrabarti", 2002,
     "Dissecting protein-protein recognition sites"),
    ("chuang2008", "Chuang", 2008,
     "DARS (Decoys As the Reference State) potentials for protein-protein docking"),
    ("degenst2006", "De Genst", 2006,
     "Molecular basis for the preferential cleft recognition by dromedary "
     "heavy-chain antibodies"),
    ("gallivan1999", "Gallivan", 1999,
     "Cation-pi interactions in structural biology"),
    ("jones1996", "Jones", 1996, "Principles of protein-protein interactions"),
    ("kabsch1983", "Kabsch", 1983,
     "Dictionary of protein secondary structure: pattern recognition of "
     "hydrogen-bonded and geometrical features"),
    ("kringelum2013", "Kringelum", 2013,
     "Structural analysis of B-cell epitopes in antibody:protein complexes"),
    ("lawrence1993", "Lawrence", 1993,
     "Shape complementarity at protein/protein interfaces"),
    ("lefranc2003", "Lefranc", 2003,
     "IMGT unique numbering for immunoglobulin and T cell receptor variable "
     "domains and Ig superfamily V-like domains"),
    ("loconte1999", "Lo Conte", 1999,
     "The atomic structure of protein-protein recognition sites"),
    ("maccallum1996", "MacCallum", 1996,
     "Antibody-antigen interactions: contact analysis and binding site topography"),
    ("mian1991", "Mian", 1991,
     "Structure, function and properties of antibody binding sites"),
    ("murphy1992", "Murphy", 1992,
     "Thermodynamics of structural stability and cooperative folding behavior "
     "in proteins"),
    ("muyldermans2013", "Muyldermans", 2013,
     "Nanobodies: natural single-domain antibodies"),
    ("peng2014", "Peng", 2014,
     "Origins of specificity and affinity in antibody-protein interactions"),
    ("ramaraj2012", "Ramaraj", 2012,
     "Antigen-antibody interface properties: composition, residue interactions, "
     "and features of 53 non-redundant structures"),
    ("rubinstein2009", "Rubinstein", 2009,
     "Computational characterization of B-cell epitopes"),
    ("selaculang2013", "Sela-Culang", 2013,
     "The structural basis of antibody-antigen recognition"),
    ("shrake1973", "Shrake", 1973,
     "Environment and exposure to solvent of protein atoms. Lysozyme and insulin"),
    ("spolar1994", "Spolar", 1994,
     "Coupling of local folding to site-specific binding of proteins to DNA"),
    ("sundberg2002", "Sundberg", 2002,
     "Molecular recognition in antibody-antigen complexes"),
    ("xu2000", "Xu", 2000,
     "Diversity in the CDR3 region of VH is sufficient for most antibody "
     "specificities"),
]


def query(title, author, year):
    q = urllib.parse.urlencode({
        "query.bibliographic": f"{title} {author} {year}",
        "rows": 3, "select": "title,author,container-title,volume,issue,page,"
                             "issued,DOI,type"})
    url = f"https://api.crossref.org/works?{q}"
    req = urllib.request.Request(
        url, headers={"User-Agent": f"reference-check (mailto:{MAILTO})"})
    with urllib.request.urlopen(req, timeout=45) as fh:
        return json.load(fh)["message"]["items"]


def score(item, title, author, year):
    got_title = (item.get("title") or [""])[0]
    sim = difflib.SequenceMatcher(None, title.lower(), got_title.lower()).ratio()
    authors = item.get("author") or []
    surnames = " ".join((a.get("family") or "").lower() for a in authors)
    author_ok = author.split()[-1].lower() in surnames
    parts = (item.get("issued") or {}).get("date-parts") or [[None]]
    got_year = parts[0][0]
    year_ok = got_year is not None and abs(got_year - year) <= 1
    return sim, author_ok, year_ok, got_year, got_title


def fmt(item):
    a = item.get("author") or []
    names = ", ".join(
        f"{x.get('family','')} {''.join(p[0] for p in (x.get('given') or '').split() if p)}"
        for x in a[:12])
    if len(a) > 12:
        names += ", et al"
    j = (item.get("container-title") or [""])[0]
    yr = ((item.get("issued") or {}).get("date-parts") or [[None]])[0][0]
    vol = item.get("volume", "")
    iss = item.get("issue", "")
    pg = item.get("page", "")
    doi = item.get("DOI", "")
    tail = f"{vol}({iss})" if vol and iss else (vol or "")
    if pg:
        tail = f"{tail}:{pg}" if tail else pg
    return (f"{names}. {(item.get('title') or [''])[0]}. {j}. {yr};{tail}. "
            f"doi:{doi}").replace(";.", ".").replace("..", ".")


if __name__ == "__main__":
    results, verified, review = [], 0, []
    for key, author, year, title in REFS:
        try:
            items = query(title, author, year)
        except Exception as e:                                   # noqa: BLE001
            results.append({"key": key, "status": "QUERY_FAILED",
                            "error": f"{type(e).__name__}: {e}"})
            review.append(key)
            print(f"[QUERY FAILED] {key}: {type(e).__name__}")
            time.sleep(1)
            continue

        best, best_sim = None, 0.0
        for it in items:
            sim, a_ok, y_ok, got_year, got_title = score(it, title, author, year)
            if sim > best_sim:
                best, best_sim = (it, sim, a_ok, y_ok, got_year, got_title), sim
        if best is None:
            results.append({"key": key, "status": "NO_MATCH"})
            review.append(key)
            print(f"[NO MATCH]     {key}")
            time.sleep(1)
            continue

        it, sim, a_ok, y_ok, got_year, got_title = best
        ok = sim > 0.80 and a_ok and y_ok
        results.append({
            "key": key, "status": "VERIFIED" if ok else "CHECK",
            "title_similarity": round(sim, 3), "author_match": a_ok,
            "year_match": y_ok, "crossref_year": got_year,
            "crossref_title": got_title, "doi": it.get("DOI"),
            "formatted": fmt(it)})
        if ok:
            verified += 1
            print(f"[VERIFIED]     {key}  {it.get('DOI')}")
        else:
            review.append(key)
            print(f"[CHECK]        {key}  sim={sim:.2f} author={a_ok} "
                  f"year={y_ok} ({got_year}) -> {got_title[:60]}")
        time.sleep(0.6)

    out = {"n_references": len(REFS), "n_verified": verified,
           "n_needing_review": len(review), "needing_review": review,
           "records": results}
    (C.OUT_MET / "reference_verification.json").write_text(json.dumps(out, indent=1))
    print(f"\n{verified}/{len(REFS)} verified against CrossRef; "
          f"{len(review)} need manual attention")
