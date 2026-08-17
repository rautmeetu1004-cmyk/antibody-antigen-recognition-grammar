"""Phase 1 - the corpus manifest: every PDB in the study and what it actually is.

Writes one row per structure across both groups, so the corpus can be read, sorted
and audited by a human without touching the analysis code.

  metadata/corpus_manifest.csv    combined, both groups, common schema
  metadata/analysis/corpus.json   composition roll-up

Antibody rows are enriched from the SAbDab summary (metadata/sabdab_summary.csv,
21,914 rows pulled from the live OPIG API). SAbDab uses extended PDB identifiers
(`pdb_000010bt`); the legacy 4-character code is the last four characters, which is
how the two are joined here.

`in_analysis_set` is the column that matters for interpretation: the curated corpus
holds 774 antibody complexes, but only those clearing the >=50-residue antigen floor
(512) reached contact extraction. Peptide-epitope complexes are kept in the manifest
and flagged, rather than deleted, so the exclusion is visible instead of invisible.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import Counter, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

FIELDS = [
    "pdb_id", "group", "title", "resolution", "method", "deposited",
    "side_a_role", "side_a_description", "side_a_length", "side_a_chains",
    "side_b_role", "side_b_description", "side_b_length", "side_b_chains",
    "organism", "cluster30", "uniprot", "cluster_size", "unit_weight",
    "antibody_type", "heavy_species", "light_species", "heavy_subclass",
    "light_subclass", "antigen_type",
    "n_contacts", "in_analysis_set", "excluded_reason",
]


def _contact_counts(arm: str) -> Counter:
    c = Counter()
    with (C.INT / f"{arm}_contacts.csv").open() as fh:
        for row in csv.DictReader(fh):
            if arm == "ppi" and row.get("direction") != "a_to_b":
                continue          # count each complex once, not once per direction
            c[row["pdb_id"]] += 1
    return c


def _sabdab_index() -> dict[str, dict]:
    """{legacy 4-char PDB -> one representative SAbDab row}.

    A PDB can hold several antibody instances; fields used here (type, species,
    subclass) are properties of the construct and agree across instances, so the
    first instance is representative. Antigen types are unioned across instances
    because a single entry can present protein + sugar + ion antigens at once.
    """
    idx: dict[str, dict] = {}
    ag_types: dict[str, set] = defaultdict(set)
    path = C.MET / "sabdab_summary.csv"
    if not path.exists():
        return {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            pdb = row["PDB"][-4:].upper()
            for t in (row.get("antigen_type") or "").split("|"):
                if t and t != "NA":
                    ag_types[pdb].add(t)
            idx.setdefault(pdb, row)
    for pdb, row in idx.items():
        row["_antigen_types"] = "|".join(sorted(ag_types.get(pdb, ())))
    return idx


def antibody_rows(counts: Counter, sabdab: dict) -> list[dict]:
    corpus = json.loads((C.MET / "antibody_corpus.json").read_text())
    units = C.units("antibody")
    out = []
    for pdb, meta in sorted(corpus.items()):
        u = units.get(pdb, {})
        sd = sabdab.get(pdb.upper(), {})
        ents = meta.get("antigen_entities") or [{}]
        # the antigen entity the units file actually scored, else the longest
        ent = max(ents, key=lambda e: e.get("length") or 0)
        length = u.get("antigen_length", ent.get("length"))

        in_set = pdb in counts
        reason = ""
        if not in_set:
            reason = ("antigen shorter than the 50-residue folded-protein floor"
                      if (length or 0) < 50 else "not in the extracted contact table")

        out.append({
            "pdb_id": pdb,
            "group": "antibody",
            "title": meta.get("title", ""),
            "resolution": meta.get("resolution", u.get("resolution", "")),
            "method": sd.get("method", ""),
            "deposited": meta.get("deposited", ""),
            "side_a_role": "antibody (paratope)",
            "side_a_description": sd.get("type", "antibody"),
            "side_a_length": "",
            "side_a_chains": "|".join(meta.get("antibody_chains", [])),
            "side_b_role": "antigen (epitope)",
            "side_b_description": u.get("antigen_description", ent.get("description", "")),
            "side_b_length": length,
            "side_b_chains": "|".join(meta.get("antigen_chains", [])),
            "organism": "|".join(u.get("antigen_organism", ent.get("organism", [])) or []),
            "cluster30": u.get("cluster30", ""),
            "uniprot": u.get("uniprot", ""),
            "cluster_size": u.get("cluster_size", ""),
            "unit_weight": u.get("weight", ""),
            "antibody_type": sd.get("type", ""),
            "heavy_species": sd.get("heavy_species", ""),
            "light_species": sd.get("light_species", ""),
            "heavy_subclass": sd.get("heavy_subclass", ""),
            "light_subclass": sd.get("light_subclass", ""),
            "antigen_type": sd.get("_antigen_types", ""),
            "n_contacts": counts.get(pdb, 0),
            "in_analysis_set": int(in_set),
            "excluded_reason": reason,
        })
    return out


def ppi_rows(counts: Counter) -> list[dict]:
    corpus = json.loads((C.MET / "ppi_corpus.json").read_text())
    units = C.units("ppi")
    out = []
    for pdb, meta in sorted(corpus.items()):
        u = units.get(pdb, {})
        a, b = meta.get("partner_a", {}), meta.get("partner_b", {})
        ua, ub = u.get("a_to_b", {}), u.get("b_to_a", {})
        in_set = pdb in counts
        out.append({
            "pdb_id": pdb,
            "group": "ppi",
            "title": meta.get("title", ""),
            "resolution": meta.get("resolution", u.get("resolution", "")),
            "method": "",
            "deposited": "",
            "side_a_role": "partner A",
            "side_a_description": a.get("description", ""),
            "side_a_length": a.get("length", ""),
            "side_a_chains": "|".join(a.get("chains", [])),
            "side_b_role": "partner B",
            "side_b_description": b.get("description", ""),
            "side_b_length": b.get("length", ""),
            "side_b_chains": "|".join(b.get("chains", [])),
            "organism": "",
            # both directions are independence units; both are recorded
            "cluster30": "|".join(filter(None, [ua.get("cluster30"), ub.get("cluster30")])),
            "uniprot": "|".join(filter(None, [ua.get("uniprot"), ub.get("uniprot")])),
            "cluster_size": "",
            "unit_weight": "",
            "antibody_type": "", "heavy_species": "", "light_species": "",
            "heavy_subclass": "", "light_subclass": "", "antigen_type": "",
            "n_contacts": counts.get(pdb, 0),
            "in_analysis_set": int(in_set),
            "excluded_reason": "" if in_set else "not in the extracted contact table",
        })
    return out


def composition(rows: list[dict]) -> dict:
    out: dict = {}
    for group in ("antibody", "ppi"):
        g = [r for r in rows if r["group"] == group]
        used = [r for r in g if r["in_analysis_set"]]
        res = [float(r["resolution"]) for r in used
               if str(r["resolution"]) not in ("", "None")]
        res.sort()
        out[group] = {
            "n_in_manifest": len(g),
            "n_in_analysis_set": len(used),
            "n_excluded": len(g) - len(used),
            "excluded_reasons": Counter(r["excluded_reason"] for r in g
                                        if not r["in_analysis_set"]).most_common(),
            "resolution_median": round(res[len(res) // 2], 2) if res else None,
            "resolution_range": [res[0], res[-1]] if res else None,
            "top_partners": Counter(
                r["side_b_description"] for r in used
                if r["side_b_description"]).most_common(12),
        }
        if group == "antibody":
            out[group]["antibody_type"] = Counter(
                r["antibody_type"] or "unannotated" for r in used).most_common()
            out[group]["method"] = Counter(
                r["method"] or "unannotated" for r in used).most_common()
            out[group]["top_organisms"] = Counter(
                r["organism"] for r in used if r["organism"]).most_common(10)
            out[group]["sabdab_annotated"] = sum(1 for r in used if r["antibody_type"])
    return out


if __name__ == "__main__":
    sabdab = _sabdab_index()
    rows = (antibody_rows(_contact_counts("antibody"), sabdab)
            + ppi_rows(_contact_counts("ppi")))

    path = C.MET / "corpus_manifest.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    comp = composition(rows)
    comp["sabdab_index_size"] = len(sabdab)
    (C.OUT_MET / "corpus.json").write_text(json.dumps(comp, indent=1))

    print(f"wrote {path.relative_to(C.ROOT)}  ({len(rows)} structures)")
    for g in ("antibody", "ppi"):
        c = comp[g]
        print(f"\n{g}: {c['n_in_manifest']} in manifest, "
              f"{c['n_in_analysis_set']} in analysis set, "
              f"{c['n_excluded']} excluded  "
              f"(resolution median {c['resolution_median']} A)")
        for reason, n in c["excluded_reasons"]:
            print(f"    excluded {n:4d}: {reason}")
        print("    top partners:", ", ".join(f"{d} x{n}" for d, n in c["top_partners"][:6]))
    print("\nantibody types:", comp["antibody"]["antibody_type"])
    print("antibody methods:", comp["antibody"]["method"])
    print("SAbDab-annotated:", comp["antibody"]["sabdab_annotated"],
          "of", comp["antibody"]["n_in_analysis_set"])
