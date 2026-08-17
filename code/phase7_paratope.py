"""Phase 7 - epitope/paratope logic: how the antibody actually divides the labour.

CDR assignment without ANARCI. SAbDab serves IMGT-numbered coordinates for the same
physical structures we already have, so instead of aligning sequences (ambiguous at
exactly the insertion-heavy loops we care about) we map local residues onto IMGT
residues by EXACT CA coordinate match. Verified on 5UMN: 1414/1414 CA atoms matched,
residue names agree, and the renumbering is real (local 64 -> IMGT 72).

IMGT regions (Lefranc): CDR1 27-38, CDR2 56-65, CDR3 105-117; everything else is
framework.

Questions answered:
  * how much of the paratope is CDR vs framework, per loop
  * is CDR-H3 dominant, and by how much
  * paratope residue composition against TWO nulls -- the general-PPI interface and
    the antibody's own surface -- because "paratopes are Tyr-rich" is only
    informative relative to a stated background
  * epitope character: how many sequence segments (discontinuous vs contiguous),
    burial, secondary structure
  * all of it computed WITHIN construct type, because 203 of 512 structures are
    nanobodies with no light chain, and pooling them would drag every L1/L2/L3
    statistic toward zero and make it look like a finding

Outputs
  metadata/analysis/paratope.json
  interactions/analysis/paratope_regions.csv
  figures/phase7_paratope.png
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
import time
import http.client
import urllib.error
import urllib.request
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import aa_common as C

SAB_DL = "https://sabdab.opig.stats.ox.ac.uk/api/download/structures"
IMGT_DIR = C.ROOT / "sabdab_imgt"
IMGT_DIR.mkdir(exist_ok=True)

_pdb = PDBParser(QUIET=True)
_cif = MMCIFParser(QUIET=True)

# IMGT unique numbering (Lefranc et al.)
REGIONS = [("FR1", 1, 26), ("CDR1", 27, 38), ("FR2", 39, 55), ("CDR2", 56, 65),
           ("FR3", 66, 104), ("CDR3", 105, 117), ("FR4", 118, 128)]


def region_of(num: int) -> str:
    for name, lo, hi in REGIONS:
        if lo <= num <= hi:
            return name
    return "other"


def download_imgt(pdb_ids):
    """Fetch SAbDab IMGT-numbered mmCIF for each entry; cached on disk."""
    todo = [p for p in pdb_ids if not (IMGT_DIR / f"{p}.cif").exists()]
    print(f"{len(pdb_ids) - len(todo)} cached, {len(todo)} to download", flush=True)
    ok = fail = 0
    for k, p in enumerate(todo, 1):
        body = json.dumps({"structures": [
            {"pdb_id": f"pdb_0000{p.lower()}", "full_structure": True}]}).encode()
        req = urllib.request.Request(SAB_DL, data=body,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=90) as fh:
                data = fh.read()
            if len(data) < 200:
                raise ValueError("response too small")
            (IMGT_DIR / f"{p}.cif").write_bytes(data)
            ok += 1
        except (urllib.error.URLError, TimeoutError, ValueError,
                http.client.HTTPException, OSError) as e:
            # a truncated response raises http.client.IncompleteRead, which is not a
            # URLError; catching only URLError let one bad download abort the run
            fail += 1
            if fail <= 5:
                print(f"  {p}: {type(e).__name__}", flush=True)
        if k % 50 == 0:
            print(f"  {k}/{len(todo)} (ok={ok} fail={fail})", flush=True)
        time.sleep(0.05)
    return ok, fail


def ca_index(model, rounded=2):
    out = {}
    for ch in model:
        for r in ch:
            if "CA" in r:
                a = r["CA"]
                if a.is_disordered():
                    a = a.disordered_get()
                out[tuple(np.round(a.coord, rounded))] = (
                    ch.id, f"{r.get_id()[1]}{r.get_id()[2].strip()}", r.get_id()[1])
    return out


def imgt_map(pdb_id):
    """{(local_chain, local_resi) -> (imgt_chain, imgt_number)} via CA coordinates."""
    cif = IMGT_DIR / f"{pdb_id}.cif"
    loc = None
    for nm in (f"{pdb_id}.pdb", f"{pdb_id}.cif"):
        p = C.ROOT / "pdb_antibody" / nm
        if p.exists():
            loc = p
            break
    if not cif.exists() or loc is None:
        return {}, 0.0
    try:
        S = ca_index(_cif.get_structure("s", str(cif))[0])
        L = ca_index((_cif if loc.suffix == ".cif" else _pdb)
                     .get_structure("l", str(loc))[0])
    except Exception:                                            # noqa: BLE001
        return {}, 0.0
    out = {}
    for coord, (lch, lresi, _) in L.items():
        hit = S.get(coord)
        if hit:
            out[(lch, lresi)] = (hit[0], hit[2])
    return out, len(out) / max(1, len(L))


def chain_roles():
    """{pdb -> {sabdab_chain -> 'H'|'L'}} from the SAbDab summary."""
    roles = defaultdict(dict)
    path = C.MET / "sabdab_summary.csv"
    if not path.exists():
        return roles
    for r in csv.DictReader(path.open()):
        p = r["PDB"][-4:].upper()
        if r.get("Hchain") and r["Hchain"] != "NA":
            roles[p][r["Hchain"]] = "H"
        if r.get("Lchain") and r["Lchain"] != "NA":
            roles[p][r["Lchain"]] = "L"
    return roles


def construct_types():
    out = {}
    for r in csv.DictReader((C.MET / "corpus_manifest.csv").open()):
        if r["group"] == "antibody" and r["in_analysis_set"] == "1":
            out[r["pdb_id"]] = r["antibody_type"] or "unannotated"
    return out


if __name__ == "__main__":
    rows = C.load_contacts("antibody")
    pdbs = sorted({r["pdb_id"] for r in rows})
    download_imgt(pdbs)

    roles = chain_roles()
    ctype = construct_types()

    maps, cover = {}, {}
    for k, p in enumerate(pdbs, 1):
        maps[p], cover[p] = imgt_map(p)
        if k % 100 == 0:
            print(f"mapped {k}/{len(pdbs)}", flush=True)

    # ---- annotate every paratope contact with its CDR / framework region
    per_contact = []
    unmapped = 0
    for r in rows:
        m = maps.get(r["pdb_id"]) or {}
        hit = m.get((r["ab_chain"], r["ab_resi"].strip()))
        if hit is None:
            unmapped += 1
            continue
        sab_chain, num = hit
        role = roles.get(r["pdb_id"], {}).get(sab_chain, "?")
        reg = region_of(num)
        label = (f"{'H' if role == 'H' else 'L' if role == 'L' else '?'}"
                 f"{reg[-1]}" if reg.startswith("CDR") else reg)
        per_contact.append({**r, "imgt": num, "region": reg, "chain_role": role,
                            "loop": label if reg.startswith("CDR") else "framework",
                            "ctype": ctype.get(r["pdb_id"], "unannotated")})

    print(f"\n{len(per_contact):,} paratope contacts mapped to IMGT "
          f"({unmapped:,} unmapped)")

    # ---- CDR vs framework share, per construct type, unit-weighted
    share = {}
    for ct in ("FAB", "FV", "SD-H", "ALL"):
        sub = [c for c in per_contact if ct == "ALL" or c["ctype"] == ct]
        if len(sub) < 200:
            continue
        w = defaultdict(float)
        for c in sub:
            w[c["loop"]] += c["weight"]
        tot = sum(w.values()) or 1.0
        share[ct] = {k: round(v / tot, 4) for k, v in
                     sorted(w.items(), key=lambda kv: -kv[1])}
        share[ct]["_n_contacts"] = len(sub)
        share[ct]["_n_structures"] = len({c["pdb_id"] for c in sub})
    # H3 dominance among CDR contacts only
    for ct, d in share.items():
        cdr = {k: v for k, v in d.items()
               if k.startswith(("H", "L")) and not k.startswith("_")}
        s = sum(cdr.values()) or 1.0
        d["_cdr_only_fractions"] = {k: round(v / s, 4) for k, v in cdr.items()}

    # ---- paratope composition vs two nulls
    def comp(rows_, key):
        w = defaultdict(float)
        for r in rows_:
            w[r[key]] += r["weight"]
        t = sum(w.values()) or 1.0
        return {a: w.get(a, 0.0) / t for a in C.AAS}

    para = comp(per_contact, "aa_probe")
    ppi_probe = comp(C.dedup_rows("ppi", C.load_contacts("ppi")), "aa_probe")
    # antibody's own surface as the second null: non-interface antibody residues
    surf = defaultdict(float)
    sres = list(csv.DictReader((C.OUT_INT / "structural_residues.csv").open()))
    for r in sres:
        if (r["arm"] == "antibody" and r["role"] == "probe"
                and r["is_interface"] == "0" and r["aa"] in C.AAS):
            try:
                if float(r["sasa_complex"] or 0) > 20:      # surface-exposed only
                    surf[r["aa"]] += 1.0
            except ValueError:
                pass
    st = sum(surf.values()) or 1.0
    ab_surface = {a: surf.get(a, 0.0) / st for a in C.AAS}

    enrich = {a: {"paratope": round(para[a], 4),
                  "vs_ppi_interface": round(float(np.log2(para[a] / ppi_probe[a]))
                                            if para[a] > 0 and ppi_probe[a] > 0
                                            else np.nan, 3),
                  "vs_antibody_surface": round(float(np.log2(para[a] / ab_surface[a]))
                                               if para[a] > 0 and ab_surface[a] > 0
                                               else np.nan, 3)}
              for a in C.AAS}

    # ---- epitope character: contiguity, burial, secondary structure
    ss_by_res = {(r["arm"], r["pdb_id"], r["chain"], r["resi"]): r for r in sres}
    epi = defaultdict(lambda: defaultdict(set))
    for r in rows:
        epi[r["pdb_id"]][r["ag_chain"]].add(r["ag_resi"].strip())
    seg_counts, sizes = [], []
    for p, chs in epi.items():
        n_seg = 0
        n_res = 0
        for ch, resis in chs.items():
            nums = sorted({int("".join(filter(str.isdigit, x)) or 0) for x in resis})
            n_res += len(nums)
            if not nums:
                continue
            n_seg += 1 + sum(1 for a, b in zip(nums, nums[1:]) if b - a > 4)
        seg_counts.append(n_seg)
        sizes.append(n_res)

    out = {
        "n_contacts_mapped": len(per_contact),
        "n_contacts_unmapped": unmapped,
        "imgt_mapping_coverage_median": round(float(np.median(list(cover.values()))), 4),
        "cdr_vs_framework": share,
        "paratope_composition": enrich,
        "epitope": {
            "median_segments": float(np.median(seg_counts)),
            "mean_segments": round(float(np.mean(seg_counts)), 2),
            "frac_discontinuous_ge2_segments": round(
                float(np.mean(np.array(seg_counts) >= 2)), 4),
            "median_epitope_residues": float(np.median(sizes)),
        },
        "note": ("CDR shares are unit-weighted and computed WITHIN construct type; "
                 "SD-H (nanobody) has no light chain, so pooling it with FAB would "
                 "depress every L1/L2/L3 statistic artefactually. Paratope "
                 "composition is reported against two nulls because enrichment is "
                 "meaningless without a stated background."),
    }
    (C.OUT_MET / "paratope.json").write_text(json.dumps(out, indent=1))

    with (C.OUT_INT / "paratope_regions.csv").open("w") as fh:
        fh.write("pdb_id,ctype,ab_chain,ab_resi,imgt,chain_role,region,loop,aa,weight\n")
        for c in per_contact:
            fh.write(f"{c['pdb_id']},{c['ctype']},{c['ab_chain']},{c['ab_resi']},"
                     f"{c['imgt']},{c['chain_role']},{c['region']},{c['loop']},"
                     f"{c['aa_probe']},{c['weight']:.6f}\n")

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.2))

    ax = axes[0]
    loops = ["H1", "H2", "H3", "L1", "L2", "L3", "framework"]
    cts = [c for c in ("FAB", "SD-H") if c in share]
    width = 0.8 / max(1, len(cts))
    for i, ct in enumerate(cts):
        v = [share[ct].get(l, 0.0) for l in loops]
        ax.bar(np.arange(len(loops)) + i * width - 0.4 + width / 2, v, width=width,
               label=f"{ct} (n={share[ct]['_n_structures']})")
    ax.set_xticks(range(len(loops))); ax.set_xticklabels(loops, fontsize=9)
    ax.set_ylabel("share of paratope contacts (unit-weighted)")
    ax.set_title("Where the paratope contacts come from", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    order = sorted(C.AAS, key=lambda a: -enrich[a]["vs_ppi_interface"]
                   if np.isfinite(enrich[a]["vs_ppi_interface"]) else 0)
    x = np.arange(len(order))
    ax.bar(x - 0.2, [enrich[a]["vs_ppi_interface"] for a in order], width=0.4,
           color="#b2182b", label="vs PPI interface")
    ax.bar(x + 0.2, [enrich[a]["vs_antibody_surface"] for a in order], width=0.4,
           color="#2166ac", label="vs antibody's own surface")
    ax.axhline(0, c="0.4", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(order, fontsize=8)
    for t, a in zip(ax.get_xticklabels(), order):
        t.set_color(C.GROUP_COLOR[C.AA_GROUP[a]])
    ax.set_ylabel("log$_2$ enrichment in paratope")
    ax.set_title("Paratope composition against two nulls", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[2]
    ax.hist(seg_counts, bins=range(1, 12), align="left", color="#1b7837", alpha=0.8)
    ax.set_xlabel("sequence segments in the epitope")
    ax.set_ylabel("structures")
    ax.set_title(f"Epitope contiguity\n"
                 f"{out['epitope']['frac_discontinuous_ge2_segments']:.0%} "
                 f"discontinuous (>=2 segments)", fontsize=10)

    fig.suptitle("Phase 7 - paratope architecture and epitope character", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(C.OUT_FIG / "phase7_paratope.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\nIMGT mapping coverage (median): "
          f"{out['imgt_mapping_coverage_median']:.3f}")
    print("\nCDR vs framework share of paratope contacts:")
    for ct, d in share.items():
        core = {k: v for k, v in d.items() if not k.startswith("_")}
        print(f"  {ct:10s} (n={d['_n_structures']:3d} structures) " +
              "  ".join(f"{k}={v:.3f}" for k, v in core.items()))
    print("\nCDR-only fractions (H3 dominance):")
    for ct, d in share.items():
        print(f"  {ct:10s} {d['_cdr_only_fractions']}")
    print("\ntop paratope enrichment vs PPI interface:")
    for a in order[:8]:
        print(f"   {a}  vs_ppi={enrich[a]['vs_ppi_interface']:+.2f}  "
              f"vs_ab_surface={enrich[a]['vs_antibody_surface']:+.2f}")
    print(f"\nepitope: median {out['epitope']['median_segments']:.0f} segments, "
          f"{out['epitope']['frac_discontinuous_ge2_segments']:.1%} discontinuous, "
          f"median {out['epitope']['median_epitope_residues']:.0f} residues")
