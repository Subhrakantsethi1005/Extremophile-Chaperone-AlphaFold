"""Comparative model-quality table and decision criteria (proposal 3.6)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

# Thresholds are conventions, recorded with the results; edit here if the group agrees on others.
CRITERIA = {
    "mean_plddt_strong": 80.0,
    "mean_plddt_min": 70.0,
    "very_low_pct_max": 10.0,
    "rama_favoured_strong": 95.0,   # built-in Top500 screen, confident residues
    "procheck_core_strong": 90.0,   # Laskowski et al. 1993: good models > 90% in core regions
    "procheck_disallowed_max": 1.0,
    "errat_acceptable": 85.0,
    "errat_poor": 70.0,
}


def assess(r: dict) -> tuple[str, list[str]]:
    """Return (evidence tier, concerns). Metrics are combined by rule, never averaged."""
    c, k = [], CRITERIA
    conf, rama_c, rama_a = r["confidence"], r["rama_confident"], r["rama_all"]
    if conf["mean_plddt"] < k["mean_plddt_min"]:
        c.append(f"low overall confidence (mean pLDDT {conf['mean_plddt']})")
    if conf["pct_very_low"] > k["very_low_pct_max"]:
        c.append(f"{conf['pct_very_low']}% residues pLDDT<50 (disordered-like or unreliable)")
    for s in conf["low_confidence_segments"]:
        if s["location"] == "internal" and s["length"] >= 10:
            c.append(f"internal low-confidence segment {s['start']}-{s['end']}")
    if conf.get("uncertain_relative_placement"):
        c.append("PAE indicates uncertain relative placement of confident regions")
    unexplained = rama_c["outliers"]
    explained = [x for x in rama_a["outliers"] if x not in unexplained]
    if unexplained:
        c.append(f"Ramachandran outliers in confident regions: {', '.join(unexplained)}")
    if explained:
        c.append(f"outliers only in low-pLDDT regions (expected): {', '.join(explained)}")
    pc, er = r.get("procheck"), r.get("errat")
    if pc:
        if pc["procheck_core"] < k["procheck_core_strong"]:
            c.append(f"PROCHECK core {pc['procheck_core']}% < {k['procheck_core_strong']}%")
        if pc["procheck_disallowed"] > k["procheck_disallowed_max"]:
            c.append(f"PROCHECK disallowed {pc['procheck_disallowed']}%")
    if er and er["errat_overall"] < k["errat_acceptable"]:
        c.append(f"ERRAT {er['errat_overall']} < {k['errat_acceptable']}")
    if er and conf["mean_plddt"] >= k["mean_plddt_strong"] and er["errat_overall"] < k["errat_acceptable"]:
        c.append("metric disagreement: high pLDDT but weak ERRAT")
    if pc and rama_c["pct_favoured"] is not None and rama_c["pct_favoured"] >= k["rama_favoured_strong"] \
            and pc["procheck_core"] < k["procheck_core_strong"]:
        c.append("metric disagreement: Top500 screen favourable but PROCHECK core low")

    weak = (conf["mean_plddt"] < k["mean_plddt_min"]
            or (er and er["errat_overall"] < k["errat_poor"])
            or (pc and pc["procheck_disallowed"] > 2 * k["procheck_disallowed_max"]))
    strong = (conf["mean_plddt"] >= k["mean_plddt_strong"]
              and (rama_c["pct_favoured"] or 0) >= k["rama_favoured_strong"]
              and not unexplained
              and (not pc or (pc["procheck_core"] >= k["procheck_core_strong"]
                              and pc["procheck_disallowed"] <= k["procheck_disallowed_max"]))
              and (not er or er["errat_overall"] >= k["errat_acceptable"]))
    tier = "weaker" if weak else "stronger" if strong else "moderate"
    pending = [n for n, v in (("PROCHECK", pc), ("ERRAT", er)) if not v]
    if pending:
        tier += f" (provisional: {' & '.join(pending)} pending)"
    return tier, c


def table_rows(results: list[dict]) -> list[dict]:
    rows = []
    for r in results:
        tier, concerns = assess(r)
        info, conf, pc, er = r["info"], r["confidence"], r.get("procheck") or {}, r.get("errat") or {}
        rows.append({
            "accession": r["accession"],
            "organism": info.get("organism", ""),
            "class": info.get("chaperone_class", ""),
            "length": conf["n_residues"],
            "route": r["meta"].get("route", ""),
            "mean_pLDDT": conf["mean_plddt"],
            "pct_pLDDT_ge90": conf["pct_very_high"],
            "residues_pLDDT_lt70": conf["n_below_70"],
            "mean_PAE": conf.get("mean_pae", "n/a"),
            "rama_favoured_pct_all": r["rama_all"]["pct_favoured"],
            "rama_outlier_pct_all": r["rama_all"]["pct_outlier"],
            "rama_favoured_pct_pLDDT70": r["rama_confident"]["pct_favoured"],
            "procheck_core_pct": pc.get("procheck_core", "pending"),
            "procheck_disallowed_pct": pc.get("procheck_disallowed", "pending"),
            "errat_overall": er.get("errat_overall", "pending"),
            "evidence_tier": tier,
            "major_concerns": "; ".join(concerns) or "none",
        })
    return rows


def write_table(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def markdown_report(results: list[dict], rows: list[dict], provenance: dict, path: Path) -> None:
    L = ["# Structural prediction and conformational validation: results", ""]
    L += ["## Provenance", "", "```json", json.dumps(provenance, indent=2), "```", ""]
    L += ["## Comparative model-quality table", ""]
    cols = ["accession", "organism", "class", "length", "mean_pLDDT", "residues_pLDDT_lt70",
            "rama_favoured_pct_all", "rama_outlier_pct_all", "procheck_core_pct",
            "procheck_disallowed_pct", "errat_overall", "evidence_tier"]
    L.append("| " + " | ".join(cols) + " |")
    L.append("|" + "---|" * len(cols))
    for row in rows:
        L.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    L += ["", "Candidates are not ranked by a single metric. The evidence tier follows the rules in "
          "`ecaf/summary.py` (thresholds below); disagreements between metrics are listed, not averaged.", "",
          "```json", json.dumps(CRITERIA, indent=2), "```", ""]
    for r, row in zip(results, rows):
        acc, conf = r["accession"], r["confidence"]
        L += [f"## {acc}: {r['info'].get('protein_name', '')}", ""]
        L.append(f"*{r['info'].get('organism', '')}*, {conf['n_residues']} residues, model route: "
                 f"{r['meta'].get('route', '?')}.")
        L += ["", f"![pLDDT/PAE]({acc}/confidence.png)", "", f"![Ramachandran]({acc}/ramachandran.png)", ""]
        L.append(f"- pLDDT bands: ≥90 {conf['pct_very_high']}%, 70–90 {conf['pct_confident']}%, "
                 f"50–70 {conf['pct_low']}%, <50 {conf['pct_very_low']}%")
        segs = ", ".join(f"{s['start']}-{s['end']} ({s['location']}, mean {s['mean_plddt']})"
                         for s in conf["low_confidence_segments"]) or "none"
        L.append(f"- Low-confidence segments (pLDDT<70, ≥5 residues): {segs}")
        L.append(f"- Ramachandran (Top500 screen) all residues: {r['rama_all']['pct_favoured']}% favoured, "
                 f"{r['rama_all']['pct_outlier']}% outliers; pLDDT≥70 only: "
                 f"{r['rama_confident']['pct_favoured']}% favoured")
        if r.get("control"):
            L.append(f"- Procedural control vs experimental structure: {json.dumps(r['control'])}")
        L.append(f"- Concerns: {row['major_concerns']}")
        L.append("")
    L += ["## Interpretation limits", "",
          "These are computational models. pLDDT is a local confidence estimate, not experimental accuracy; "
          "PROCHECK, ERRAT and the Ramachandran screen assess geometry and interaction patterns only and do not "
          "establish function or prove that a fold exists experimentally.", ""]
    path.write_text("\n".join(L))
