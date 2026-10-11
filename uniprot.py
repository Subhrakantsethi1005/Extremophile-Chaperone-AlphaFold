"""Candidate retrieval and curation from UniProtKB/TrEMBL (proposal 3.2)."""

from __future__ import annotations

import csv
import datetime as dt
import re
from dataclasses import dataclass, field, asdict
from pathlib import Path

from Bio import Align

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
TSV_FIELDS = [
    "accession", "id", "protein_name", "organism_name", "organism_id",
    "length", "annotation_score", "protein_existence", "keyword", "sequence",
]
AMBIGUOUS = set("XBZJUO")


@dataclass
class Candidate:
    accession: str
    entry_name: str
    protein_name: str
    organism: str
    taxon_id: str
    length: int
    annotation_score: str
    protein_existence: str
    keywords: str
    sequence: str
    organism_category: str = ""
    chaperone_class: str = ""
    exclusion_reason: str = ""
    notes: list[str] = field(default_factory=list)


def build_query(cfg: dict) -> str:
    chap = " OR ".join(f"({t})" for t in cfg["chaperone_query_terms"])
    orgs = " OR ".join(f'(organism_name:"{o["name"]}")' for o in cfg["organisms"])
    parts = ["(reviewed:false)", "(fragment:false)", f"({chap})", f"({orgs})"]
    if cfg.get("exclude_with_experimental_structure", True):
        parts.append("(structure_3d:false)")
    q = " AND ".join(parts)
    if cfg.get("exclude_protein_level_evidence", True):
        q += " NOT (existence:1)"
    return q


def fetch_uniprot(cfg: dict, max_records: int = 5000, session=None):
    """Query the UniProt REST API. Returns (candidates, provenance dict)."""
    import requests

    s = session or requests.Session()
    query = build_query(cfg)
    params = {"query": query, "format": "tsv", "fields": ",".join(TSV_FIELDS), "size": 500}
    url, rows, headers = UNIPROT_SEARCH, [], {}
    while url and len(rows) < max_records:
        r = s.get(url, params=params, timeout=120)
        r.raise_for_status()
        headers = headers or dict(r.headers)
        lines = r.text.splitlines()
        if lines:
            reader = csv.reader(lines[1:], delimiter="\t")
            rows.extend(reader)
        url, params = _next_link(r.headers.get("Link", "")), None
    cands = [_row_to_candidate(row) for row in rows if len(row) >= len(TSV_FIELDS)]
    prov = {
        "uniprot_query": query,
        "uniprot_release": headers.get("X-UniProt-Release", "unknown"),
        "uniprot_release_date": headers.get("X-UniProt-Release-Date", "unknown"),
        "access_date_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "records_returned": len(cands),
    }
    return cands, prov


def _next_link(link_header: str):
    m = re.search(r'<([^>]+)>;\s*rel="next"', link_header or "")
    return m.group(1) if m else None


def _row_to_candidate(row) -> Candidate:
    d = dict(zip(TSV_FIELDS, row))
    return Candidate(
        accession=d["accession"], entry_name=d["id"], protein_name=d["protein_name"],
        organism=d["organism_name"], taxon_id=d["organism_id"], length=int(d["length"] or 0),
        annotation_score=d["annotation_score"], protein_existence=d["protein_existence"],
        keywords=d["keyword"], sequence=d["sequence"].strip(),
    )


def read_tsv(path: Path) -> list[Candidate]:
    with open(path, newline="") as fh:
        return [_row_to_candidate([r[f] for f in TSV_FIELDS]) for r in csv.DictReader(fh, delimiter="\t")]


def classify(c: Candidate, cfg: dict) -> None:
    name = c.protein_name.lower()
    for cls in cfg["chaperone_classes"]:
        if any(p in name for p in cls["patterns"]):
            c.chaperone_class = cls["name"]
            break
    for o in cfg["organisms"]:
        if o["name"].lower() in c.organism.lower():
            c.organism_category = o["category"]
            break


def apply_filters(cands: list[Candidate], cfg: dict) -> list[Candidate]:
    """Annotate each candidate with an exclusion reason (empty = retained)."""
    classes = {c["name"]: c for c in cfg["chaperone_classes"]}
    for c in cands:
        classify(c, cfg)
        seq = c.sequence.upper()
        amb = sum(ch in AMBIGUOUS for ch in seq) / max(len(seq), 1)
        if not c.chaperone_class:
            c.exclusion_reason = "chaperone class not assignable from name"
        elif not c.organism_category:
            c.exclusion_reason = "organism not in extremophile list"
        elif amb > cfg["max_ambiguous_fraction"]:
            c.exclusion_reason = f"ambiguous residues {amb:.1%}"
        elif cfg.get("require_initial_met", True) and not seq.startswith("M"):
            c.exclusion_reason = "no initial Met (possible fragment)"
        elif not classes[c.chaperone_class]["min_len"] <= len(seq) <= classes[c.chaperone_class]["max_len"]:
            lo, hi = classes[c.chaperone_class]["min_len"], classes[c.chaperone_class]["max_len"]
            c.exclusion_reason = f"length {len(seq)} outside {lo}-{hi} for class"
    return cands


def identity(a: str, b: str) -> float:
    """Global-alignment identity over the shorter sequence length."""
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aln = al.align(a, b)[0]
    matches = sum(
        a[i] == b[j]
        for (s1, e1), (s2, e2) in zip(*aln.aligned)
        for i, j in zip(range(s1, e1), range(s2, e2))
    )
    return matches / min(len(a), len(b))


def select(cands: list[Candidate], cfg: dict) -> list[Candidate]:
    """Greedy, deterministic, diversity-first selection of non-redundant candidates.

    Round-robin over chaperone classes and organism categories so the final set
    spans as many classes/environments as possible; within a bucket the entry
    with the lowest annotation score (least characterized) and then the
    accession order wins. A candidate is rejected if it exceeds the identity
    threshold to any already-selected sequence.
    """
    pool = sorted(
        (c for c in cands if not c.exclusion_reason),
        key=lambda c: (float(c.annotation_score or 0), c.accession),
    )
    chosen: list[Candidate] = []
    used_classes, used_orgs = set(), set()
    while len(chosen) < cfg["target_n"] and pool:
        def priority(c):
            return (c.chaperone_class in used_classes, c.organism in used_orgs,
                    c.organism_category in {x.organism_category for x in chosen})
        pool.sort(key=priority)  # stable sort keeps the annotation-score order
        c = pool.pop(0)
        hit = next((x for x in chosen if identity(c.sequence, x.sequence) > cfg["max_pairwise_identity"]), None)
        if hit:
            c.exclusion_reason = f"redundant with {hit.accession}"
            continue
        chosen.append(c)
        used_classes.add(c.chaperone_class)
        used_orgs.add(c.organism)
    for c in pool:
        c.exclusion_reason = c.exclusion_reason or "not needed (target_n reached)"
    return chosen


def write_outputs(all_cands: list[Candidate], chosen: list[Candidate], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cols = [k for k in asdict(all_cands[0]).keys() if k not in ("sequence", "notes")] if all_cands else []
    with open(outdir / "screening_log.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols + ["selected"])
        sel = {c.accession for c in chosen}
        for c in all_cands:
            d = asdict(c)
            w.writerow([d[k] for k in cols] + [c.accession in sel])
    with open(outdir / "selected_candidates.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "entry_name", "protein_name", "organism", "organism_category",
                    "chaperone_class", "length", "annotation_score", "protein_existence"])
        for c in chosen:
            w.writerow([c.accession, c.entry_name, c.protein_name, c.organism, c.organism_category,
                        c.chaperone_class, c.length, c.annotation_score, c.protein_existence])
    with open(outdir / "selected_candidates.fasta", "w") as fh:
        for c in chosen:
            fh.write(f">tr|{c.accession}|{c.entry_name} {c.protein_name} OS={c.organism} OX={c.taxon_id}\n")
            fh.write("\n".join(c.sequence[i:i + 60] for i in range(0, len(c.sequence), 60)) + "\n")
