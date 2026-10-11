"""Obtaining AlphaFold models (proposal 3.3) and normalising them to one layout.

Every route writes ``structures/<ACCESSION>/`` containing:
  model.pdb   top-ranked model, PDB format (pLDDT in the B-factor column)
  pae.json    {"pae": [[...]]} when the route provides PAE, otherwise absent
  meta.json   route, source files, software/model version, ranking metrics
"""

from __future__ import annotations

import json
import pickle
import zipfile
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.SeqUtils import seq1

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Return {accession: (header, sequence)}; accession taken from tr|ACC|NAME or first token."""
    out, header, seq = {}, None, []
    for line in Path(path).read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if header is not None:
                tok = header.split()[0]
                acc = tok.split("|")[1] if tok.count("|") >= 2 else tok
                out[acc] = (header, "".join(seq))
            header, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    return out


def write_server_jobs(fasta: Path, out_json: Path, seeds: list[int] | None = None) -> int:
    """AlphaFold Server batch-upload JSON (one monomer job per accession, default settings)."""
    jobs = []
    for acc, (_, seq) in read_fasta(fasta).items():
        jobs.append({
            "name": acc,
            "modelSeeds": [str(s) for s in (seeds or [])],
            "sequences": [{"proteinChain": {"sequence": seq, "count": 1}}],
            "dialect": "alphafoldserver",
            "version": 1,
        })
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(jobs, indent=2))
    return len(jobs)


def model_sequence(pdb_path: Path) -> str:
    s = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
    chain = next(s[0].get_chains())
    return "".join(seq1(r.get_resname()) for r in chain if r.id[0] == " ")


def _write_pdb(structure, path: Path) -> None:
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path))


def _finish(outdir: Path, meta: dict, expected_seq: str | None) -> dict:
    seq = model_sequence(outdir / "model.pdb")
    meta["model_length"] = len(seq)
    if expected_seq is not None:
        meta["sequence_matches_input"] = seq == expected_seq
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2))
    return meta


def normalise_pae(obj) -> list[list[float]] | None:
    """Accept the AFDB (v1-v4), AF2 pickle and AF3/Server PAE layouts."""
    if isinstance(obj, list) and obj and isinstance(obj[0], dict):
        obj = obj[0]
    if isinstance(obj, dict):
        if "predicted_aligned_error" in obj:
            return [list(map(float, r)) for r in obj["predicted_aligned_error"]]
        if "pae" in obj:
            return [list(map(float, r)) for r in obj["pae"]]
        if "distance" in obj and "residue1" in obj:  # AFDB v1/v2 flat layout
            n = max(obj["residue1"])
            m = [[0.0] * n for _ in range(n)]
            for i, j, d in zip(obj["residue1"], obj["residue2"], obj["distance"]):
                m[i - 1][j - 1] = float(d)
            return m
    if hasattr(obj, "tolist"):
        return obj.tolist()
    return None


def fetch_afdb(acc: str, outdir: Path, expected_seq: str | None = None, session=None) -> dict:
    """Download the precomputed AlphaFold DB model + PAE for an accession."""
    import requests

    s = session or requests.Session()
    r = s.get(AFDB_API.format(acc=acc), timeout=60)
    r.raise_for_status()
    entry = r.json()[0]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "model.pdb").write_bytes(s.get(entry["pdbUrl"], timeout=120).content)
    meta = {"route": "AlphaFold DB", "entry_id": entry.get("entryId"),
            "model_version": entry.get("latestVersion"), "model_created": entry.get("modelCreatedDate"),
            "afdb_global_plddt": entry.get("globalMetricValue"), "source_url": entry["pdbUrl"]}
    if entry.get("paeDocUrl"):
        pae = normalise_pae(s.get(entry["paeDocUrl"], timeout=120).json())
        (outdir / "pae.json").write_text(json.dumps({"pae": pae}))
    if expected_seq and entry.get("uniprotSequence") and entry["uniprotSequence"] != expected_seq:
        meta["warning"] = "AFDB sequence differs from curated FASTA; re-predict with AlphaFold Server/local AF2"
    return _finish(outdir, meta, expected_seq)


def import_server(src: Path, acc: str, outdir: Path, expected_seq: str | None = None) -> dict:
    """Import an AlphaFold Server download (zip or extracted folder); keeps model_0 (top ranked)."""
    src = Path(src)
    files = {}
    if src.suffix == ".zip":
        with zipfile.ZipFile(src) as z:
            for n in z.namelist():
                files[Path(n).name] = z.read(n)
    else:
        for p in src.rglob("*"):
            if p.is_file():
                files[p.name] = p.read_bytes()

    def pick(suffix):
        hits = sorted(n for n in files if n.endswith(suffix))
        return hits[0] if hits else None

    cif = pick("_model_0.cif")
    if not cif:
        raise FileNotFoundError(f"no *_model_0.cif in {src}")
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = outdir / "model_0.cif"
    tmp.write_bytes(files[cif])
    _write_pdb(MMCIFParser(QUIET=True).get_structure(acc, str(tmp)), outdir / "model.pdb")
    meta = {"route": "AlphaFold Server", "model_file": cif}
    full = pick("_full_data_0.json")
    if full:
        pae = normalise_pae(json.loads(files[full]))
        if pae:
            (outdir / "pae.json").write_text(json.dumps({"pae": pae}))
    summ = pick("_summary_confidences_0.json")
    if summ:
        sc = json.loads(files[summ])
        meta.update({k: sc.get(k) for k in ("ptm", "ranking_score", "fraction_disordered", "has_clash")})
    return _finish(outdir, meta, expected_seq)


def import_af2(src: Path, acc: str, outdir: Path, expected_seq: str | None = None) -> dict:
    """Import a local AlphaFold v2.3.2 output directory (ranked_0.pdb + ranking_debug.json)."""
    src = Path(src)
    outdir.mkdir(parents=True, exist_ok=True)
    structure = PDBParser(QUIET=True).get_structure(acc, str(src / "ranked_0.pdb"))
    _write_pdb(structure, outdir / "model.pdb")
    meta = {"route": "AlphaFold v2.3.2 (local)", "model_file": "ranked_0.pdb"}
    rd = src / "ranking_debug.json"
    if rd.exists():
        rank = json.loads(rd.read_text())
        top = rank["order"][0]
        meta["top_model"] = top
        metric = "plddts" if "plddts" in rank else "iptm+ptm"
        meta["ranking_metric"] = metric
        meta["ranking_value"] = rank[metric][top]
        pkl = src / f"result_{top}.pkl"
        if pkl.exists():
            try:
                with open(pkl, "rb") as fh:
                    res = pickle.load(fh)
                if "predicted_aligned_error" in res:
                    (outdir / "pae.json").write_text(json.dumps({"pae": normalise_pae(res["predicted_aligned_error"])}))
                if "ptm" in res:
                    meta["ptm"] = float(res["ptm"])
            except Exception as e:  # pickles can need jax/haiku to unpickle
                meta["pae_note"] = f"could not read {pkl.name}: {e}"
    return _finish(outdir, meta, expected_seq)
