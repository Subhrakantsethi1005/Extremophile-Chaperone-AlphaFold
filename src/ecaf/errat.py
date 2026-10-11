"""ERRAT non-bonded interaction quality (Colovos & Yeates, 1993; proposal 3.5).

ERRAT is distributed by UCLA-DOE as a standalone program and via the SAVES
server; its statistical reference matrices are not redistributed here, so this
module does not re-implement it. It (a) runs a local ERRAT binary when one is
configured, (b) parses ERRAT log output, and (c) reads a manual CSV for scores
taken from the SAVES report.
"""

from __future__ import annotations

import csv
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

OVERALL_RE = re.compile(r"Overall\s+quality\s+factor\**\s*[:=]?\s*([\d.]+)", re.I)


def parse_errat_log(text: str) -> float | None:
    m = OVERALL_RE.search(text)
    return float(m.group(1)) if m else None


def run_errat(pdb_path: Path, outdir: Path) -> float | None:
    """Run ERRAT if ECAF_ERRAT_CMD is set, e.g. ``ECAF_ERRAT_CMD="errat {stem}"``.

    The command runs in a scratch directory holding ``model.pdb``; ``{pdb}`` and
    ``{stem}`` are substituted. Any ``*.logf``/``*.log``/``*.ps`` produced is kept.
    """
    template = os.environ.get("ECAF_ERRAT_CMD")
    if not template:
        return None
    with tempfile.TemporaryDirectory() as td:
        shutil.copy(pdb_path, Path(td) / "model.pdb")
        cmd = template.format(pdb="model.pdb", stem="model")
        proc = subprocess.run(shlex.split(cmd), cwd=td, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, check=False)
        outdir.mkdir(parents=True, exist_ok=True)
        texts = [proc.stdout]
        for f in Path(td).iterdir():
            if f.suffix in (".logf", ".log", ".ps", ".pdf"):
                shutil.copy(f, outdir / f.name)
                if f.suffix in (".logf", ".log"):
                    texts.append(f.read_text(errors="ignore"))
        (outdir / "errat_stdout.txt").write_text(proc.stdout)
        for t in texts:
            v = parse_errat_log(t)
            if v is not None:
                return v
    return None


def load_errat(accession: str, validation_dir: Path) -> dict | None:
    """Look for validation/<ACC>/errat/*.logf|*.log|*.txt, then validation/errat_manual.csv."""
    d = validation_dir / accession / "errat"
    for f in sorted(list(d.glob("*.logf")) + list(d.glob("*.log")) + list(d.glob("*.txt"))):
        v = parse_errat_log(f.read_text(errors="ignore"))
        if v is not None:
            return {"errat_overall": v, "errat_source": f.name}
    manual = validation_dir / "errat_manual.csv"
    if manual.exists():
        with open(manual, newline="") as fh:
            for r in csv.DictReader(fh):
                if r["accession"] == accession:
                    return {"errat_overall": float(r["overall_quality_factor"]),
                            "errat_source": "manual (SAVES)",
                            "errat_flagged_regions": r.get("flagged_regions", "")}
    return None
