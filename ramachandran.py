"""Backbone conformational validation (proposal 3.5).

Two independent sources of Ramachandran statistics are supported:

1. PROCHECK v3.5 (the method named in the proposal). PROCHECK is licensed
   software; run it locally (``procheck model.pdb 2.0``) or through the SAVES
   server, then this module parses the ``*.sum`` file or a manual CSV.
2. A built-in screen that classifies every residue against the Lovell et al.
   (2003) Top500 phi/psi distributions (as packaged in PyRAMA): favoured
   (>= 2% contour), allowed, or outlier. It always runs, maps every unusual
   residue back to the sequence, and lets outliers be cross-referenced with
   pLDDT. It is a MolProbity-style cross-check, *not* a substitute for
   PROCHECK's core/allowed/generous/disallowed regions.
"""

from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PPBuilder

TABLES = {"General": "pref_general.data", "Gly": "pref_glycine.data",
          "Pro": "pref_proline.data", "PrePro": "pref_preproline.data"}
FAVOURED = 0.02
ALLOWED = {"General": 0.0005, "Gly": 0.002, "Pro": 0.002, "PrePro": 0.002}


@lru_cache(maxsize=None)
def _table(kind: str) -> np.ndarray:
    from importlib.resources import files

    grid = np.zeros((180, 180))
    path = files("pyrama") / "data" / TABLES[kind]
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        phi, psi, val = map(float, line.split()[:3])
        grid[_bin(phi), _bin(psi)] = val
    return grid


def _bin(angle: float) -> int:
    return int(math.floor((angle + 180.0) / 2.0)) % 180


def residue_class(resname: str, next_resname: str | None) -> str:
    if resname == "GLY":
        return "Gly"
    if resname == "PRO":
        return "Pro"
    if next_resname == "PRO":
        return "PrePro"
    return "General"


def phi_psi_table(pdb_path: Path) -> list[dict]:
    s = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
    rows = []
    for chain in s[0]:
        for pp in PPBuilder().build_peptides(chain):
            angles = pp.get_phi_psi_list()
            for i, (res, (phi, psi)) in enumerate(zip(pp, angles)):
                if phi is None or psi is None:
                    continue
                nxt = pp[i + 1].get_resname() if i + 1 < len(pp) else None
                kind = residue_class(res.get_resname(), nxt)
                phi_d, psi_d = math.degrees(phi), math.degrees(psi)
                dens = float(_table(kind)[_bin(phi_d), _bin(psi_d)])
                region = "favoured" if dens >= FAVOURED else "allowed" if dens >= ALLOWED[kind] else "outlier"
                rows.append({"chain": chain.id, "resnum": res.id[1], "resname": res.get_resname(),
                             "class": kind, "phi": round(phi_d, 1), "psi": round(psi_d, 1),
                             "density": dens, "region": region,
                             "plddt": round(res["CA"].get_bfactor(), 2) if "CA" in res else None})
    return rows


def summarise(rows: list[dict], min_plddt: float | None = None) -> dict:
    sel = [r for r in rows if min_plddt is None or (r["plddt"] is not None and r["plddt"] >= min_plddt)]
    n = len(sel)
    count = {k: sum(r["region"] == k for r in sel) for k in ("favoured", "allowed", "outlier")}
    out = {"n_residues": n}
    for k, v in count.items():
        out[f"pct_{k}"] = round(100 * v / n, 2) if n else None
    out["outliers"] = [f"{r['resname']}{r['resnum']}" for r in sel if r["region"] == "outlier"]
    return out


def write_csv(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["resnum"])
        w.writeheader()
        w.writerows(rows)


def plot(rows: list[dict], title: str, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors

    fig, axes = plt.subplots(2, 2, figsize=(9, 8.5))
    for ax, kind in zip(axes.flat, ["General", "Gly", "Pro", "PrePro"]):
        g = _table(kind).T  # rows = psi, cols = phi
        lv = colors.BoundaryNorm([0, ALLOWED[kind], FAVOURED, 1.01], 3)
        ax.imshow(g, origin="lower", extent=(-180, 180, -180, 180),
                  cmap=colors.ListedColormap(["#FFFFFF", "#CFE8FC", "#7FB8E8"]), norm=lv)
        pts = [r for r in rows if r["class"] == kind]
        ok = [r for r in pts if r["region"] != "outlier"]
        bad = [r for r in pts if r["region"] == "outlier"]
        ax.scatter([r["phi"] for r in ok], [r["psi"] for r in ok], s=6, c="#222")
        ax.scatter([r["phi"] for r in bad], [r["psi"] for r in bad], s=14, c="#D62728")
        for r in bad:
            ax.annotate(f"{r['resname']}{r['resnum']}", (r["phi"], r["psi"]), fontsize=6, color="#D62728")
        ax.set_title(f"{kind} (n={len(pts)})")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xticks(range(-180, 181, 90))
        ax.set_yticks(range(-180, 181, 90))
        ax.axhline(0, c="#999", lw=0.5)
        ax.axvline(0, c="#999", lw=0.5)
        ax.set_xlabel("φ (°)")
        ax.set_ylabel("ψ (°)")
    fig.suptitle(f"{title}: Ramachandran (Top500 favoured/allowed contours)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- PROCHECK --
SUM_RE = re.compile(
    r"Ramachandran plot:\s*([\d.]+)%\s*core\s*([\d.]+)%\s*allow\s*([\d.]+)%\s*gener\s*([\d.]+)%\s*disall",
    re.I)


def parse_procheck_sum(text: str) -> dict | None:
    m = SUM_RE.search(text)
    if not m:
        return None
    core, allow, gener, disall = map(float, m.groups())
    return {"procheck_core": core, "procheck_allowed": allow,
            "procheck_generous": gener, "procheck_disallowed": disall}


def run_procheck(pdb_path: Path, outdir: Path, resolution: str = "2.0") -> dict | None:
    """Run a local PROCHECK v3.5 install if `procheck` is on PATH; returns parsed .sum or None."""
    exe = shutil.which("procheck")
    if not exe:
        return None
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "model.pdb"
        shutil.copy(pdb_path, tmp)
        subprocess.run([exe, tmp.name, resolution], cwd=td, check=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        outdir.mkdir(parents=True, exist_ok=True)
        for f in Path(td).iterdir():
            if f.suffix in (".sum", ".ps", ".out", ".log") or f.name.endswith(".rin"):
                shutil.copy(f, outdir / f.name)
        sums = list(Path(td).glob("*.sum"))
        return parse_procheck_sum(sums[0].read_text()) if sums else None


def load_procheck(accession: str, validation_dir: Path) -> dict | None:
    """Look for PROCHECK results: validation/<ACC>/procheck/*.sum, then validation/procheck_manual.csv."""
    for f in sorted((validation_dir / accession / "procheck").glob("*.sum")):
        parsed = parse_procheck_sum(f.read_text())
        if parsed:
            return parsed
    manual = validation_dir / "procheck_manual.csv"
    if manual.exists():
        with open(manual, newline="") as fh:
            for r in csv.DictReader(fh):
                if r["accession"] == accession:
                    return {f"procheck_{k}": float(r[k]) for k in ("core", "allowed", "generous", "disallowed")}
    return None
