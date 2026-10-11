"""Primary confidence analysis from AlphaFold outputs (proposal 3.4)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

# AlphaFold DB / EMBL-EBI confidence bands
BANDS = [("very_high", 90, 101), ("confident", 70, 90), ("low", 50, 70), ("very_low", 0, 50)]


def residue_plddt(pdb_path: Path) -> tuple[list[int], list[str], np.ndarray]:
    """Per-residue pLDDT read from the CA B-factor (AF2, AF3/Server and AFDB all store it there)."""
    s = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
    chain = next(s[0].get_chains())
    nums, names, vals = [], [], []
    for r in chain:
        if r.id[0] != " " or "CA" not in r:
            continue
        nums.append(r.id[1])
        names.append(r.get_resname())
        vals.append(r["CA"].get_bfactor())
    v = np.array(vals, dtype=float)
    if v.size and v.max() <= 1.0:  # 0-1 scaled confidence (e.g. ESMFold)
        v = v * 100
    return nums, names, v


def segments(mask: np.ndarray, min_len: int) -> list[tuple[int, int]]:
    """0-based inclusive (start, end) runs where mask is True and run length >= min_len."""
    out, start = [], None
    for i, m in enumerate(list(mask) + [False]):
        if m and start is None:
            start = i
        elif not m and start is not None:
            if i - start >= min_len:
                out.append((start, i - 1))
            start = None
    return out


def analyse(pdb_path: Path, pae_path: Path | None = None, low_cut: float = 70.0,
            min_low_run: int = 5, min_domain: int = 30) -> dict:
    nums, names, p = residue_plddt(pdb_path)
    n = len(p)
    res = {
        "n_residues": n,
        "mean_plddt": round(float(p.mean()), 2),
        "median_plddt": round(float(np.median(p)), 2),
    }
    for name, lo, hi in BANDS:
        res[f"pct_{name}"] = round(100 * float(((p >= lo) & (p < hi)).mean()), 1)
    res["n_below_70"] = int((p < low_cut).sum())
    low = []
    for a, b in segments(p < low_cut, min_low_run):
        where = "N-terminus" if a == 0 else "C-terminus" if b == n - 1 else "internal"
        low.append({"start": nums[a], "end": nums[b], "length": b - a + 1,
                    "mean_plddt": round(float(p[a:b + 1].mean()), 1), "location": where})
    res["low_confidence_segments"] = low

    if pae_path and Path(pae_path).exists():
        pae = np.array(json.loads(Path(pae_path).read_text())["pae"], dtype=float)
        res["pae_available"] = True
        res["mean_pae"] = round(float(pae.mean()), 2)
        doms = segments(p >= low_cut, min_domain)
        res["confident_regions"] = [{"start": nums[a], "end": nums[b]} for a, b in doms]
        pairs = []
        for i in range(len(doms)):
            for j in range(i + 1, len(doms)):
                (a1, b1), (a2, b2) = doms[i], doms[j]
                inter = 0.5 * (pae[a1:b1 + 1, a2:b2 + 1].mean() + pae[a2:b2 + 1, a1:b1 + 1].mean())
                pairs.append({"region_a": f"{nums[a1]}-{nums[b1]}", "region_b": f"{nums[a2]}-{nums[b2]}",
                              "mean_inter_pae": round(float(inter), 2)})
        res["inter_region_pae"] = pairs
        res["uncertain_relative_placement"] = any(x["mean_inter_pae"] > 15 for x in pairs)
    else:
        res["pae_available"] = False
    res["_profile"] = {"resnum": nums, "resname": names, "plddt": [round(float(x), 2) for x in p]}
    return res


def plot(result: dict, title: str, out_png: Path, pae_path: Path | None = None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prof = result["_profile"]
    has_pae = pae_path is not None and Path(pae_path).exists()
    fig, axes = plt.subplots(1, 2 if has_pae else 1, figsize=(12 if has_pae else 8, 3.8),
                             gridspec_kw={"width_ratios": [1.6, 1]} if has_pae else None)
    ax = axes[0] if has_pae else axes
    for (name, lo, hi), col in zip(BANDS, ["#0053D6", "#65CBF3", "#FFDB13", "#FF7D45"]):
        ax.axhspan(lo, min(hi, 100), color=col, alpha=0.15, lw=0)
    ax.plot(prof["resnum"], prof["plddt"], color="#222", lw=1.2)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue")
    ax.set_ylabel("pLDDT")
    ax.set_title(f"{title}: mean pLDDT {result['mean_plddt']}")
    if has_pae:
        pae = np.array(json.loads(Path(pae_path).read_text())["pae"])
        im = axes[1].imshow(pae, cmap="Greens_r", vmin=0, vmax=30)
        axes[1].set_title("Predicted aligned error (Å)")
        axes[1].set_xlabel("Scored residue")
        axes[1].set_ylabel("Aligned residue")
        fig.colorbar(im, ax=axes[1], fraction=0.046)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
