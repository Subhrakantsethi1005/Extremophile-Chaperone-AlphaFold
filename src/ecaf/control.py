"""Procedural control (proposal 3.7): compare a predicted model with an experimental structure."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.SeqUtils import seq1


def _ca_chain(path: Path, chain_id: str | None):
    parser = MMCIFParser(QUIET=True) if str(path).endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    model = parser.get_structure("x", str(path))[0]
    chain = model[chain_id] if chain_id else next(model.get_chains())
    res = [r for r in chain if r.id[0] == " " and "CA" in r]
    return res, "".join(seq1(r.get_resname()) for r in res)


def compare(model_pdb: Path, exp_path: Path, exp_chain: str | None = None, plddt_cut: float = 70.0) -> dict:
    m_res, m_seq = _ca_chain(model_pdb, None)
    e_res, e_seq = _ca_chain(exp_path, exp_chain)
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aln = al.align(m_seq, e_seq)[0]
    pairs = [(i, j) for (s1, e1), (s2, e2) in zip(*aln.aligned)
             for i, j in zip(range(s1, e1), range(s2, e2))]
    ident = sum(m_seq[i] == e_seq[j] for i, j in pairs) / max(len(pairs), 1)

    def rmsd(sel):
        if len(sel) < 3:
            return None
        sup = Superimposer()
        sup.set_atoms([e_res[j]["CA"] for _, j in sel], [m_res[i]["CA"] for i, _ in sel])
        return round(float(sup.rms), 3)

    confident = [(i, j) for i, j in pairs if m_res[i]["CA"].get_bfactor() >= plddt_cut]
    out = {"aligned_residues": len(pairs), "sequence_identity": round(ident, 3),
           "ca_rmsd_all": rmsd(pairs), "confident_residues": len(confident),
           "ca_rmsd_confident": rmsd(confident)}
    try:  # TM-score is optional (pip install tmtools)
        import tmtools
        a = np.array([m_res[i]["CA"].coord for i, _ in pairs], dtype=float)
        b = np.array([e_res[j]["CA"].coord for _, j in pairs], dtype=float)
        s = "".join(m_seq[i] for i, _ in pairs)
        r = tmtools.tm_align(a, b, s, "".join(e_seq[j] for _, j in pairs))
        out["tm_score"] = round(float(r.tm_norm_chain2), 3)
    except ImportError:
        pass
    return out
