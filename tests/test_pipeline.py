"""Offline tests. 1A8O (HIV capsid C-terminal domain, 1.7 Å crystal structure) is a real
experimental structure used as a geometry reference; AlphaFold-style inputs are built from it
with synthetic pLDDT values so every code path runs without network access."""

import json
import zipfile
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBIO, PDBParser

from ecaf import confidence, control, errat, ramachandran, structures, summary, uniprot
from ecaf.__main__ import main

DATA = Path(__file__).parent / "data"
XTAL = DATA / "1A8O.pdb"
CFG = json.loads((Path(__file__).parents[1] / "config" / "selection.json").read_text())


def af_like(tmp_path, low_first=8):
    """Copy of 1A8O with pLDDT in B-factors: first `low_first` residues 40, rest 92."""
    s = PDBParser(QUIET=True).get_structure("x", str(XTAL))
    std = [r for r in next(s[0].get_chains()) if r.id[0] == " "]
    for i, r in enumerate(std):
        for at in r:
            at.set_bfactor(40.0 if i < low_first else 92.0)
    io = PDBIO()
    io.set_structure(s)
    out = tmp_path / "model.pdb"
    io.save(str(out))
    return out, s


def test_ramachandran_crystal_structure_is_mostly_favoured():
    rows = ramachandran.phi_psi_table(XTAL)
    s = ramachandran.summarise(rows)
    assert s["n_residues"] > 55
    assert s["pct_favoured"] >= 95
    assert s["pct_outlier"] <= 2


def test_ramachandran_known_regions():
    t = ramachandran._table("General")
    assert t[ramachandran._bin(-63), ramachandran._bin(-41)] > 0.9          # alpha helix
    assert t[ramachandran._bin(-119), ramachandran._bin(131)] > 0.5         # beta strand
    assert t[ramachandran._bin(60), ramachandran._bin(-150)] < 0.02         # left-handed, disfavoured
    assert ramachandran.residue_class("ALA", "PRO") == "PrePro"
    assert ramachandran.residue_class("GLY", "PRO") == "Gly"


def test_confidence_bands_and_segments(tmp_path):
    model, _ = af_like(tmp_path)
    n = len(confidence.residue_plddt(model)[0])
    pae = np.full((n, n), 3.0)
    pae_path = tmp_path / "pae.json"
    pae_path.write_text(json.dumps({"pae": pae.tolist()}))
    r = confidence.analyse(model, pae_path, min_domain=20)
    assert r["n_below_70"] == 8
    assert r["low_confidence_segments"][0]["location"] == "N-terminus"
    assert pytest.approx(r["pct_very_low"], abs=0.1) == round(100 * 8 / n, 1)
    assert r["pae_available"] and not r["uncertain_relative_placement"]
    confidence.plot(r, "t", tmp_path / "c.png", pae_path)
    assert (tmp_path / "c.png").stat().st_size > 1000


def test_pae_layouts():
    assert structures.normalise_pae([{"predicted_aligned_error": [[0, 1], [1, 0]]}]) == [[0, 1], [1, 0]]
    flat = {"residue1": [1, 1, 2, 2], "residue2": [1, 2, 1, 2], "distance": [0, 5, 6, 0]}
    assert structures.normalise_pae(flat) == [[0, 5], [6, 0]]
    assert structures.normalise_pae({"pae": [[1]]}) == [[1.0]]


def test_procheck_and_errat_parsers():
    line = "| Ramachandran plot:   91.7% core    7.6% allow    0.7% gener    0.0% disall  |"
    assert ramachandran.parse_procheck_sum(line) == {
        "procheck_core": 91.7, "procheck_allowed": 7.6, "procheck_generous": 0.7, "procheck_disallowed": 0.0}
    assert errat.parse_errat_log("...\n# Overall quality factor**:  87.931\n") == 87.931
    assert errat.parse_errat_log("nothing") is None


def _cand(acc, name, org, seq, score="1.0"):
    return uniprot.Candidate(acc, acc + "_X", name, org, "1", len(seq), score, "3: Inferred from homology",
                             "", seq)


def test_filters_and_redundancy():
    rng = np.random.default_rng(0)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    s1 = "M" + "".join(rng.choice(aa, 150))
    s2 = "M" + "".join(rng.choice(aa, 150))
    near_s1 = s1[:140] + "".join(rng.choice(aa, 11))
    cands = [
        _cand("A1", "Hsp20/alpha crystallin family protein", "Pyrococcus sp. 1", s1),
        _cand("A2", "Hsp20 family protein", "Thermococcus sp. 2", near_s1),          # redundant
        _cand("A3", "Prefoldin subunit beta", "Sulfolobus sp. 3", s2[:140]),
        _cand("A4", "Hsp20 family protein", "Escherichia coli", s2),                 # not extremophile
        _cand("A5", "Hsp20 family protein", "Pyrococcus sp. 5", s2[:60]),            # too short
        _cand("A6", "Hsp20 family protein", "Pyrococcus sp. 6", "A" + s2[1:]),       # no Met
        _cand("A7", "Hsp20 family protein", "Pyrococcus sp. 7", "M" + "X" * 20 + s2[21:]),
    ]
    uniprot.apply_filters(cands, CFG)
    reasons = {c.accession: c.exclusion_reason for c in cands}
    assert reasons["A4"].startswith("organism")
    assert reasons["A5"].startswith("length")
    assert reasons["A6"].startswith("no initial Met")
    assert reasons["A7"].startswith("ambiguous")
    chosen = uniprot.select(cands, dict(CFG, target_n=5))
    assert {c.accession for c in chosen} == {"A1", "A3"}
    assert reasons["A2"] == "" and cands[1].exclusion_reason == "redundant with A1"
    q = uniprot.build_query(CFG)
    assert "reviewed:false" in q and "KW-0143" in q and "Pyrococcus" in q


def test_server_import_and_jobs(tmp_path):
    _, s = af_like(tmp_path)
    io = MMCIFIO()
    io.set_structure(s)
    cif = tmp_path / "fold_x_model_0.cif"
    io.save(str(cif))
    n = len(list(next(s[0].get_chains())))
    zp = tmp_path / "fold.zip"
    with zipfile.ZipFile(zp, "w") as z:
        z.write(cif, "fold_x/fold_x_model_0.cif")
        z.writestr("fold_x/fold_x_full_data_0.json", json.dumps({"pae": np.ones((n, n)).tolist()}))
        z.writestr("fold_x/fold_x_summary_confidences_0.json", json.dumps({"ptm": 0.8, "ranking_score": 0.85}))
    meta = structures.import_server(zp, "TEST1", tmp_path / "out")
    assert meta["ptm"] == 0.8 and (tmp_path / "out" / "pae.json").exists()
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">tr|TEST1|T_X test\nMKV\n")
    out = tmp_path / "jobs.json"
    assert structures.write_server_jobs(fasta, out) == 1
    assert json.loads(out.read_text())[0]["sequences"][0]["proteinChain"]["sequence"] == "MKV"


def test_control_self_comparison(tmp_path):
    model, _ = af_like(tmp_path)
    r = control.compare(model, XTAL)
    assert r["sequence_identity"] == 1.0
    assert r["ca_rmsd_all"] < 1e-3
    if "tm_score" in r:
        assert r["tm_score"] > 0.99


def test_assess_rules():
    base = {"confidence": {"mean_plddt": 88, "pct_very_low": 2, "low_confidence_segments": []},
            "rama_all": {"outliers": ["GLY5"], "pct_favoured": 97}, "rama_confident": {"outliers": [], "pct_favoured": 98}}
    tier, c = summary.assess(dict(base, procheck={"procheck_core": 92, "procheck_disallowed": 0},
                                  errat={"errat_overall": 93}))
    assert tier == "stronger" and any("low-pLDDT" in x for x in c)
    tier, c = summary.assess(dict(base, procheck=None, errat={"errat_overall": 60}))
    assert tier.startswith("weaker") and any("disagreement" in x for x in c)


def test_cli_end_to_end(tmp_path, monkeypatch):
    monkeypatch.delenv("ECAF_ERRAT_CMD", raising=False)
    model, _ = af_like(tmp_path)
    sd = tmp_path / "structures" / "TEST1"
    sd.mkdir(parents=True)
    model.rename(sd / "model.pdb")
    (sd / "meta.json").write_text(json.dumps({"route": "test fixture"}))
    vd = tmp_path / "validation"
    vd.mkdir()
    (vd / "errat_manual.csv").write_text("accession,overall_quality_factor,flagged_regions\nTEST1,91.2,\n")
    (vd / "procheck_manual.csv").write_text("accession,core,allowed,generous,disallowed\nTEST1,93.1,6.9,0,0\n")
    res = tmp_path / "results"
    main(["--structures", str(tmp_path / "structures"), "--validation", str(vd), "--results", str(res),
          "--data", str(tmp_path / "data"), "analyze"])
    table = (res / "model_quality_table.csv").read_text()
    assert "TEST1" in table and "91.2" in table and "93.1" in table
    assert (res / "TEST1" / "ramachandran.png").exists()
    assert "Comparative model-quality table" in (res / "REPORT.md").read_text()
