"""Command-line entry point: ``python -m ecaf <command>`` (run from the repo root with PYTHONPATH=src)."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from pathlib import Path

from . import __version__, confidence, control, errat, ramachandran, structures, summary, uniprot

ROOT = Path.cwd()


def _load_cfg(path):
    return json.loads(Path(path).read_text())


def _provenance_update(path: Path, **kw):
    prov = json.loads(path.read_text()) if path.exists() else {}
    prov.update(kw)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(prov, indent=2))


def cmd_select(a):
    cfg = _load_cfg(a.config)
    if a.from_tsv:
        cands, prov = uniprot.read_tsv(Path(a.from_tsv)), {"uniprot_query": "offline TSV: " + a.from_tsv}
    else:
        cands, prov = uniprot.fetch_uniprot(cfg)
    uniprot.apply_filters(cands, cfg)
    chosen = uniprot.select(cands, cfg)
    uniprot.write_outputs(cands, chosen, Path(a.out))
    _provenance_update(Path(a.results) / "provenance.json", selection=prov, selection_config=cfg)
    print(f"{len(cands)} records screened, {len(chosen)} selected -> {a.out}/selected_candidates.fasta")
    for c in chosen:
        print(f"  {c.accession}\t{c.chaperone_class}\t{c.organism}\t{c.length} aa")


def cmd_server_jobs(a):
    n = structures.write_server_jobs(Path(a.fasta), Path(a.out))
    print(f"wrote {n} AlphaFold Server jobs -> {a.out} (upload at https://alphafoldserver.com)")


def cmd_fetch_afdb(a):
    seqs = structures.read_fasta(Path(a.fasta)) if a.fasta and Path(a.fasta).exists() else {}
    if a.accessions:
        seqs = {acc: seqs.get(acc, (None, None)) for acc in a.accessions}
    for acc, (_, seq) in seqs.items():
        try:
            meta = structures.fetch_afdb(acc, Path(a.structures) / acc, seq)
            print(f"{acc}: {meta['route']} v{meta.get('model_version')} "
                  f"{'(sequence OK)' if meta.get('sequence_matches_input') else '(CHECK SEQUENCE)'}")
        except Exception as e:
            print(f"{acc}: not available from AlphaFold DB ({e}); use server-jobs/import-server")


def cmd_import(a):
    seqs = structures.read_fasta(Path(a.fasta)) if a.fasta else {}
    seq = seqs.get(a.accession, (None, None))[1]
    fn = structures.import_server if a.kind == "server" else structures.import_af2
    meta = fn(Path(a.src), a.accession, Path(a.structures) / a.accession, seq)
    print(json.dumps(meta, indent=2))


def _candidate_info(sel_tsv: Path) -> dict:
    if not sel_tsv.exists():
        return {}
    with open(sel_tsv, newline="") as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def cmd_analyze(a):
    sdir, vdir, rdir = Path(a.structures), Path(a.validation), Path(a.results)
    info = _candidate_info(Path(a.data) / "selected_candidates.tsv")
    accs = a.accessions or sorted(p.name for p in sdir.iterdir() if (p / "model.pdb").exists())
    cfg = _load_cfg(a.config) if Path(a.config).exists() else {}
    ctrl = cfg.get("procedural_control", {})
    results = []
    for acc in accs:
        model, pae = sdir / acc / "model.pdb", sdir / acc / "pae.json"
        out = rdir / acc
        out.mkdir(parents=True, exist_ok=True)
        meta = json.loads((sdir / acc / "meta.json").read_text()) if (sdir / acc / "meta.json").exists() else {}
        conf = confidence.analyse(model, pae)
        confidence.plot(conf, acc, out / "confidence.png", pae)
        rows = ramachandran.phi_psi_table(model)
        ramachandran.write_csv(rows, out / "ramachandran_residues.csv")
        ramachandran.plot(rows, acc, out / "ramachandran.png")
        pc = ramachandran.run_procheck(model, vdir / acc / "procheck") or ramachandran.load_procheck(acc, vdir)
        er_v = errat.run_errat(model, vdir / acc / "errat")
        er = {"errat_overall": er_v, "errat_source": "local run"} if er_v is not None else errat.load_errat(acc, vdir)
        r = {"accession": acc, "info": info.get(acc, {"protein_name": meta.get("label", "")}), "meta": meta,
             "confidence": conf, "rama_all": ramachandran.summarise(rows),
             "rama_confident": ramachandran.summarise(rows, 70.0), "procheck": pc, "errat": er}
        exp = Path(a.validation) / acc / "experimental"
        if ctrl.get("accession") == acc and exp.exists():
            ref = next(iter(sorted(exp.glob("*.pdb")) + sorted(exp.glob("*.cif"))), None)
            if ref:
                r["control"] = control.compare(model, ref, ctrl.get("chain"))
        (out / "confidence_profile.csv").write_text(
            "resnum,resname,plddt\n" + "\n".join(f"{n},{m},{p}" for n, m, p in zip(*conf["_profile"].values())))
        dump = {k: v for k, v in r.items()}
        dump["confidence"] = {k: v for k, v in conf.items() if k != "_profile"}
        (out / "metrics.json").write_text(json.dumps(dump, indent=2))
        results.append(r)
        print(f"{acc}: mean pLDDT {conf['mean_plddt']}, Rama favoured {r['rama_all']['pct_favoured']}%, "
              f"PROCHECK {'ok' if pc else 'pending'}, ERRAT {'ok' if er else 'pending'}")
    if not results:
        sys.exit("no models found; populate structures/<ACC>/model.pdb first")
    rows = summary.table_rows(results)
    summary.write_table(rows, rdir / "model_quality_table.csv")
    _provenance_update(rdir / "provenance.json", software={
        "ecaf": __version__, "python": platform.python_version(),
        "biopython": __import__("Bio").__version__, "numpy": __import__("numpy").__version__},
        models={acc: r["meta"] for acc, r in zip(accs, results)})
    prov = json.loads((rdir / "provenance.json").read_text())
    summary.markdown_report(results, rows, prov, rdir / "REPORT.md")
    print(f"table -> {rdir}/model_quality_table.csv ; report -> {rdir}/REPORT.md")


def main(argv=None):
    p = argparse.ArgumentParser(prog="ecaf", description=__doc__)
    p.add_argument("--config", default="config/selection.json")
    p.add_argument("--data", default="data")
    p.add_argument("--structures", default="structures")
    p.add_argument("--validation", default="validation")
    p.add_argument("--results", default="results")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("select", help="query UniProtKB/TrEMBL, filter and pick non-redundant candidates")
    s.add_argument("--from-tsv", help="screen a previously downloaded UniProt TSV instead of querying")
    s.add_argument("--out", default="data")
    s.set_defaults(fn=cmd_select)

    s = sub.add_parser("server-jobs", help="write AlphaFold Server batch JSON")
    s.add_argument("--fasta", default="data/selected_candidates.fasta")
    s.add_argument("--out", default="data/alphafold_server_jobs.json")
    s.set_defaults(fn=cmd_server_jobs)

    s = sub.add_parser("fetch-afdb", help="download precomputed AlphaFold DB models where they exist")
    s.add_argument("--fasta", default="data/selected_candidates.fasta")
    s.add_argument("accessions", nargs="*", help="extra/only accessions, e.g. the procedural control Q57733")
    s.set_defaults(fn=cmd_fetch_afdb)

    s = sub.add_parser("import", help="import an AlphaFold Server zip or a local AF v2.3.2 output folder")
    s.add_argument("kind", choices=["server", "af2"])
    s.add_argument("accession")
    s.add_argument("src")
    s.add_argument("--fasta", default="data/selected_candidates.fasta")
    s.set_defaults(fn=cmd_import)

    s = sub.add_parser("analyze", help="confidence + Ramachandran + PROCHECK/ERRAT + summary table/report")
    s.add_argument("accessions", nargs="*")
    s.set_defaults(fn=cmd_analyze)

    a = p.parse_args(argv)
    if a.cmd == "import" and not Path(a.fasta).exists():
        a.fasta = None
    a.fn(a)


if __name__ == "__main__":
    main()
