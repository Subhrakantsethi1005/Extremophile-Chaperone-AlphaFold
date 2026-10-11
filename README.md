# Extremophile Chaperone AlphaFold

**Structural prediction and conformational validation of uncharacterized extremophile chaperone proteins using AlphaFold.**
UBT521 *A.I. in Biotechnology*, Group 15, Thapar Institute of Engineering and Technology.

This repository implements the workflow from the group's research proposal as a reproducible
Python pipeline (`ecaf`). Each stage matches a section of the proposal:

| Proposal section | Stage | Command / module |
|---|---|---|
| 3.2 Materials and data sources | Query UniProtKB/TrEMBL, apply the exclusion criteria, pick 3–5 non-redundant candidates | `select` · `uniprot.py` |
| 3.3 Structure prediction | AlphaFold Server batch JSON, import Server or local AlphaFold v2.3.2 output, or download AlphaFold DB models | `server-jobs`, `import`, `fetch-afdb` · `structures.py` |
| 3.4 Primary confidence analysis | Mean/median pLDDT, confidence bands, low-confidence segments, PAE between confident regions | `analyze` · `confidence.py` |
| 3.5 Conformational validation | PROCHECK v3.5 (local run or SAVES `.sum`), a built-in Ramachandran screen with residue mapping, ERRAT (local run or SAVES score) | `analyze` · `ramachandran.py`, `errat.py` |
| 3.6 Analysis plan and decision criteria | Comparative table and rule-based evidence tier. Metric disagreements are listed, not averaged | `analyze` · `summary.py` |
| 3.7 Controls and reproducibility | Known-structure control (*M. jannaschii* Hsp16.5, Q57733 / PDB 1SHS): CA-RMSD and TM-score; `provenance.json` with release, dates and versions | `control.py` |

## Setup

```bash
pip install -r requirements.txt
export PYTHONPATH=src          # or: pip install -e .
python -m pytest               # offline test suite (10 tests)
```

## Running the project

Run every command from the repository root.

### 1. Curate candidates (3.2)

```bash
python -m ecaf select
```

This queries UniProtKB with `reviewed:false`, `fragment:false`, no experimental 3D structure, no
protein-level evidence, chaperone terms (keyword KW-0143, *chaperone*, *chaperonin*, *heat shock*,
*prefoldin*, *thermosome*) and the extremophile genera listed in `config/selection.json`. It then:

* assigns a chaperone class from the protein name and excludes records whose length falls outside
  that class's range (for example, sHSP 110–220 aa and Hsp60/thermosome 500–600 aa)
* excludes records with more than 1 % ambiguous residues or no initial Met (a possible fragment)
* selects up to `target_n` = 5 candidates, diversity first (different classes, organisms and stress
  types), and rejects any sequence with more than 40 % identity to one already selected

Outputs:

* `data/selected_candidates.fasta`
* `data/selected_candidates.tsv`
* `data/screening_log.tsv`: every record screened and why it was kept or excluded
* `results/provenance.json`: UniProt release, release date, access date and the exact query

To screen a UniProt TSV you downloaded in a browser instead, use `--from-tsv file.tsv`. The file needs
the columns `accession, id, protein_name, organism_name, organism_id, length, annotation_score,
protein_existence, keyword, sequence`.

### 2. Obtain AlphaFold models (3.3)

Use any of these routes. Each one writes `structures/<ACC>/model.pdb`, plus `pae.json` and `meta.json`.

* **AlphaFold Server** (the proposal's fallback route):
  ```bash
  python -m ecaf server-jobs        # -> data/alphafold_server_jobs.json; upload at alphafoldserver.com
  python -m ecaf import server <ACC> path/to/fold_<acc>.zip
  ```
  The import keeps `model_0` (the top-ranked model), converts it to PDB format and keeps PAE, pTM
  and the ranking score.
* **Local AlphaFold v2.3.2** (default settings, one FASTA per protein):
  ```bash
  python -m ecaf import af2 <ACC> path/to/af2_output/<ACC>/
  ```
  The import reads `ranked_0.pdb` and `ranking_debug.json`. It also reads PAE from `result_<top>.pkl`
  when a pTM model was used.
* **AlphaFold DB precomputed models** (many TrEMBL entries already have one):
  ```bash
  python -m ecaf fetch-afdb          # all selected candidates
  python -m ecaf fetch-afdb Q57733   # procedural control
  ```
  If the AFDB sequence differs from the curated FASTA, the import records a warning. Re-predict that
  candidate with one of the other routes.

### 3. External validation: PROCHECK v3.5 and ERRAT (3.5)

Both tools are licensed or server-hosted programs, so they are not bundled here.

* **SAVES server** (https://saves.mbi.ucla.edu): upload `structures/<ACC>/model.pdb` and run
  PROCHECK and ERRAT. Then either:
  * save the PROCHECK `*.sum` file to `validation/<ACC>/procheck/` and the ERRAT log to
    `validation/<ACC>/errat/`, or
  * type the numbers into `validation/procheck_manual.csv` (`core, allowed, generous, disallowed` %)
    and `validation/errat_manual.csv` (`overall_quality_factor`, plus any flagged regions).
* **Local installs**: if `procheck` is on your `PATH`, `analyze` runs it automatically. For ERRAT,
  set a command template, for example `export ECAF_ERRAT_CMD="errat {stem}"`.

For the procedural control, put the experimental structure (for example `1SHS.pdb`) in
`validation/Q57733/experimental/`.

### 4. Analyse and report (3.4–3.7)

```bash
python -m ecaf analyze
```

For each candidate this writes the following to `results/<ACC>/`:

* `confidence.png`: pLDDT profile and PAE heatmap
* `confidence_profile.csv`
* `ramachandran.png`
* `ramachandran_residues.csv`: φ/ψ, region and pLDDT for every residue
* `metrics.json`

It also writes these to `results/`:

* `model_quality_table.csv`: one row per candidate with accession, organism, length, mean pLDDT,
  residues with pLDDT < 70, Ramachandran favoured and outlier %, PROCHECK core and disallowed %,
  ERRAT score, evidence tier and major concerns
* `REPORT.md`: the full results write-up with figures
* `provenance.json`

## Decision criteria (3.6)

Candidates are not ranked by a single metric. `summary.py` assigns one of three evidence tiers:

* **stronger**: mean pLDDT ≥ 80, at least 95 % of confident residues (pLDDT ≥ 70) in favoured
  Ramachandran regions, no Ramachandran outliers among confident residues, PROCHECK core ≥ 90 % with
  disallowed ≤ 1 %, and ERRAT ≥ 85
* **weaker**: mean pLDDT < 70, ERRAT < 70, or PROCHECK disallowed > 2 %
* **moderate**: everything else

Concerns are reported explicitly. Ramachandran outliers that fall only in low-pLDDT regions are
marked as expected. Conflicts are flagged, for example high pLDDT with weak ERRAT. Until PROCHECK
and ERRAT results are supplied, the tier is labelled *provisional*. The thresholds are conventions:
they are listed in `summary.CRITERIA` and copied into the report.

## About the built-in Ramachandran screen

The screen classifies every residue against the Lovell et al. (2003) Top500 φ/ψ distributions
(Favoured ≥ 2 % contour; Allowed ≥ 0.05 % for general residues and ≥ 0.2 % for Gly/Pro/pre-Pro).
It is a MolProbity-style cross-check that maps each unusual residue back to the sequence and to its
pLDDT. It does **not** replace PROCHECK's core/allowed/generous/disallowed statistics, which the
proposal specifies. Report PROCHECK as the primary Ramachandran result.

Tested on the 1.7 Å crystal structure 1A8O (`tests/data`), the screen gives 98.3 % favoured and
0 % outliers, as expected for a well-refined experimental structure.

## Limitations

AlphaFold models are computational hypotheses. pLDDT is a local confidence estimate, not
experimental accuracy. PROCHECK, ERRAT and the Ramachandran screen assess geometry and interaction
patterns only. No functional or experimental structural claims are made.

## Layout

```
config/selection.json   selection criteria, extremophile genera, chaperone classes, control
src/ecaf/               pipeline package (python -m ecaf ...)
data/                   curated FASTA/TSV and screening log (generated)
structures/<ACC>/       model.pdb, pae.json, meta.json (generated or imported)
validation/             PROCHECK/ERRAT outputs or manual CSVs, experimental control structure
results/                per-candidate figures/metrics, comparative table, REPORT.md
tests/                  offline tests (data: PDB 1A8O)
```
