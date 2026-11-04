# Nominating Selective and Diverse Targets: Multi-Objective Active Learning and Portfolio Hedging

Reproduction code for the MLCB 2026 paper. The framework treats target nomination as
a **selective multi-objective** problem: maximize efficacy in a target cell population
while minimizing activity in non-target populations, under a fixed experimental
budget. It then **hedges** the nominated cohort against mechanistic redundancy using
prior biological knowledge (CORUM protein complexes).

All experiments use the DepMap CRISPR gene-effect resource as an offline oracle.
Ground-truth lethality is hidden and revealed only when the active-learning loop
"assays" a batch of genes.

This repository has two uses:

1. **Check the numbers.** `results/manuscript/` holds the raw outputs of the run
   reported in the paper. `scripts/reproduce_tables.py` recomputes every table from
   them, with no data download and no GPU, in about a minute.
2. **Re-run everything.** The rest of this README covers downloading the data,
   preprocessing it, building the embeddings and similarity graphs, and running the
   active-learning experiments from scratch.

---

## Repository layout

```
src/geneal/              library: GP surrogates, EHVI, baselines, hedging, metrics, report
scripts/
  download_data.py       fetch + MD5-verify the raw inputs
  preprocess_depmap.py   DepMap CSV -> gene_effect.parquet + genome-wide panel
  build_graph_caches.py  CORUM + STRING edge lists and CORUM membership
  embed_pubmedbert.py    PubMedBERT gene-text embeddings
  run_analysis.py        the experiments (Section 2 = "A", Section 3 = "B")
  aggregate_ablation.py  combine per-cell-line shards into one result + report
  merge_ab_report.py     combine a separate A run and B run into one result + report
  launch_ablation_sweep.sh  optional SLURM launcher (one job per target cell line)
  reproduce_tables.py    recompute the manuscript tables from a result directory
jobs/                    SLURM job scripts used by the launcher
inputs/                  gene-text snapshot used for the embeddings (see Data)
results/manuscript/      outputs of the manuscript run + recomputed tables
tests/                   pytest suite
```

---

## Installation

```bash
conda create -n geneal -c conda-forge python=3.12
conda activate geneal
pip install -e ".[embed,dev]"
```

`embed` installs `mygene` and `sentence-transformers` (only needed to rebuild the
embeddings). Everything runs on CPU. Versions used for the manuscript run: Python
3.12, numpy 2.5, pandas 3.0, scipy 1.18, scikit-learn 1.9, torch 2.11, gpytorch 1.15,
networkx 3.6, sentence-transformers 5.5, mygene 3.2.

```bash
pytest          # unit tests
```

---

## 1. Check the manuscript numbers (no download needed)

```bash
python scripts/reproduce_tables.py          # reads results/manuscript/
```

It prints, and writes as CSV to `results/manuscript/tables/`:

| output | manuscript |
|---|---|
| `table1.csv` | Table 1, selectivity (mean ± SEM over 12 lines × 2 seeds × 3 contrasts) |
| `tableA1.csv` | Table A.1, change in selectivity ratio per round (mean ± 95% CI over lines × seeds; contrasts averaged within each) |
| `table2_corum.csv` | Table 2, CORUM hedging evaluated on CORUM (mean ± SEM) |
| `tableA2_string.csv` | Table A.2, CORUM hedging evaluated on STRING (mean ± SEM) |
| `in_text.csv` | GP 5-fold CV R² / Spearman, and the target vs non-target efficacy correlation |

The result files are described in [`results/README.md`](results/README.md).

---

## 2. Data

| input | version | how |
|---|---|---|
| DepMap `CRISPRGeneEffect.csv`, `Model.csv` | DepMap Public **26Q1** | **manual**: download both from the [DepMap portal](https://depmap.org/portal/data_page/?tab=allData) into `data/raw/depmap/` |
| STRING human `protein.info` + `protein.links` | **v12.0** | automatic |
| CORUM human complexes (`humanComplexes.txt`) | **5.3** | tried automatically through the official CORUM API; otherwise download from [CORUM](https://mips.helmholtz-muenchen.de/corum/) into `data/corum_dl/` |
| gene text (symbol, name, summary) | MyGene.info snapshot | shipped in `inputs/gene_text_mygene.json.gz` |

```bash
python scripts/download_data.py
```

This fetches what it can and checks every file against the MD5 of the copy used for
the manuscript. It prints manual-download instructions for anything missing.
CORUM only serves its current release, so a newer CORUM version fails the MD5 check.
The pipeline still runs, but the hedging results may differ slightly.

MyGene.info annotations change over time. The exact gene texts embedded for the paper
are therefore shipped in `inputs/`, and `embed_pubmedbert.py` reads them by default.

---

## 3. Preprocess, embed, build graphs

Run all commands from the repository root.

```bash
# DepMap CSV -> genes x cell-lines parquet + genome-wide panel (18,531 genes)
python scripts/preprocess_depmap.py

# CORUM / STRING edge lists + CORUM complex membership
python scripts/build_graph_caches.py

# PubMedBERT (NeuML/pubmedbert-base-embeddings) gene-text embeddings, 18,531 x 768
python scripts/embed_pubmedbert.py
```

---

## 4. Run the experiments

Fixed configuration (manuscript):

| setting | value |
|---|---|
| target cell lines | 12 (fewest missing values, contrast lines excluded) |
| non-target (contrast) lines | ACH-002462 (RPE1-ss48), ACH-001310 (HA1E), ACH-000133 (Hs 729) |
| seeds | 0, 1 |
| initial set / rounds / batch | 40 / 8 / 10 (120 genes assayed) |
| nominated cohort K | 30 |
| non-target effect | absolute knockout effect \|t\| in the contrast line (`--tox-transform abs`, the default) |
| threshold τ | 0.5 on \|t\| (absolute Chronos scale; two-sided on the signed effect) |
| surrogate | joint multi-task GP (Matérn-5/2) over (efficacy, non-target effect) |
| candidate pool | full genome for both sections |

**Single machine** (runs all 12 lines sequentially; slow, many hours on CPU):

```bash
python scripts/run_analysis.py \
  --joint-gp --n-cell-lines 12 --tox-transform abs --analyses both \
  --contrast-lines ACH-002462 ACH-001310 ACH-000133 \
  --seeds 0 1 --K 30 --n-initial 40 --n-rounds 8 --batch 10 --tau 0.5 --export-figs
# -> res/runs_ablation/<timestamp>/  (ablation.parquet, rounds.parquet, ..., report.html, figs/)
```

To run one target line (for example as a parallel job), add `--cell-lines <ModelID>
--no-report --run-name shard_<i>` and keep the shard directories under one parent. Then
combine them with `python scripts/aggregate_ablation.py <parent_dir>`.

**SLURM** (one array task per target line, then a dependent aggregation job):

```bash
cp slurm_env.template.sh slurm_env.sh     # fill in partition / account / conda paths
CONTRASTS=3 CONTRAST_IDS="ACH-002462 ACH-001310 ACH-000133" JOINT=1 SEEDS="0 1" \
  ANALYSES=both TOX_TRANSFORM=abs PANEL_A=none PANEL_B=none EXPORTFIGS=1 \
  TAG=mlcb bash scripts/launch_ablation_sweep.sh 12
```

Then recompute the tables from the new run:

```bash
python scripts/reproduce_tables.py --run res/runs_ablation/<run_dir>
```

Figures are written to `<run_dir>/figs/` (`--export-figs` / `EXPORTFIGS=1`):
`gene_cloud_contrast.{pdf,png}` is Figure 1 and `ds_jointeff_contrast.{pdf,png}` is
Figure A.1.

### How the manuscript result was produced

All manuscript tables and figures come from **one** run of the SLURM command above
(both sections in the same run, `--tox-transform abs`). Section 3 reuses Section 2's
acquisition histories, so the unhedged rows of Table 2 are the same nominations as the
corresponding rows of Table 1. `results/manuscript/` holds that run's outputs unchanged.

The GP surrogates seed their own (forked) torch random number generator, and every
acquisition and nomination step uses an explicit, seeded numpy generator, so a run is
deterministic for given inputs. Re-running on the same inputs with the library versions
listed above should reproduce `results/manuscript/`; floating-point differences across
hardware or library versions can change individual picks and hence the last digits.

---

## Method summary (where things live)

| component | code |
|---|---|
| non-target effect definition (\|t\| or signed) | `src/geneal/data/selective.py` (`toxicity_vector`) |
| exact GP (Matérn-5/2) and joint multi-task GP (seeded, deterministic fits) | `src/geneal/models/surrogate.py`, `models/multitask.py` |
| Monte-Carlo EHVI, hypervolume, Pareto fronts | `src/geneal/models/multiobjective.py` |
| EHVI active-learning loop (kriging-believer batches; per round: top-150 shortlist by μ+σ summed over objectives, reference point and front from assayed labels) | `src/geneal/runner/bivariate.py` |
| greedy / random / cluster / info-diverse acquisition, per-round filter | `src/geneal/runner/ablation.py` (`acquire`, `_safe_candidates`) |
| nomination: τ filter, EHVI-EHVI hypervolume selection, k-DPP hedging | `src/geneal/runner/ablation.py` (`nominate`), `models/hedged_selection.py` |
| concentration, robustness, risk, N_eff | `src/geneal/metrics/portfolio.py` |
| HTML / LaTeX report | `src/geneal/report/ablation_report.py` |

---

## License

Code: Apache License 2.0, Copyright 2026 Genentech, Inc. (see [`LICENSE`](LICENSE)).

### Data licenses

DepMap data are released under CC BY 4.0, STRING under CC BY 4.0, and CORUM under its
own terms. The files in `results/manuscript/` contain DepMap-derived gene-effect
values and are redistributed under the DepMap terms. Please cite the original
resources when using them.
