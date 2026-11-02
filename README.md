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
| threshold τ | 0.5 (absolute, on the Chronos lethality scale) |
| surrogate | joint multi-task GP (Matérn-5/2) over (efficacy, non-target effect) |
| candidate pool | full genome for both sections |

**Single machine** (runs all 12 lines sequentially; slow, many hours on CPU):

```bash
python scripts/run_analysis.py \
  --joint-gp --n-cell-lines 12 \
  --contrast-lines ACH-002462 ACH-001310 ACH-000133 \
  --seeds 0 1 --K 30 --n-initial 40 --n-rounds 8 --batch 10 --tau 0.5
# -> res/runs_ablation/<timestamp>/  (ablation.parquet, rounds.parquet, ..., report.html)
```

To run one target line (for example as a parallel job), add `--cell-lines <ModelID>
--no-report --run-name shard_<i>` and keep the shard directories under one parent. Then
combine them with `python scripts/aggregate_ablation.py <parent_dir>`.

**SLURM** (one array task per target line, then a dependent aggregation job):

```bash
cp slurm_env.template.sh slurm_env.sh     # fill in partition / account / conda paths
CONTRASTS=3 CONTRAST_IDS="ACH-002462 ACH-001310 ACH-000133" JOINT=1 SEEDS="0 1" \
  TAG=mlcb bash scripts/launch_ablation_sweep.sh 12
```

Then recompute the tables from the new run:

```bash
python scripts/reproduce_tables.py --run res/runs_ablation/<run_dir>
```

### How the manuscript result was produced

Section 2 (Table 1, Table A.1) and Section 3 (Tables 2, A.2) come from two runs of the
same code and configuration: the first with `--analyses both`, the second with
`--analyses b` (genome-wide hedging). `scripts/merge_ab_report.py` combined them
without recomputing anything:

```bash
python scripts/merge_ab_report.py --a-run <run_1> --b-run <run_2> --out <merged_dir>
```

The two runs fit their GP surrogates independently, and GP training is not bit-exact
across runs. As a result, the unhedged rows of Table 2 differ slightly from the
corresponding rows of Table 1 (for example, EHVI target efficacy 1.395 vs 1.259). A
fresh re-run should match the manuscript within that run-to-run variation, not
digit for digit.

---

## Method summary (where things live)

| component | code |
|---|---|
| exact GP (Matérn-5/2) and joint multi-task GP | `src/geneal/models/surrogate.py`, `models/multitask.py` |
| Monte-Carlo EHVI, hypervolume, Pareto fronts | `src/geneal/models/multiobjective.py` |
| EHVI active-learning loop (kriging-believer batches, top-150 UCB shortlist per round) | `src/geneal/runner/bivariate.py` |
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
