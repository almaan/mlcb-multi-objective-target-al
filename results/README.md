# Manuscript outputs

`manuscript/` holds the raw outputs of the run reported in the paper, as written by
`scripts/run_analysis.py` and combined with `scripts/merge_ab_report.py`.
`manuscript/tables/` holds the tables recomputed from them by
`scripts/reproduce_tables.py`.

Conventions used in every file:
- **efficacy** = lethality in the target line (−Chronos; higher = more lethal).
- **toxicity** = lethality in the non-target (contrast) line. The manuscript calls this
  "non-target efficacy" or "non-target activity".
- **tox_source** names the non-target population. `contrast_1/2/3` are ACH-002462,
  ACH-001310 and ACH-000133. `aggregate` is an extra source that the manuscript
  tables do not use.
- All metrics are computed on **true** DepMap values of the nominated genes, not on
  surrogate predictions.

| file | one row per | contents |
|---|---|---|
| `meta.json` | run | configuration (lines, seeds, K, τ, rounds, batch, contrasts, joint GP) |
| `ablation.parquet` | analysis × source × method × line × seed | final nominated-set metrics. `analysis` = `A` (Section 2, Table 1) or `B` (Section 3, Tables 2 / A.2) |
| `rounds.parquet` | source × method × line × seed × round | per-round metrics of the assayed set and of the set nominated after that round (Table A.1) |
| `assayed.parquet` | source × acquisition × line × seed | metrics of the 120 assayed genes |
| `gp_cv.parquet` | line × fold | GP 5-fold cross-validation (R², Spearman, RMSE) |
| `group_counts.parquet` | source × base × operator × line × seed × complex | CORUM complex counts of nominees (diversity barplot) |
| `scatter.parquet` | gene × line × source | true efficacy / toxicity of every candidate gene, GP-predicted toxicity, and a `pick_<method>` flag per Section-2 method (Figure 1, Figure A.1) |
| `multicontrast.parquet` | line × seed × method | exploratory multi-contrast EHVI; not reported in the paper |

Method keys in `ablation.parquet` / `rounds.parquet` map to manuscript labels as follows:

| key | label | key | label |
|---|---|---|---|
| `random` | random | `ehvi` | EHVI |
| `cluster` | cluster | `ehvi_pareto` | EHVI-EHVI |
| `info_div` | info-diverse | `ehvi_trunc` | EHVI·nom |
| `greedy` | greedy | `ehvi_safe` | EHVI·per-round (·RT) |
| `trunc_pred` | greedy·nom | `trunc_known`, `known_safe` | oracle variants (not reported) |
| `greedy_safe` | greedy·per-round (·RT) | | |

In Section B, `<base>+none` is the base method without hedging and `<base>+kdpp_corum`
is the base method with CORUM k-DPP hedging. `<base>+kdpp_string` is an additional
STRING-similarity hedge that the manuscript does not report.
