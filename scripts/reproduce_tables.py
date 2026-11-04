# scripts/reproduce_tables.py
"""Recompute the manuscript's tables and in-text statistics from a run directory.

Defaults to the shipped manuscript outputs (results/manuscript/), so the numbers
in the paper can be checked without re-running the active-learning loop. Point
--run at a fresh run (res/runs_ablation/<run>/) to compare a re-run.

  Table 1    selectivity (Section 2): mean +/- SEM over 12 lines x 2 seeds x 3 contrasts
  Table A.1  change in selectivity ratio per round: mean +/- 95% CI over lines x seeds
             (the 3 contrasts are averaged within each line x seed, n = 24)
  Table 2    CORUM hedging, evaluated on CORUM:  mean +/- SEM
  Table A.2  CORUM hedging, evaluated on STRING: mean +/- SEM
  In text    GP 5-fold CV (R^2, Spearman; mean +/- 95% CI over lines x folds) and the
             target / non-target efficacy correlation over all candidate genes.

Writes one CSV per table to --out and prints them.
Run from the repository root:  python scripts/reproduce_tables.py
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

FLOOR = 0.05   # T/NT ratio denominator floor (same as the report)

T1 = [("random", "random"), ("cluster", "cluster"), ("info_div", "info-diverse"),
      ("greedy", "greedy"), ("trunc_pred", "greedy·nom"), ("greedy_safe", "greedy·per-round"),
      ("ehvi", "EHVI"), ("ehvi_pareto", "EHVI-EHVI"), ("ehvi_trunc", "EHVI·nom"),
      ("ehvi_safe", "EHVI·per-round")]
T2_BASES = [("greedy", "greedy"), ("trunc_pred", "greedy·nom"), ("greedy_safe", "greedy·RT"),
            ("ehvi", "EHVI"), ("ehvi_trunc", "EHVI·nom"), ("ehvi_safe", "EHVI·RT")]


def _fmt(m, s):
    return f"{m:.3f}±{s:.3f}"


def _table(df, key, order, cols, spread="sem"):
    rows = []
    for k, label in order:
        g = df[df[key] == k]
        row = {"method": label, "n": len(g)}
        for c, name in cols:
            x = g[c].dropna().to_numpy(float)
            se = x.std(ddof=1) / np.sqrt(len(x))
            row[name] = _fmt(x.mean(), se if spread == "sem" else 1.96 * se)
        rows.append(row)
    return pd.DataFrame(rows)


def table1(ab):
    a = ab[(ab.analysis == "A") & ab.tox_source.str.startswith("contrast")].copy()
    a["ratio"] = a.mean_efficacy / a.mean_toxicity.clip(lower=FLOOR)
    return _table(a, "method", T1, [
        ("mean_efficacy", "Target Efficacy"), ("mean_toxicity", "Non-Target Efficacy"),
        ("ratio", "T/NT ratio"), ("useful_efficacy", "Realizable target efficacy"),
        ("n_safe", "# safe (of K)")])


def table_a1(rounds):
    d = rounds[rounds.tox_source.str.startswith("contrast")].copy()
    d["ratio"] = d.nom_mean_efficacy / d.nom_mean_toxicity.clip(lower=FLOOR)
    rows = []
    for k, label in T1:
        piv = d[d.method == k].pivot_table(index=["cell_line", "seed"], columns="round",
                                           values="ratio")   # averages the contrasts
        r0 = piv.columns.min()
        row = {"method": label, "n": len(piv)}
        for r in [c for c in piv.columns if c % 2 == 0]:
            x = (piv[r] - piv[r0]).to_numpy(float)
            row[f"Δ r{int(r)}"] = _fmt(x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x)))
        rows.append(row)
    return pd.DataFrame(rows)


def _hedge_order():
    order = []
    for op, suffix in [("none", ""), ("kdpp_corum", " + CORUM")]:
        for base, label in T2_BASES:
            order.append((f"{base}+{op}", label + suffix))
    # manuscript groups greedy rows before EHVI rows
    return sorted(order, key=lambda x: (not x[1].startswith("greedy"), "CORUM" in x[1]))


def table2(ab, graph):
    b = ab[(ab.analysis == "B") & ab.tox_source.str.startswith("contrast")]
    sfx = "" if graph == "corum" else "_string"
    return _table(b, "method", _hedge_order(), [
        (f"concentration{sfx}", "Conc."), (f"robustness{sfx}", "Robust."),
        (f"n_pathways{sfx}", "Avg. Units"), (f"risk_{graph}", "Risk"),
        (f"neff_{graph}", "Neff"), ("mean_efficacy", "T. Eff."),
        ("mean_toxicity", "NT. Eff.")])


def in_text(run):
    cv = pd.read_parquet(run / "gp_cv.parquet")
    out = {}
    for c in ["r2", "spearman"]:
        x = cv[c].to_numpy(float)
        out[f"GP CV {c}"] = _fmt(x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x)))
    sc = run / "scatter.parquet"
    if sc.exists():
        s = pd.read_parquet(sc, columns=["efficacy", "toxicity", "tox_source"])
        s = s[s.tox_source.str.startswith("contrast")]
        out["target vs non-target efficacy r"] = f"{np.corrcoef(s.efficacy, s.toxicity)[0, 1]:.4f}"
    return pd.DataFrame([{"statistic": k, "value": v} for k, v in out.items()])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="results/manuscript")
    ap.add_argument("--out", default=None, help="CSV output dir (default <run>/tables)")
    args = ap.parse_args()
    run = Path(args.run)
    out = Path(args.out) if args.out else run / "tables"
    out.mkdir(parents=True, exist_ok=True)

    ab = pd.read_parquet(run / "ablation.parquet")
    tables = {"table1": table1(ab),
              "tableA1": table_a1(pd.read_parquet(run / "rounds.parquet")),
              "table2_corum": table2(ab, "corum"),
              "tableA2_string": table2(ab, "string"),
              "in_text": in_text(run)}
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        for name, t in tables.items():
            t.to_csv(out / f"{name}.csv", index=False)
            print(f"\n== {name} ==")
            print(t.to_string(index=False))
    print(f"\nCSVs written to {out}/")


if __name__ == "__main__":
    main()
