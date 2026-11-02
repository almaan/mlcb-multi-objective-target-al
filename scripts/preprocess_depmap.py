# scripts/preprocess_depmap.py
"""Curate the raw DepMap CRISPR gene-effect CSV into the pipeline's inputs.

  data/raw/depmap/CRISPRGeneEffect.csv  (cell lines x genes, Chronos)
    -> data/processed/depmap/gene_effect.parquet  (genes x cell lines;
       index 'gene' = 'SYMBOL (Entrez)', columns 'cell_line' = ModelID)
    -> data/processed/depmap/panel_all.txt        (Entrez ids, one per line,
       in gene order; the genome-wide candidate panel)

Run from the repository root:  python scripts/preprocess_depmap.py
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from geneal.data.depmap import parse_entrez


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="data/raw/depmap/CRISPRGeneEffect.csv")
    ap.add_argument("--out-dir", default="data/processed/depmap")
    args = ap.parse_args()

    raw = pd.read_csv(args.raw, index_col=0)
    ge = raw.T
    ge.index.name = "gene"
    ge.columns.name = "cell_line"
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ge.to_parquet(out / "gene_effect.parquet")
    (out / "panel_all.txt").write_text("\n".join(str(parse_entrez(g)) for g in ge.index))
    print(f"gene_effect {ge.shape[0]} genes x {ge.shape[1]} cell lines -> {out}")


if __name__ == "__main__":
    main()
