import argparse, os, sys
from pathlib import Path

import numpy as np
import pandas as pd


pval_threshold = 5e-8
root           = "path/to/MVP"

READ_COLS = ["chrom","pos","af","ref","alt","beta","sebeta","pval","num_samples","num_cases","ci", "or"]
DTYPES = {"chrom":"string","ref":"string","alt":"string", "ci":"string"}

def read_table_auto(path: str) -> pd.DataFrame | None:
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
        compression = "gzip" if magic == b"\x1f\x8b" else None
    except Exception:
        compression = None

    try:
        present = pd.read_csv(path, sep="\t", nrows=0, compression=compression).columns
        usecols = [c for c in READ_COLS if c in present]
        return pd.read_csv(path, sep="\t", compression=compression,
                           usecols=usecols, dtype=DTYPES, engine="pyarrow")
    except Exception:
        try:
            return pd.read_csv(path, sep="\t", compression=compression,
                               usecols=usecols, dtype=DTYPES, low_memory=False)
        except Exception as e:
            print(f"[SKIP] {path}: {e}", file=sys.stderr)
            return None

def process_file(path: str, pop: str, summary_file: str) -> list[dict]:
    df = read_table_auto(path)

    df = df[df["pval"] <= pval_threshold]
    
    if df is None or df.empty:
        return []

    base = os.path.basename(path)
    suffix = f".{pop}.GIA.dbGaP.txt.gz"
    trait = base[:-len(suffix)] if base.endswith(suffix) else base.split(".")[0]

    df["trait"] = trait

    df.to_csv(
        summary_file,
        sep="\t",
        mode="a",
        header=not os.path.exists(summary_file),
        index=False,
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", type=str, required=True)
    args = parser.parse_args()
    pop = args.pop

    out_root = Path("MVP_significant_variants")
    out_root.mkdir(exist_ok=True)
    summary_file = out_root / f"{pop}_significant_variants.tsv"

    files = [os.path.join(root, x) for x in os.listdir(root)
             if x.endswith(f"{pop}.GIA.dbGaP.txt.gz")]
    if not files:
        print(f"[INFO] No files for pop={pop}", file=sys.stderr)
        return

    for f in files:
        process_file(f, pop, summary_file)

if __name__ == "__main__":
    main()