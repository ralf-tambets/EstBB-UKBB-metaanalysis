import argparse, os, sys
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.special import log_ndtr, erfcinv

from collections import defaultdict
from bisect import bisect_left, insort

from pyfaidx import Fasta

fasta = Fasta("hg38.fa")

effect_prior   = 0.15
bayes_threshold= 5
window         = 1_000_000
pval_threshold = 5e-8
root           = "/gpfs/space/projects/genomic_references/summary_stats/MVP"


LOG10 = np.log(10.0)

def ref(chrom: str, pos: int):
    return str(fasta[f"chr{chrom}"][pos-1:pos]).upper()

def neglog10p_from_or_ci(or_, ci_low, ci_high):
    or_ = pd.to_numeric(or_, errors="coerce")
    L   = pd.to_numeric(ci_low, errors="coerce")
    U   = pd.to_numeric(ci_high, errors="coerce")
    zstar = 1.959963984540054  # 95% CI

    beta = np.log(or_)
    se   = (np.log(U) - np.log(L)) / (2.0*zstar)
    z    = beta / se

    logp = np.log(2.0) + log_ndtr(-np.abs(z))
    return -(logp / LOG10)  

def approx_bf_from_log10p(x_log10, F, typ, N, s):
    r0 = 0.15**2 if typ=="quant" else 0.2**2
    V  = 1/(2*N*F*(1-F)) if typ=="quant" else 1/(2*N*F*(1-F)*s*(1-s))
    r  = r0/(r0+V)

    x = np.asarray(x_log10, float)

    z = np.empty_like(x)
    mask = x <= 300
    if np.any(mask):
        z[mask] = np.sqrt(2.0) * erfcinv(10.0**(-x[mask]))
    if np.any(~mask):
        y = x[~mask] * LOG10              
        z0 = np.sqrt(2*y)          
        f  = -y + np.log(z0) + 0.5*np.log(2*np.pi) + 0.5*z0*z0
        fp = 1/z0 + z0
        z[~mask] = z0 - f/fp

    return 0.5*(np.log1p(-r) + r*(z*z))

def approx_bf_beta(beta: pd.Series, se: pd.Series, prior: float) -> pd.Series:
    v = se*se
    r = prior**2 / (prior**2 + v)
    z = beta / se
    return 0.5*(np.log1p(-r) + r*z*z)

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

def pick_leads_by_pval(possible: pd.DataFrame, window: int) -> pd.Series:
    chrom = possible["chrom"].astype("string").to_numpy()
    pos   = possible["pos"].to_numpy(np.int64)
    pval  = pd.to_numeric(possible["pval"], errors="coerce").to_numpy()
    order = np.lexsort((pos, chrom, pval))  

    kept_mask = np.zeros(len(possible), dtype=bool)
    kept_pos  = defaultdict(list) 

    for i in order:
        c, p = chrom[i], int(pos[i])
        arr = kept_pos[c]
        j = bisect_left(arr, p)
        if (j > 0 and p - arr[j-1] <= window) or (j < len(arr) and arr[j] - p <= window):
            continue
        kept_mask[i] = True
        insort(arr, p)
    return pd.Series(kept_mask, index=possible.index)

def make_chr_index(df: pd.DataFrame):
    return {c: sub.sort_values("pos").reset_index(drop=True)
            for c, sub in df.groupby("chrom", sort=False)}

def region_slice(sub: pd.DataFrame, center: int, w: int) -> pd.DataFrame:
    pos = sub["pos"].to_numpy()
    lo, hi = center - w, center + w
    i0 = np.searchsorted(pos, lo, side="left")
    i1 = np.searchsorted(pos, hi, side="right")
    return sub.iloc[i0:i1]

def process_file(path: str, pop: str, signals_dir: str) -> list[dict]:
    df = read_table_auto(path)
    if df is None or df.empty:
        return []

    base = os.path.basename(path)
    suffix = f".{pop}.GIA.dbGaP.txt.gz"
    trait = base[:-len(suffix)] if base.endswith(suffix) else base.split(".")[0]

    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].where(~df["chrom"].isin(["23","23.0"]), "X")
    df["pos"]   = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df.dropna(subset=["pos"], inplace=True)
    df["pos"]   = df["pos"].astype(int)
    df["af"]    = pd.to_numeric(df["af"], errors="coerce")
    df["MAF"]   = df["af"].where(df["af"] <= 0.5, 1 - df["af"])
    df["ref"]   = df["ref"].astype(str)
    df["alt"]   = df["alt"].astype(str)

    df["is_ref"] = [
        False if pd.isna(c) or pd.isna(p) or pd.isna(r) else ref(c, p) == r.upper()
        for c, p, r in zip(df["chrom"], df["pos"].astype(int), df["ref"])
    ]

    df["is_alt"] = [
        False if pd.isna(c) or pd.isna(p) or pd.isna(r) else ref(c, p) == r.upper()
        for c, p, r in zip(df["chrom"], df["pos"].astype(int), df["alt"])
    ]

    df.loc[df["is_alt"], ["ref", "alt"]] = df.loc[df["is_alt"], ["alt", "ref"]].values
    df.dropna(subset=["ref", "alt"], inplace=True)

    if {"beta","sebeta"}.issubset(df.columns):
        df["beta"]   = pd.to_numeric(df["beta"], errors="coerce")
        df["sebeta"] = pd.to_numeric(df["sebeta"], errors="coerce")
        df = df.loc[df["sebeta"].notna() & df["beta"].notna() & (df["sebeta"] > 0)]
        df["lbf"] = approx_bf_beta(df["beta"], df["sebeta"], effect_prior)
    else:
        req = ["pval","MAF","num_samples","num_cases","or","ci"]
        if not all(c in df.columns for c in req):
            print(f"error req not found in columns for {path}")
            return []
        df["pval"]        = pd.to_numeric(df["pval"], errors="coerce")
        df["num_samples"] = pd.to_numeric(df["num_samples"], errors="coerce")
        df["num_cases"]   = pd.to_numeric(df["num_cases"], errors="coerce")
        df["or"]   = df["or"].astype(float)
        df["ci"] = df["ci"].astype(str)
        cc_s = df["num_cases"] / df["num_samples"]
        ci = (df["ci"].astype(str)
                        .str.strip("[](){} ")
                        .str.split(",", n=1, expand=True)
                        .astype(float))
        df[["ci_low","ci_high"]] = ci

        df["-log10(p)"] = neglog10p_from_or_ci(df["or"].to_numpy(),
                                df["ci_low"].to_numpy(),
                                df["ci_high"].to_numpy())

        df["lbf"] = approx_bf_from_log10p(df["-log10(p)"].to_numpy(),
                                df["MAF"].to_numpy(),
                                "cc",
                                df["num_samples"].to_numpy(),
                                cc_s.to_numpy())

    if "pval" in df.columns:
        possible = df[(df["lbf"] >= bayes_threshold) & (df["pval"] <= pval_threshold)].copy()
    else:
        possible = df[df["lbf"] >= bayes_threshold].copy()
    if possible.empty:
        return []

    lead_mask = pick_leads_by_pval(possible, window)
    leads = possible.loc[lead_mask]

    by_chr = make_chr_index(df)
    summaries: list[dict] = []

    for _, lead in leads.iterrows():
        c = lead["chrom"]
        sub = by_chr.get(c)
        if sub is None:
            continue
        region = region_slice(sub, int(lead["pos"]), window)
        if region.empty:
            continue

        strength = float(region["lbf"].to_numpy().max())
        loc_min  = int(region["pos"].min())
        loc_max  = int(region["pos"].max())
        signal_id = f"{pop}_{trait}_chr{c}:{loc_min}-{loc_max}"

        lead_variant = f"chr{c}_{int(lead['pos'])}_{lead['ref']}_{lead['alt']}"

        summaries.append({
            "signal":          signal_id,
            "chromosome":      c,
            "location_min":    loc_min,
            "location_max":    loc_max,
            "signal_strength": strength,
            "lead_variant":    lead_variant,
        })

        vstr = ("chr" + region["chrom"].astype(str) + "_" +
                region["pos"].astype(str) + "_" +
                region["ref"].astype(str) + "_" +
                region["alt"].astype(str))
        mat = pd.DataFrame(region["lbf"].to_numpy()[None, :],
                           columns=vstr.to_numpy(), index=[signal_id])
        mat.to_pickle(os.path.join(signals_dir, f"{signal_id}.pickle"))

    return summaries

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", type=str, required=True)
    args = parser.parse_args()
    pop = args.pop

    out_root = Path("MVP")
    out_root.mkdir(exist_ok=True)
    signals_dir = out_root / f"{pop}_signals"
    signals_dir.mkdir(parents=True, exist_ok=True)
    summary_fn = out_root / f"{pop}_summary.tsv"

    files = [os.path.join(root, x) for x in os.listdir(root)
             if x.endswith(f"{pop}.GIA.dbGaP.txt.gz")]
    if not files:
        print(f"[INFO] No files for pop={pop}", file=sys.stderr)
        return

    all_rows: list[dict] = []
    for f in files:
        rows = process_file(f, pop, str(signals_dir))
        if rows:
            all_rows.extend(rows)

    if all_rows:
        df_sum = pd.DataFrame.from_records(all_rows)
        df_sum.to_csv(summary_fn, sep="\t",
                      mode="a", header=not summary_fn.exists(), index=False)

if __name__ == "__main__":
    main()