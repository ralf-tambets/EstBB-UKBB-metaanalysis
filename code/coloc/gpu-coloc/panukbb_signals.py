import argparse, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict
from bisect import bisect_left, insort

BAYES_THRESHOLD = 5.0
WINDOW_BP       = 1_000_000
NEGLOG10_P_MIN  = 7.3
ROOT            = "/path/to/PAN_UKBB_lifted"

def approx_bf_beta(beta: pd.Series, se: pd.Series, prior: float) -> pd.Series:

    v = se * se
    r = (prior * prior) / (prior * prior + v)  
    z = beta / se
    return 0.5 * (np.log1p(-r) + r * z * z)

def pick_leads_by_score(possible: pd.DataFrame, window: int, score_col: str) -> pd.Series:
    chrom = possible["chr"].astype("string").to_numpy()
    pos   = possible["pos"].to_numpy(np.int64)
    score = pd.to_numeric(possible[score_col], errors="coerce").to_numpy()

    order = np.lexsort((pos, chrom, score))

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

def make_chr_index(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    return {c: sub.sort_values("pos").reset_index(drop=True)
            for c, sub in df.groupby("chr", sort=False)}

def region_slice(sub: pd.DataFrame, center: int, w: int) -> pd.DataFrame:
    pos = sub["pos"].to_numpy()
    lo, hi = center - w, center + w
    i0 = np.searchsorted(pos, lo, side="left")
    i1 = np.searchsorted(pos, hi, side="right")
    return sub.iloc[i0:i1]

def _column_triplet_for_pop(pop: str) -> tuple[str, str, str]:
    p = pop.upper()
    if p in {"AFR","AMR","CSA","EAS","EUR","MID"}:
        return (f"beta_{p}", f"se_{p}", f"neglog10_pval_{p}")
    if p == "META":
        return ("beta_meta", "se_meta", "neglog10_pval_meta")
    if p in {"META_HQ", "META-HQ", "METAHQ"}:
        return ("beta_meta_hq", "se_meta_hq", "neglog10_pval_meta_hq")
    raise ValueError(f"Unsupported pop '{pop}'")

def process_file(path: str, pop: str, signals_dir: str) -> list[dict]:
    beta_col, se_col, p_col = _column_triplet_for_pop(pop)
    need_cols = ["chr","pos","ref","alt", beta_col, se_col, p_col]

    df = pd.read_parquet(path, columns=need_cols)
    if df is None or df.empty:
        print(f"[WARN] Empty: {path}", file=sys.stderr)
        return []

    df.dropna(subset=need_cols, inplace=True)
    if df.empty:
        return []

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df.dropna(subset=["pos"], inplace=True)
    df["pos"] = df["pos"].astype(np.int64)

    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)

    df["beta"]   = pd.to_numeric(df[beta_col], errors="coerce")
    df["sebeta"] = pd.to_numeric(df[se_col],   errors="coerce")
    df[p_col]    = pd.to_numeric(df[p_col],    errors="coerce")

    df = df.loc[df["sebeta"].notna() & df["beta"].notna() & (df["sebeta"] > 0)]
    if df.empty:
        return []

    base  = os.path.basename(path)
    trait = Path(base).stem.split(".")[0]         
    trait_type = base.split("-")[0]               
    prior = 0.15 if trait_type in {"biomarkers", "continuous"} else 0.20
    df["lbf"] = approx_bf_beta(df["beta"], df["sebeta"], prior=prior)

    possible = df[(df["lbf"] >= BAYES_THRESHOLD) & (df[p_col] >= NEGLOG10_P_MIN)].copy()
    if possible.empty:
        return []

    possible["_score_asc"] = -possible[p_col].to_numpy()
    possible.dropna(subset=["chr","pos","_score_asc"], inplace=True)
    if possible.empty:
        return []

    lead_mask = pick_leads_by_score(possible, WINDOW_BP, "_score_asc")
    leads = possible.loc[lead_mask]
    if leads.empty:
        return []

    by_chr = make_chr_index(df)
    summaries: list[dict] = []

    for _, lead in leads.iterrows():
        c = lead["chr"]
        sub = by_chr.get(c)
        if sub is None or sub.empty:
            continue

        region = region_slice(sub, int(lead["pos"]), WINDOW_BP)
        if region.empty:
            continue

        strength = float(region["lbf"].to_numpy().max())
        loc_min  = int(region["pos"].min())
        loc_max  = int(region["pos"].max())
        signal_id = f"{pop}_{trait}_{c}:{loc_min}-{loc_max}"
        lead_variant = f"{c}_{int(lead['pos'])}_{lead['ref']}_{lead['alt']}"

        summaries.append({
            "signal":          signal_id,
            "chromosome":      c,
            "location_min":    loc_min,
            "location_max":    loc_max,
            "signal_strength": strength,
            "lead_variant":    lead_variant,
        })

        vstr = (region["chr"].astype(str) + "_" +
                region["pos"].astype(str) + "_" +
                region["ref"].astype(str) + "_" +
                region["alt"].astype(str))

        mat = pd.DataFrame(region["lbf"].to_numpy()[None, :],
                           columns=vstr.to_numpy(), index=[signal_id])

        outp = Path(signals_dir) / f"{signal_id}.pickle"
        outp.parent.mkdir(parents=True, exist_ok=True)
        mat.to_pickle(outp)

    return summaries

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", type=str, required=True,
                        help="One of AFR, AMR, CSA, EAS, EUR, MID, META, META_HQ")
    parser.add_argument("--root", type=str, default=ROOT)
    parser.add_argument("--out",  type=str, default="PanUKBB")
    args = parser.parse_args()
    pop = args.pop

    out_root = Path(args.out)
    out_root.mkdir(exist_ok=True)
    signals_dir = out_root / f"{pop}_signals"
    signals_dir.mkdir(parents=True, exist_ok=True)
    summary_fn = out_root / f"{pop}_summary.tsv"

    files = [str(Path(args.root) / x) for x in os.listdir(args.root) if x.endswith(".parquet")]
    if not files:
        print(f"[INFO] No .parquet files in {args.root} for pop={pop}", file=sys.stderr)
        return

    all_rows: list[dict] = []
    for f in files:
        try:
            rows = process_file(f, pop, str(signals_dir))
            if rows:
                all_rows.extend(rows)
        except Exception as e:
            print(f"[ERROR] {f}: {e}", file=sys.stderr)

    if all_rows:
        pd.DataFrame.from_records(all_rows).to_csv(
            summary_fn, sep="\t", mode="a", header=not summary_fn.exists(), index=False
        )

if __name__ == "__main__":
    main()
