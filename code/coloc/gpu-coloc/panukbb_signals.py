import argparse, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict
from bisect import bisect_left, insort
import pyarrow.parquet as pq

BAYES_THRESHOLD = 5.0
WINDOW_BP       = 1_000_000
NEGLOG10_P_MIN  = 7.3  # ~5e-8
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

def _columns_for_pop(pop: str) -> tuple[str, str, str, str, str, str]:
    p = pop.upper()
    if p in {"AFR","AMR","CSA","EAS","EUR","MID"}:
        return (f"beta_{p}", f"se_{p}", f"neglog10_pval_{p}", f"af_{p}", f"af_cases_{p}", f"af_controls_{p}")
    if p == "META":
        return ("beta_meta", "se_meta", "neglog10_pval_meta", "af_meta", "af_cases_meta", "af_controls_meta")
    if p in {"META_HQ", "META-HQ", "METAHQ"}:
        return ("beta_meta_hq", "se_meta_hq", "neglog10_pval_meta_hq", "af_meta_hq", "af_cases_meta_hq", "af_controls_meta_hq")
    raise ValueError(f"Unsupported pop '{pop}'")

def _get_or_none(s: pd.Series, col: str):
    val = s.get(col, None)
    return None if pd.isna(val) else val

def process_file(path: str, pop: str, signals_dir: str) -> list[dict]:
    beta_col, se_col, p_col, af_col, af_cases_col, af_control_col = _columns_for_pop(pop)
    need_cols  = ["chr","pos","ref","alt", beta_col, se_col, p_col]
    af_cols    = [af_col, af_cases_col, af_control_col]

    schema_cols = set(pq.read_schema(path).names)

    missing_need = [c for c in need_cols if c not in schema_cols]
    if missing_need:
        print(f"[WARN] Missing required cols in {path}: {missing_need}", file=sys.stderr)
        return []

    present_af_cols = [c for c in af_cols if c in schema_cols]

    df = pd.read_parquet(path, columns=need_cols + present_af_cols)
    if df is None or df.empty:
        print(f"[WARN] Empty: {path}", file=sys.stderr)
        return []

    df.dropna(subset=need_cols, inplace=True)
    if df.empty:
        return []

    if present_af_cols:
        has_af = (
            df[af_col].notna() if af_col in df.columns
            else pd.Series(False, index=df.index)
        )
        has_pair = (
            (df[af_cases_col].notna() & df[af_control_col].notna())
            if (af_cases_col in df.columns and af_control_col in df.columns)
            else pd.Series(False, index=df.index)
        )
        df = df[has_af | has_pair]
        if df.empty:
            return []
        
    df["chr"] = df["chr"].astype("string")

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df.dropna(subset=["pos"], inplace=True)
    df["pos"] = df["pos"].astype(np.int64)

    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)

    # 4) Numeric coercions for effect fields
    df["beta"]   = pd.to_numeric(df[beta_col], errors="coerce")
    df["sebeta"] = pd.to_numeric(df[se_col],   errors="coerce")
    df[p_col]    = pd.to_numeric(df[p_col],    errors="coerce")

    if af_col in df.columns:
        df[af_col] = pd.to_numeric(df[af_col], errors="coerce")
    if af_cases_col in df.columns:
        df[af_cases_col] = pd.to_numeric(df[af_cases_col], errors="coerce")
    if af_control_col in df.columns:
        df[af_control_col] = pd.to_numeric(df[af_control_col], errors="coerce")

    # df["maf"] = df[af_col].apply(lambda x: min(x, 1 - x) if pd.notnull(x) else np.nan)

    df = df.loc[df["sebeta"].notna() & df["beta"].notna() & (df["sebeta"] > 0)]
    if df.empty:
        return []

    base  = os.path.basename(path)
    # trait = Path(base).stem.split(".")[0]      
    # trait = Path(base).stem.removesuffix(".parquet")
    trait = base.removesuffix(".parquet")
    trait_type = base.split("-")[0]               
    prior = 0.15 if trait_type in {"biomarkers", "continuous"} else 0.20
    df["lbf"] = approx_bf_beta(df["beta"], df["sebeta"], prior=prior)

    possible = df[(df["lbf"] >= BAYES_THRESHOLD) & (df[p_col] >= NEGLOG10_P_MIN)].copy()
    if possible.empty:
        return []

    # score so that "smaller is better" for lexsort (use -neglog10p → highest first)
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

        loc_min  = int(region["pos"].min())
        loc_max  = int(region["pos"].max())
        signal_id = f"{pop}_{trait}_{c}:{loc_min}-{loc_max}"
        lead_variant = f"{c}_{int(lead['pos'])}_{lead['ref']}_{lead['alt']}"

        summaries.append({
            "signal":          signal_id,
            "chromosome":      c.removeprefix("chr"),
            "location_min":    loc_min,
            "location_max":    loc_max,
            "signal_strength": lead["lbf"],
            "lead_variant":    lead_variant,
            "af":         _get_or_none(lead, af_col),
            "af_cases":   _get_or_none(lead, af_cases_col),
            "af_control": _get_or_none(lead, af_control_col),
            "-log10p":lead[p_col],
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
