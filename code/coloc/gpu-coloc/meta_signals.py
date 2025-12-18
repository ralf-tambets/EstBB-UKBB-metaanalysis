import argparse, os, gc  
import numpy as np
import pandas as pd

effect_prior    = 0.15
bayes_threshold = 7
window          = 1_000_000

root       = "path/to/meta-analysis/root"
leads_path = "path/to/not_MAF_filtered/02_loci_merged/"

pop_dir_tag = {
    "EstBB":"EstBB","UKBB_AFR":"AFR","UKBB_AMR":"AMR","UKBB_CSA":"CSA",
    "UKBB_EAS":"EAS","UKBB_EUR":"EUR","UKBB_MID":"MID",
}
leads_dict = {
    "EstBB":"EST_loci_merged_full.tsv","AFR":"AFR_loci_merged_full.tsv",
    "AMR":"AMR_loci_merged_full.tsv","CSA":"CSA_loci_merged_full.tsv",
    "EAS":"EAS_loci_merged_full.tsv","EUR":"EUR_loci_merged_full.tsv",
    "MID":"MID_loci_merged_full.tsv",
}

def calc_lbf(beta: pd.Series, se: pd.Series, prior: float) -> pd.Series:
    v = se**2
    r = prior**2 / (prior**2 + v)
    z = beta / se
    return 0.5 * (np.log1p(-r) + r * z**2)

def norm_chr(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip()
    s = s.str.replace(r'^chr', '', regex=True, case=False)
    s = s.str.replace(r'\.0$', '', regex=True)
    s = s.str.replace(r'^(23|x)$', 'X', regex=True, case=False)
    return s.str.upper()

def build_variant(chrom: pd.Series, pos: pd.Series, a0: pd.Series, a1: pd.Series) -> pd.Series:
    return (
        "chr" + chrom.astype(str) + "_" +
        pos.astype(int).astype(str) + "_" +
        a0.astype(str).str.upper() + "_" +
        a1.astype(str).str.upper()
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", type=str, required=True)
    args = parser.parse_args()

    pop = args.pop
    tag = pop_dir_tag[pop]

    pop_path   = f"{root}/{pop}/STEP2"
    summary_fn = f"EST_UK_META/{pop}_summary.tsv"
    os.makedirs("EST_UK_META", exist_ok=True)
    os.makedirs(f"EST_UK_META/{pop}_signals", exist_ok=True)

    leads_df = pd.read_csv(f"{leads_path}{leads_dict[tag]}", sep="\t")
    leads_df.dropna(subset=["POS"], inplace=True)
    leads_df = leads_df[pd.to_numeric(leads_df["POS"], errors="coerce").notna()]
    leads_df["CHR"] = norm_chr(leads_df["CHR"])
    leads_df["POS"] = leads_df["POS"].astype(int)
    leads_df["ALL0"] = leads_df["ALL0"].astype(str).str.upper()
    leads_df["ALL1"] = leads_df["ALL1"].astype(str).str.upper()
    leads_df["variant"] = build_variant(leads_df["CHR"], leads_df["POS"], leads_df["ALL0"], leads_df["ALL1"])

    for metabolite, lead_block in leads_df.groupby("metabolite"):
        print(f"Processing {pop} {metabolite}...")
        parquet_fn = f"{pop_path}/{metabolite}/{metabolite}_{tag}.parquet"
        if not os.path.exists(parquet_fn):
            continue

        df = pd.read_parquet(parquet_fn, columns=["CHROM","GENPOS","ALLELE0","ALLELE1","BETA","SE"])

        se = df["SE"].to_numpy(copy=False)
        beta = df["BETA"].to_numpy(copy=False)
        mask = np.isfinite(se) & np.isfinite(beta) & (se > 0)
        df = df.loc[mask]

        df["CHROM"]  = norm_chr(df["CHROM"])
        df["ALLELE0"] = df["ALLELE0"].astype("string").str.upper().astype("category")
        df["ALLELE1"] = df["ALLELE1"].astype("string").str.upper().astype("category")

        df["lbf"] = calc_lbf(df["BETA"], df["SE"], effect_prior)

        df["variant"] = build_variant(df["CHROM"], df["GENPOS"], df["ALLELE0"], df["ALLELE1"])

        lbf_by_key = df.set_index(["CHROM","GENPOS","ALLELE0","ALLELE1"])["lbf"]

        for _, lead in lead_block.iterrows():
            lead_key = (lead.CHR, int(lead.POS), lead.ALL0, lead.ALL1)  # NEW
            lead_lbf = lbf_by_key.get(lead_key, np.nan)                 # NEW
            if not np.isfinite(lead_lbf) or lead_lbf < bayes_threshold:
                continue

            region = df[
                (df["CHROM"] == lead.CHR) &
                (df["GENPOS"].between(int(lead.POS) - window, int(lead.POS) + window))
            ]
            if region.empty:
                continue

            strength  = float(region["lbf"].max())
            signal_id = f"{pop}_{metabolite}_chr{lead.CHR}:{int(region['GENPOS'].min())}-{int(region['GENPOS'].max())}"

            pd.DataFrame([{
                "signal":          signal_id,
                "chromosome":      lead.CHR,
                "location_min":    int(region["GENPOS"].min()),
                "location_max":    int(region["GENPOS"].max()),
                "signal_strength": strength,
                "lead_variant":    lead.variant,
            }]).to_csv(
                summary_fn, sep="\t", mode="a",
                header=not os.path.exists(summary_fn), index=False,
            )

            mat = pd.DataFrame(region.set_index('variant')['lbf']).T

            mat.to_pickle(os.path.join(f"EST_UK_META/{pop}_signals", f"{signal_id}.pickle"))

        del lbf_by_key, df, se, beta, mask
        gc.collect()

if __name__ == "__main__":
    main()
