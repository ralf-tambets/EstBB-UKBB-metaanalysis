import os
import numpy as np
import pandas as pd

effect_prior = 0.2
bayes_threshold = 5      
lead_threshold = 5e-8
window = 1_000_000     
root = "/path/to/FinnGen+MVP+UKBB" #where the sumstats are located

def calc_lbf(beta: pd.Series, se: pd.Series, prior: float) -> pd.Series:
    v = se**2
    r = prior**2 / (prior**2 + v)
    z = beta / se
    return 0.5 * (np.log1p(-r) + r * z**2)  

summary_fn = f"FinnGen+MVP+UKBB_summary.tsv"
os.makedirs(f"FinnGen+MVP+UKBB_signals", exist_ok=True)

for gwas in os.listdir(root):
    df = pd.read_csv(f"/path/to/{gwas}", sep="\t")

    df["CHR"] = (
        df["#CHR"].astype(int)
        .astype(str)
    )

    df["CHR"] = df["CHR"].where(~df["CHR"].isin(["23", 23]), "X")
    df["POS"] = df["POS"].astype(int)

    df["variant"] = (
        "chr" + df["CHR"] + "_"
        + df["POS"].astype(int).astype(str) + "_"
        + df["REF"].astype(str) + "_"
        + df["ALT"].astype(str)
    )

    df["lbf"] = calc_lbf(df["all_inv_var_meta_beta"], df["all_inv_var_meta_sebeta"], effect_prior)
    
    possible = df.loc[(df["all_inv_var_meta_p"] <= lead_threshold)&(df["lbf"] >= bayes_threshold)].sort_values("all_inv_var_meta_p", ascending=False)

    leads = []
    for _, row in possible.iterrows():
        if not any(
            (abs(lead.POS - row.POS) <= window) and (lead.CHR == row.CHR)
            for lead in leads
        ):
            leads.append(row)

    for _, lead in pd.DataFrame(leads).iterrows():

        region = df[
            (df["CHR"] == lead.CHR)
            & df["POS"].between(lead.POS - window, lead.POS + window)
        ]
        if region.empty:
            continue

        strength   = region["lbf"].max()

        gwas_name = gwas.removesuffix(".tsv") 

        signal_id  = (
            f"FinnGen+MVP+UKBB_{gwas_name}_chr{lead.CHR}:"
            f"{int(region['POS'].min())}-{int(region['POS'].max())}"
        )

        summary_row = pd.DataFrame([{
            "signal":          signal_id,
            "chromosome":      lead.CHR,
            "location_min":    int(region["POS"].min()),
            "location_max":    int(region["POS"].max()),
            "signal_strength": strength,
            "lead_variant":    lead.variant,
        }])
        summary_row.to_csv(
            summary_fn,
            sep="\t",
            mode="a",
            header=not os.path.exists(summary_fn),
            index=False,
        )

        mat = pd.DataFrame(region.set_index('variant')['lbf']).T

        mat.to_pickle(os.path.join(f"FinnGen+MVP+UKBB/FinnGen+MVP+UKBB_signals", f"{signal_id}.pickle"))

