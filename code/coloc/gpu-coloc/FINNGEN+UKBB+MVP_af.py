import os
import numpy as np
import pandas as pd

lead_threshold = 5e-8
root = "/gpfs/space/projects/genomic_references/summary_stats/FinnGen+MVP+UKBB"

summary_fn = f"FinnGen+MVP+UKBB_important_variants.tsv"

for gwas in os.listdir(root):
    df = pd.read_csv(f"/gpfs/space/projects/genomic_references/summary_stats/FinnGen+MVP+UKBB/{gwas}", sep="\t", usecols=['#CHR', 'POS', 'REF', 'ALT','MVP_EUR_af_alt', "fg_af_alt","all_inv_var_meta_p"])

    df = df[df["all_inv_var_meta_p"] <= lead_threshold]

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

    df["trait"] = gwas.removesuffix(".tsv")

    df.to_csv(
        summary_fn,
        sep="\t",
        mode="a",
        header=not os.path.exists(summary_fn),
        index=False,
    )
