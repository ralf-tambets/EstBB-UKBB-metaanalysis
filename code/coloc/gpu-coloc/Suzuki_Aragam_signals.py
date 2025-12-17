import argparse
import os
import numpy as np
import pandas as pd
from tqdm import tqdm

def approx_bf_estimates(variant, z, V, type, suffix=None, sdY=1, effect_priors={'quant': 0.15, 'cc': 0.2}):
    sd_prior = effect_priors['quant'] * sdY if type == "quant" else effect_priors['cc']
    
    r = sd_prior**2 / (sd_prior**2 + V)
    lABF = 0.5 * (np.log(1 - r) + (r * z**2))
    
    ret = pd.DataFrame({"variant": variant, 'lbf': lABF})
    
    if suffix is not None:
        ret.columns = [f"{col}.{suffix}" for col in ret.columns]

    ret.loc[:, 'lbf'] = pd.to_numeric(ret['lbf'], errors='coerce')

    ret_mat = pd.DataFrame(ret.set_index('variant')['lbf']).T

    return ret_mat

parser = argparse.ArgumentParser(description="Seperate and map GWAS signals")

parser.add_argument("--summary", type=str, required=True, help="Path where to write summary, e.g., 'met_summary.tsv'.")
parser.add_argument("--output", type=str, required=True, help="Dir where to write the signals, e.g., 'diasease_signals'.")
parser.add_argument("--input_table", type=str, required=True, help="table of inputs, e.g., 'gwas_table.tsv'.")

args = parser.parse_args()

summary_file_path = args.summary 
signals_dir = args.output
os.makedirs(signals_dir, exist_ok=True)

input_table = pd.read_csv(args.input_table, sep="\t")

for gwas in tqdm(input_table.itertuples(index=False), desc=f"Processing GWAS's"):

    lead_variants = pd.read_csv(gwas.leads, sep="\t")

    diasease = gwas.GWAS

    # lead_variants = lead_variants[lead_variants["diasease"]==diasease]

    summary_rows = []

    diasease_df = pd.read_parquet(
        gwas.GWAS_path
    )

    for region in tqdm(lead_variants.itertuples(index=False), desc=f"Processing {diasease} regions"):
        chrom = region.CHROM
        region_start = region.GENPOS - 1_000_000
        region_end = region.GENPOS + 1_000_000

        region_diasease = diasease_df[diasease_df["CHROM"] == chrom]

        chrom_label = "X" if chrom == 23 else str(chrom)

        region_signal = region_diasease[(region_diasease["GENPOS"] >= region_start) & (region_diasease["GENPOS"] <= region_end)]
        # region_signal.loc[:, "MAF"] = np.minimum(
        #     region_signal["AC"] / (2 * region_signal["N"]),
        #     1 - region_signal["AC"] / (2 * region_signal["N"])
        # )
        # region_signal = region_signal[region_signal["MAF"] >= 0.01]

        # if region_signal.empty:
        #     continue

        region_signal["variant"] = (
            "chr" + chrom_label + "_" +
            region_signal["GENPOS"].astype(int).astype(str) + "_" +
            region_signal["REF"].astype(str) + "_" +
            region_signal["ALT"].astype(str)
        )

        region_signal["z"] = region_signal["BETA"] / region_signal["SE"]
        region_signal["V"] = region_signal["SE"] ** 2

        mat = approx_bf_estimates(region_signal["variant"],region_signal["z"], region_signal["V"], "cc") 
        
        signal_strength = mat.T["lbf"].max()
        variant_id = mat.T["lbf"].idxmax()

        if signal_strength < 5:
            continue

        signal = f"{diasease}_chr{chrom_label}:{region_start}-{region_end}"

        summary_data = pd.DataFrame([{
            'signal': signal,
            'chromosome': chrom_label,
            'location_min': region_start,
            'location_max': region_end,
            'signal_strength': signal_strength,
            'lead_variant': variant_id
        }])

        header_needed = not os.path.exists(summary_file_path)  
        summary_data.to_csv(summary_file_path, sep='\t', mode='a', header=header_needed, index=False)

        output_file_name = f"{signal}.pickle"
        output_file_path = os.path.join(signals_dir, output_file_name)

        mat.to_pickle(output_file_path)

