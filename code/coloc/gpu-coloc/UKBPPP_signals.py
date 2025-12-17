import argparse
import os
import tarfile
import gzip
import numpy as np
import pandas as pd
from tqdm import tqdm

def calc_lbf(beta: pd.Series, se: pd.Series, prior: float) -> pd.Series:
    v = se**2
    r = prior**2 / (prior**2 + v)
    z = beta / se
    return 0.5 * (np.log1p(-r) + r * z**2)

parser = argparse.ArgumentParser()
parser.add_argument("--pop", type=str, required=True)
args = parser.parse_args()

pop = args.pop

md3 = pd.read_csv("olink_protein_map_3k_v1.tsv", sep="\t")

md3["chr"] = md3["chr"].astype(str)

#African   'Central_South Asian'  'East Asian'            'Middle East' American   Combined              'European (discovery)'
os.makedirs("UKBPPP_gpu", exist_ok=True)

summary_fn = f"UKBPPP_gpu/{pop}_summary.tsv"
os.makedirs(f"UKBPPP_gpu/{pop}_signals", exist_ok=True)

path = f"UKBPPP/UKB-PPP pGWAS summary statistics/{pop}"

signals_dir = f"UKBPPP_gpu/{pop}_signals"
for _, prot in tqdm(md3.iterrows()):
    tar_path = path + "/" + prot['UKBPPP_ProteinID'].replace(":", "_") + "_" + prot['Panel'] + ".tar"

    try:
        with tarfile.open(tar_path, "r:*") as tar:
            for member in tar.getmembers():
                if f"chr{prot['chr']}_" in member.name:
                    with tar.extractfile(member.name) as raw:   
                        with gzip.open(raw, mode="rt") as f:       
                            df = pd.read_csv(f, sep=" ", usecols=["CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "LOG10P", "BETA","SE"])

                            df["CHROM"] = df["CHROM"].astype(str)
                            df = df[(df["CHROM"]==prot["chr"])&(df["GENPOS"]>= int(prot["gene_start"])-1_000_000) & (df["GENPOS"]<= int(prot["gene_end"])+1_000_000)]

                            if df.empty:
                                break

                            df["variant"] =  ("chr" + df["CHROM"] + "_" +
                                df["GENPOS"].astype(int).astype(str) + "_" +
                                df["ALLELE0"].astype(str).str.upper() + "_" +
                                df["ALLELE1"].astype(str).str.upper()
                            )

                            df["maf"] = df["A1FREQ"].where(df["A1FREQ"] <= 0.5, 1 - df["A1FREQ"])

                            df["lbf"] = calc_lbf(df["BETA"], df["SE"], 0.15)

                            signal_id = pop + "_" + member.name.split("/")[-1].removesuffix(".gz") + f"_chr{prot['chr']}:{df['GENPOS'].min()}-{df['GENPOS'].max()}"

                            idx = df["lbf"].idxmax() 
                            signal_strength = df.loc[idx, "lbf"]

                            if signal_strength <= 5:
                                print(member.name, " no strong cis signal")
                                break
                            lead_variant   = df.loc[idx, "variant"] 

                            pd.DataFrame([{
                                "signal":          signal_id,
                                "chromosome":      prot["chr"],
                                "location_min":    df["GENPOS"].min(),
                                "location_max":    df["GENPOS"].max(),
                                "signal_strength": signal_strength,
                                "lead_variant":    lead_variant,
                                "maf": df.loc[idx, "maf"],
                                "-log10p": df.loc[idx, "LOG10P"]
                            }]).to_csv(
                                summary_fn, sep="\t", mode="a",
                                header=not os.path.exists(summary_fn), index=False,
                            )

                            mat = pd.DataFrame(df["lbf"].to_numpy()[None, :],
                                columns=df["variant"].to_numpy(), index=[signal_id])
                            mat.to_pickle(os.path.join(signals_dir, f"{signal_id}.pickle"))
                            break
    except Exception as e:
        print(prot['UKBPPP_ProteinID'],e)