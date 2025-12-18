import os
import pandas as pd
import duckdb
from tqdm import tqdm

PATH = ... #Path to meta_EUR metabolite parquets
OUTDIR = "metabolite_variants"
os.makedirs(OUTDIR, exist_ok=True)

lead = pd.read_csv("path/to/meta_EUR_loci_merged_full.tsv", sep="\t", usecols=["snp_alternate"])
lead["variant"] = "chr" + lead["snp_alternate"].astype(str)
lead = lead[["variant"]].drop_duplicates()

con = duckdb.connect()
con.register("lead", lead) 

for file in tqdm(os.listdir(PATH)):
    if not file.endswith(".parquet"):
        continue

    metabolite = file.removesuffix("_EstBB_UKBB_EUR_metaanalysis.parquet")
    parquet_path = f"{PATH}/{file}"
    out_path = os.path.join(OUTDIR, f"{metabolite}.parquet")

    con.execute(f"""
        COPY (
            WITH fixed AS (
                SELECT
                    CASE
                        WHEN CAST(CHROM AS INT) = 23 THEN 'X'
                        ELSE CAST(CAST(CHROM AS INT) AS VARCHAR)
                    END                                  AS chr_str,
                    CAST(GENPOS AS BIGINT)               AS pos_int,
                    CAST(ALLELE0 AS VARCHAR)             AS a0,
                    CAST(ALLELE1 AS VARCHAR)             AS a1,
                    Effect                                AS beta,
                    StdErr                                AS se
                FROM read_parquet('{parquet_path}')
            ),
            keyed AS (
                SELECT
                    'chr' || chr_str || '_' || CAST(pos_int AS VARCHAR) || '_' || a0 || '_' || a1 AS variant,
                    beta,
                    se,
                    beta / NULLIF(se, 0) AS z
                FROM fixed
            )
            SELECT k.variant, k.beta, k.se, k.z
            FROM keyed k
            INNER JOIN lead l ON k.variant = l.variant
        )
        TO '{out_path}' (FORMAT PARQUET);
    """)
