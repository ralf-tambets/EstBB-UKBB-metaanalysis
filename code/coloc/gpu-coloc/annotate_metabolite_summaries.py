import os
import pandas as pd
import requests, pandas as pd, math, threading
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from tqdm import tqdm

_SESSION_LOCK = threading.Lock()
_SESSION = None

def _get_session():
    global _SESSION
    if _SESSION is None:
        with _SESSION_LOCK:
            if _SESSION is None:
                s = requests.Session()
                retries = Retry(
                    total=6, backoff_factor=0.5,
                    status_forcelist=[429, 500, 502, 503, 504],
                    allowed_methods={"GET"}
                )
                adapter = HTTPAdapter(pool_connections=64, pool_maxsize=64, max_retries=retries)
                s.mount("https://", adapter)
                s.headers.update({"Accept": "application/json", "User-Agent": "nearest-genes/1.0"})
                _SESSION = s
    return _SESSION

@lru_cache(maxsize=None)
def _get_chrom_lengths():
    url = "https://rest.ensembl.org/info/assembly/homo_sapiens"
    r = _get_session().get(url, timeout=20, headers={"Accept": "application/json"})
    r.raise_for_status()
    data = r.json()
    return {c["name"]: int(c["length"]) for c in data["top_level_region"] if c["coord_system"] == "chromosome"}

def _clamp_region(chrom, start, end):
    chrom_lengths = _get_chrom_lengths()
    if chrom not in chrom_lengths:
        raise ValueError(f"Unknown chromosome: {chrom}")
    length = chrom_lengths[chrom]
    return max(1, start), min(length, end)

def _fetch_genes(chrom: str, pos: int, win: int = 1_500_000):
    start, end = _clamp_region(chrom, pos - win, pos + win)
    url = f"https://rest.ensembl.org/overlap/region/homo_sapiens/{chrom}:{start}-{end}?feature=gene"
    r = _get_session().get(url, timeout=20)
    r.raise_for_status()
    return r.json()

@lru_cache(maxsize=None)
def get_nearest_genes(chrom: str, pos: int, win: int = 1_500_000):
    genes = _fetch_genes(str(chrom), int(pos), int(win))
    if not genes:
        return []

    rows = []
    for g in genes:
        strand = g.get("strand", 1)
        start, end = int(g["start"]), int(g["end"])
        tss = start if strand == 1 else end
        tes = end if strand == 1 else start
        d_tss = (pos - tss) * strand
        d_tes = (pos - tes) * strand
        if g.get("external_name") is not None and g.get("biotype")=="protein_coding":
            rows.append({
                "gene_id": g["id"],
                "symbol": g.get("external_name"),
                "d_TSS": int(d_tss),
                "d_TES": int(d_tes),
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("d_TSS", key=lambda s: s.abs()).reset_index(drop=True)
    return df.head(3).to_dict(orient="records")

def annotate_nearest_genes(df: pd.DataFrame, chr_col="chr", pos_col="pos", win=1_500_000, workers=10):
    def _norm_chr(c):
        c = str(c)
        return c[3:] if c.lower().startswith("chr") else c

    tasks = [( _norm_chr(df.at[i, chr_col]), int(df.at[i, pos_col]) ) for i in range(len(df))]

    out = [None] * len(df)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(get_nearest_genes, c, p, win): idx for idx, (c, p) in enumerate(tasks)}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Nearest genes"):
            idx = futures[fut]
            try:
                out[idx] = fut.result()
            except Exception as e:
                print(f"No genes {e}")
                out[idx] = []
                
    df = df.copy()
    df["nearest_genes"] = out
    return df

pop_dir_tag = {
    "EstBB":"Est","UKBB_AFR":"AFR","UKBB_AMR":"AMR","UKBB_CSA":"CSA",
    "UKBB_EAS":"EAS","UKBB_EUR":"EUR","UKBB_MID":"MID", "meta_EUR":"meta_EUR"
}

os.makedirs("annotated", exist_ok=True)

for prefix in pop_dir_tag.keys():
    print(prefix)
    df = pd.read_csv(f"meta_summaries/{prefix}_summary.tsv", sep="\t")

    df["chr"] = df["lead_variant"].str.split("_").str[0].str.replace("chr", "")
    df["pos"] = df["lead_variant"].str.split("_").str[1]

    df["trait"] = df["signal"].str.replace(r"_chr.*", "", regex=True).str.replace(f"{prefix}_", "", regex=True)

    lead = pd.read_csv(f'leadid/{pop_dir_tag[prefix]}_loci_merged_full.tsv', sep="\t")

    lead.dropna(subset=["CHR"], inplace=True)

    lead["maf"] = lead["AF_ALL1"].apply(lambda x: min(x, 1 - x))

    lead.rename(columns={"LOG10P": "log10p"}, inplace=True)

    lead["CHR"] = lead["CHR"].apply(lambda x: "X" if int(x) == 23 else int(x))

    lead["lead_variant"] = "chr" + lead["CHR"].astype(str) + "_" + lead["SNP"].str.split(":").str[1].astype(str)

    lead = lead[["metabolite", "lead_variant", "maf", "log10p"]]

    df["metabolite"] = df["signal"].str.replace(r"_chr.*", "", regex=True).str.replace(f"{prefix}_", "", regex=True)

    df = df.merge(lead, on=("metabolite", "lead_variant"), how="inner")

    annotate_nearest_genes(df, workers=10).to_csv(f"annotated_meta_summaries/{prefix}_summary_annotated.tsv", sep="\t", index=False)