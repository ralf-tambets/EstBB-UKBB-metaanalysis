import argparse, os, sys
import pandas as pd
import numpy as np
from functools import lru_cache
from pyliftover import LiftOver
from pyfaidx import Fasta

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-file",  required=True)
    ap.add_argument("--out-dir",  required=True)
    ap.add_argument("--chain",    default="GRCh37_to_GRCh38.chain.gz")
    ap.add_argument("--fasta",    default="hg38.fa")
    return ap.parse_args()

def main():
    a = parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    out_path = os.path.join(a.out_dir, os.path.basename(a.in_file))

    if os.path.exists(out_path):
        return

    df = pd.read_parquet(a.in_file)

    df.rename({"pos":"Pos","chr":"Chr","ref":"Ref","alt":"Alt"}, axis=1, inplace=True, errors="ignore")
    need = {"Chr","Pos","Ref","Alt"}
    if not need.issubset(df.columns):
        missing = need - set(df.columns)
        raise ValueError(f"Missing columns {missing} in {a.in_file}")

    lo    = LiftOver(a.chain)
    fasta = Fasta(a.fasta)

    Chr = df["Chr"].astype("string")
    Pos = pd.to_numeric(df["Pos"], errors="coerce").astype("Int64")
    mask_valid = Chr.notna() & Pos.notna()
    chr_pos = pd.DataFrame({"Chr": Chr[mask_valid], "Pos": Pos[mask_valid].astype("int64")})

    uniq_src = chr_pos.drop_duplicates().to_records(index=False)

    @lru_cache(maxsize=None)
    def lift_one(chrom: str, pos: int):
        hit = lo.convert_coordinate(chrom, pos-1)
        if not hit:
            return (None, None)
        c38, p38 = hit[0][0], hit[0][1] + 1
        return ("chr"+str(c38), int(p38))

    lifted_map = { (c,p): lift_one(str(c), int(p)) for c,p in uniq_src }

    tgt = chr_pos.apply(lambda r: lifted_map.get((r["Chr"], int(r["Pos"])), (None, None)), axis=1)
    tgt = pd.DataFrame(tgt.tolist(), index=chr_pos.index, columns=["chr","pos"])
    df.loc[:, ["chr","pos"]] = np.nan
    df.loc[mask_valid, ["chr","pos"]] = tgt

    df = df.dropna(subset=["chr","pos"])
    df["pos"] = df["pos"].astype("Int64")

    Ref = df["Ref"].astype("string").str.upper()
    Alt = df["Alt"].astype("string").str.upper()

    chr_pos2 = df[["chr","pos"]].copy()
    chr_pos2["pos"] = chr_pos2["pos"].astype("int64")
    uniq_tgt = chr_pos2.drop_duplicates().to_records(index=False)

    def ref_base(chrom: str, pos: int) -> str:
        return str(fasta[chrom][pos-1:pos]).upper()

    ref_map = { (c,int(p)): ref_base(str(c), int(p)) for c,p in uniq_tgt }

    ref_at_site = chr_pos2.apply(lambda r: ref_map[(r["chr"], int(r["pos"]))], axis=1)
    is_ref = (ref_at_site.values == Ref.values)
    is_alt = (ref_at_site.values == Alt.values)

    df["ref"] = np.where(is_ref, Ref, np.where(is_alt, Alt, pd.NA))
    df["alt"] = np.where(is_ref, Alt, np.where(is_alt, Ref, pd.NA))
    df = df.dropna(subset=["ref","alt"])

    df.drop(["Chr","Pos","Ref","Alt"], axis=1, errors="ignore", inplace=True)

    df.to_parquet(out_path, index=False, engine="pyarrow", compression="zstd")

if __name__ == "__main__":
    main()
