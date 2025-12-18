from __future__ import annotations

import gzip
import io
from pathlib import Path
from urllib.parse import urlparse

import boto3
import pandas as pd
from botocore import UNSIGNED
from botocore.client import Config

MANIFEST = "phenotype_manifest.tsv" # the panukbb manifest accessible at https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?gid=1450719288#gid=1450719288     
DEST     = Path(
    "/path/to/PAN_UKBB"
)
DEST.mkdir(parents=True, exist_ok=True)

s3 = boto3.client(
    "s3",
    config=Config(
        signature_version=UNSIGNED,
        read_timeout=3000,                   
        connect_timeout=30,             
        retries={"max_attempts": 10, "mode": "standard"},  
    ),
)

manifest_df = pd.read_csv(MANIFEST, sep="\t", usecols=["aws_path"])

for uri in manifest_df["aws_path"]:
    p = urlparse(uri)                   
    bucket, key = p.netloc, p.path.lstrip("/")
    out_file = DEST / (Path(key).stem.replace(".tsv", "") + ".parquet")

    if out_file.exists():
        print(f"✓ {out_file.name} already done, skipping")
        continue

    obj = s3.get_object(Bucket=bucket, Key=key)
    with gzip.GzipFile(fileobj=obj["Body"]) as gz:
        df = pd.read_csv(
            io.TextIOWrapper(gz),
            sep="\t",
            low_memory=False, 
        )

    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].astype("string")

    df.to_parquet(out_file, index=False)
