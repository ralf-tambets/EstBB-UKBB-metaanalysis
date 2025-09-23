#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --job-name=panukbb_liftover
#SBATCH --partition=amd
#SBATCH --array=0-902

module load python/3.12.3
source /path/to/venv/bin/activate

IN_DIR="/path/to/PAN_UKBB"
OUT_DIR="/path/to/PAN_UKBB_lifted"
CHAIN="GRCh37_to_GRCh38.chain.gz"
FASTA="hg38.fa"

mapfile -t FILES < <(find "$IN_DIR" -maxdepth 1 -type f -name "*.parquet" | sort)
TOTAL=${#FILES[@]}

CHUNK=8   
TASK=${SLURM_ARRAY_TASK_ID:-0}
START=$(( TASK * CHUNK ))
END=$(( START + CHUNK - 1 ))

if (( START >= TOTAL )); then
  echo "Task $TASK start index $START out of range (0..$((TOTAL-1)))"
  exit 0
fi
(( END >= TOTAL )) && END=$(( TOTAL - 1 ))

for IDX in $(seq $START $END); do
  FILE=${FILES[$IDX]}
  echo "[$(date)] Task $TASK processing $FILE"
  python3 panukbb_liftover.py \
    --in-file "$FILE" \
    --out-dir "$OUT_DIR" \
    --chain "$CHAIN" \
    --fasta "$FASTA"
done
