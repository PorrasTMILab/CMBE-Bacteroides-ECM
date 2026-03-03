#!/bin/bash
#SBATCH --job-name=merops_phmmer_bacteroides          # Job name
#SBATCH --mail-type=END,FAIL                           # Adjust email notifications as needed
#SBATCH --mail-user=your.email@example.com             # Your email here
#SBATCH --ntasks=1                                    # Single CPU task
#SBATCH --cpus-per-task=2                             # CPU cores per task
#SBATCH --mem=10gb                                    # Memory per task
#SBATCH --time=12:00:00                               # Time limit
#SBATCH --output=logs/phmmer_prokka_%A_%a.log         # Output log (array job)
#SBATCH --array=0-63                                  # Adjust to genome array size

# === INPUT FILES ===
INPUT_DIR="../results/prokka_runs/combined_faa"       # Folder with *.faa files from Prokka
OUTPUT_DIR="../results/merops_phmmer_out"              # Output folder for MEROPS phmmer results
MEROPS_DB="../databases/merops_pepunit.fa"             # MEROPS database fasta
MEROPS_META="../databases/merops_meta_processed.csv"   # MEROPS metadata file

# Get .faa files from input directory sorted
FILES=($(find "$INPUT_DIR" -name "*.faa" | sort))
QUERY="${FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$QUERY" .faa)

TBL_FILE="$OUTPUT_DIR/${BASENAME}_phmmer.tbl"
CSV_FILE="$OUTPUT_DIR/${BASENAME}_phmmer_annotated.csv"
mkdir -p "$OUTPUT_DIR"

# Skip if output already exists
if [[ -f "$TBL_FILE" ]]; then
  echo "⏩ Skipping $BASENAME — already annotated."
  exit 0
fi

# === Run phmmer ===
phmmer --cpu 2 --tblout "$TBL_FILE" "$QUERY" "$MEROPS_DB"

# === Annotate top hits ===
python3 - << END_SCRIPT
import pandas as pd
meta = pd.read_csv("$MEROPS_META", sep="\t")
meta = meta.set_index("MEROPS_ID").to_dict("index")
hits = []
with open("$TBL_FILE") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split()
        if len(fields) < 13:
            continue
        target, query, evalue, score = fields[0], fields[2], float(fields[4]), float(fields[5])
        meta_info = meta.get(target, {"Description": "NA", "Peptidase_Family": "NA", "Peptidase_Class": "NA"})
        hits.append({"Query_ID": query, "MEROPS_ID": target, "E-value": evalue, "Score": score, **meta_info, "Genome": "$BASENAME"})
if hits:
    df = pd.DataFrame(hits)
    df['rank'] = df.groupby("Query_ID")['E-value'].rank(method='first')
    df[df['rank'] == 1].drop(columns=['rank']).to_csv("$CSV_FILE", index=False)
END_SCRIPT

echo "✅ Done: $BASENAME"
