#==========================================
#Prokka Annotation
#==========================================


# Load UFRC Prokka module
module purge
module load prokka/1.14.6

# Define input and output directories
INPUT_DIR="../genomes"                  # Folder with input .fna files
OUTPUT_BASE="../results/prokka_runs"   # Output directory for Prokka results

# Make sure output directory exists
mkdir -p "$OUTPUT_BASE"

# Loop through each .fna file and run Prokka
for FILE in "$INPUT_DIR"/*.fna; do
    BASENAME=$(basename "$FILE" .fna)
    OUTDIR="$OUTPUT_BASE/${BASENAME}"
    echo "Running Prokka on $BASENAME"

    mkdir -p "$OUTDIR"

    prokka \
        --outdir "$OUTDIR" \
        --prefix "$BASENAME" \
        --cpus 8 \
        --force \
        "$FILE"
done


#==========================================
#eggNOG
#==========================================

# Load eggnog-mapper
module load eggnog-mapper/2.1.12

# Define input and output directories
IN_DIR="../results/prokka_runs"                   # Folder with Prokka *.faa files
OUT_DIR="../results/eggnog_annotations"           # Output folder for eggNOG results

# Make sure output directory exists
mkdir -p "$OUT_DIR"

# Loop through each Prokka-generated .faa file
for faa in "$IN_DIR"/*/*.faa; do
    BASENAME=$(basename "$faa" .faa)
    echo "🧬 Processing $BASENAME..."
    emapper.py -i "$faa" \
        -o "$BASENAME" \
        --output_dir "$OUT_DIR" \
        --cpu 16 \
        --itype proteins
done


#==========================================
#run_dbcan
#==========================================

# === Load module ===
module spider run_dbcan

# === Define directories ===
GENOME_DIR="../genomes"              # Folder containing genome .fna files
OUT_BASE="../results/dbcan_runs"     # Output folder for dbCAN results
DB_DIR="../databases/dbcan_db"       # dbCAN database directory

# === Loop over genomes ===
for fna_file in "$GENOME_DIR"/*.fna; do
    BASENAME=$(basename "$fna_file" .fna)
    OUTDIR="${OUT_BASE}/${BASENAME}_dbCAN"
    echo "🔍 Processing $BASENAME"
    mkdir -p "$OUTDIR"
    run_dbcan easy_substrate \
        --input_raw_data "$fna_file" \
        --mode prok \
        --output_dir "$OUTDIR" \
        --db_dir "$DB_DIR" \
        --gff_type prodigal
done



#==========================================
#MEROPS
#==========================================

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
