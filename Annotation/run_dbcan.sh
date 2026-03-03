#!/bin/bash
#SBATCH --job-name=dbcan_bacteroides
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24gb
#SBATCH --time=48:00:00
#SBATCH --output=logs/run_dbcan_all_%j.log

# === Load module if needed ===
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
