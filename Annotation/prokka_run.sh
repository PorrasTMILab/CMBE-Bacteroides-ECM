#!/bin/bash
#SBATCH --job-name=prokka_bacteroides                   # job name
#SBATCH --mail-type=END,FAIL                             # adjust email notifications as needed
#SBATCH --mail-user=your.email@example.com               # your email here
#SBATCH --ntasks=1                                       # single CPU task
#SBATCH --cpus-per-task=8                                # Number of CPU cores
#SBATCH --mem=32gb                                       # memory
#SBATCH --time=24:00:00                                  # time limit
#SBATCH --output=logs/run_prokka_%j.log                  # output log file

# Deactivate conda and Java interference
conda deactivate 2>/dev/null || true
unset JAVA_HOME
hash -r

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
