#!/bin/bash
#SBATCH --job-name=eggnog_bacteroides             # Job name
#SBATCH --mail-type=END,FAIL                       # Adjust email notifications as needed
#SBATCH --mail-user=your.email@example.com        # Your email here
#SBATCH --ntasks=1                                # Single task
#SBATCH --cpus-per-task=16                         # Number of cores
#SBATCH --mem=60gb                                 # Memory request
#SBATCH --time=48:00:00                            # Time limit
#SBATCH --output=logs/eggnog_prokka_%j.log        # Output log

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
