#!/bin/bash
#SBATCH --job-name=run_prokka.slurm                         #job name
#SBATCH --mail-type=ALL
#SBATCH --mail-user=                                        #where to send email    
#SBATCH --ntasks=1                                          #run on a single CPU
#SBATCH --cpus-per-task=8                                   #Number of cpu cores to use (adjust as needed)
#SBATCH --mem=32gb                                          #job memory request
#SBATCH --time=24:00:00                                     #time limit hrs:min:sec
#SBATCH --output=run_prokka_%j.log                          #standard output and error log


# Deactivate conda and Java interference
conda deactivate 2>/dev/null || true
unset JAVA_HOME
hash -r

# Load UFRC Prokka module
module purge
module load prokka/1.14.6

# Define input and output directories
INPUT_DIR="/home/username/project/data/"
OUTPUT_BASE="/home/username/project/data/"

# Make output directory if needed
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
