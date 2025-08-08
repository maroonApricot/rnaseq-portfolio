#!/bin/bash

usage() {
    echo "Usage: $0 -i <salmon_index> -f <fastq_dir> -o <output_directory>"
    echo "  -i: Path to Salmon index directory"
    echo "  -f: Directory containing single-end FASTQ files"
    echo "  -o: Path to place output files"
    exit 1
}

THREADS=8

# Parse command line arguments
while getopts ":i:f:o:" opt; do
  case $opt in
    i) SALMON_INDEX="$OPTARG" ;;
    f) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    *) usage ;;
  esac
done

# Check required parameters
if [ -z "$SALMON_INDEX" ] || [ -z "$FASTQ_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

for R in "$FASTQ_DIR"/*.fastq.gz; do
    SAMPLE=$(basename "$R" ".fastq.gz")
    # R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"

    # if [[ ! -f "$R2" ]]; then
    #     echo "Warning: Missing R2 for $SAMPLE, skipping."
    #     continue
    # fi

    echo "Running Salmon quant for sample: $SAMPLE"
    # i run without gc bias because it's experimental and provides little benefit
    salmon quant -i "$SALMON_INDEX" \
                -l A \
                -r "$R" \
                --validateMappings \
                --seqBias \
                -p "$THREADS" \
                -o "${OUTPUT_DIR}/${SAMPLE}"
done

echo "All Salmon quantification completed."
