#!/bin/bash

INPUT_FILE="accession.txt"

# Create output directory if it doesn't exist
 OUTPUT_DIR="./"
# mkdir -p "$OUTPUT_DIR"

# Read SRR IDs from input file line by line
while read -r SRR_ID || [[ -n "$SRR_ID" ]]; do

    echo "processing $SRR_ID..."

    #Download SRA file
    prefetch "$SRR_ID"
    if [ $? -ne 0 ]; then
        echo "what is this? can't download $SRR_ID?? skipping..."
        continue
    fi

    # Convert to FASTQ and split paired reads
    fasterq-dump "$SRR_ID" -O "$OUTPUT_DIR" --split-files --threads 8
    if [ $? -ne 0 ]; then
        echo "what is this? $SRR_ID can't convert? skipping..."
        continue
    else
        gzip "${SRR_ID}_1.fastq"
        gzip "${SRR_ID}_2.fastq"
        echo "yayayayay fastq down, time to delete $SRR_ID sra!"
        rm -r "$SRR_ID"
    fi

    # Delete the .sra file to save space
    # SRA_FILE_PATH=$(vdb-validate --list "$SRR_ID" 2>/dev/null | grep "$SRR_ID\.sra" | head -1)
    # if [ -z "$SRA_FILE_PATH" ]; then
    #     # If above doesn't work, try default SRA location in user cache
    #     SRA_FILE_PATH="$HOME/ncbi/public/sra/${SRR_ID}.sra"
    # # fi

    # if [ -f "$SRA_FILE_PATH" ]; then
    #     echo "Deleting $SRA_FILE_PATH to save space."
    # #     rm "$SRA_FILE_PATH"
    # # else
    #     echo "Could not find .sra file for $SRR_ID to delete."
    # fi

    echo "done with $SRR_ID"
    echo "-------------------------"

done < "$INPUT_FILE"

echo "All done! War is over for now >:)"
