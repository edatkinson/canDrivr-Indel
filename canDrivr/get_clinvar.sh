#!/bin/bash

CLINVAR_GENOME="variant_summary.txt"
OUTPUT_DIR="/Users/edatkinson/Repos/canDrivr-Indel/canDrivr"
# Ensure we are in the correct directory
if [ -d "$OUTPUT_DIR" ]; then
    cd $OUTPUT_DIR
else
    echo "Directory $OUTPUT_DIR does not exist."
    exit 1
fi
ls

# Ensure we are in the correct directory
cd $OUTPUT_DIR


# 1. Filter for coding variants
echo "Filtering for coding variants..."
cat "variant_summary.txt" | awk -F"\t" '{
    if ($2 ~ /(Indel|Insertion|Deletion)/ && $7 ~ /Benign/ && $17 ~ /GRCh38/)
        print $0 
}' > clinvar_coding.tsv


# Extract specific columns and filter out unwanted chromosomes
echo "Processing columns and filtering chromosomes..."
awk -F"\t" 'BEGIN {OFS="\t"} {print $19, $20, $21, $33, $34}' clinvar_coding.tsv | \
    awk -F"\t" 'BEGIN {OFS="\t"} {if ($1 != "X" && $1 != "Y" && $1 != "MT" && $1 != "na" && $4 != "na" && $5 != "na") print $0}' > clinvar_filtered.tsv

echo "Filtering for short sequences..."
awk -F"\t" 'BEGIN {OFS="\t"} {if (length($4) < 20 && length($5) < 20) print $0}' clinvar_filtered.tsv > clinvar_short_sequences.tsv

echo "Cleaning up redundant files..."
rm clinvar_filtered.tsv

echo "Process completed!"
