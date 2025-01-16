#!/bin/bash


CADD_GENOME="simulation_InDels.tsv"
OUTPUT_DIR="/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/CADD_Data"
# Ensure we are in the correct directory
if [ -d "$OUTPUT_DIR" ]; then
    cd $OUTPUT_DIR
else
    echo "Directory $OUTPUT_DIR does not exist."
    exit 1
fi

# Extract specific columns, filter out unwanted chromosomes, and filter for short sequences
echo "Processing columns, filtering chromosomes, and filtering for short sequences..."
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' $CADD_GENOME | \
    awk -F"\t" 'BEGIN {OFS="\t"} {if ($7 ~ /CodingTranscript/ && $7 !~ /NonCodingTranscript/) {print $0; print $0 > "cadd_training_data.tsv"}}'

echo "Process completed!"
