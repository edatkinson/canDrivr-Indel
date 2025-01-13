#!/bin/bash

COSMIC_GENOME="Cosmic_GenomeScreensMutant_v100_GRCh38.tsv"
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

# gunzip -c $COSMIC_GENOME > cosmic_genome.tsv

# 1. Filter for coding variants
echo "Filtering for coding variants..."
cat "cosmic_genome.tsv" | awk -F"\t" '{
    if ($12 ~ /frameshift_variant|inframe_deletion|inframe_insertion/) 
        print $0 
}' > cosmic_coding.tsv

# 2. Filter for SNVs only (Single Nucleotide Variants)
echo "Filtering for Indels..."
awk -F"\t" '{if(length($24) != length($25) print $0}' cosmic_coding.tsv > cosmic_coding_snv.tsv

echo "Process completed!"