#!/bin/bash

COSMIC_GENOME="Cosmic_GenomeScreensMutant_v100_GRCh38.tsv"
OUTPUT_DIR="/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/data"
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
awk -F"\t" '{if(length($24) != length($25) && $26 ~ /Confirmed somatic variant/) print $0}' cosmic_coding.tsv > cosmic_coding_indel.tsv

echo "Filtering for short sequences..."
awk -F"\t" '{if(length($24) < 20 && length($25) < 20) print $0}' cosmic_coding_indel.tsv > cosmic_coding_indel_short.tsv

echo "Cleaning up some unwanted columns..."
awk -F"\t" 'BEGIN {OFS="\t"} {print $15, $16, $17, $24, $25}' cosmic_coding_indel_short.tsv > cosmic_coding_indel_short_cleaned.tsv

echo "Filtering for specific chromosomes..."
awk -F"\t" 'BEGIN {OFS="\t"} {if ($1 != "X" && $1 != "Y" && $1 != "MT" && $1 != "na" && $4 != "na" && $5 != "na") print $0}' cosmic_coding_indel_short_cleaned.tsv > cosmic_coding_indel_short_filtered.tsv

echo "Removing files..."
rm cosmic_coding.tsv cosmic_coding_indel.tsv cosmic_coding_indel_short_cleaned.tsv


echo "Process completed!"