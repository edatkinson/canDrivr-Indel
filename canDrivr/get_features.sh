#!/bin/bash

# Features to automate:
# DNA_shape, dinucleotide properties, gc CpG, conservation


# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <variantDir> <variantFile> <outputDir>"
    exit 1
fi

# Assign arguments
variantDir=$1
variantFile=$2
outputDir=$3

# Ensure output directory exists
mkdir -p "$outputDir"


# Navigate to variant directory
cd "$variantDir" || { echo "Error: Cannot change to directory $variantDir"; exit 1; }

# Construct new file name by appending "_no_driver_stat"
# variantFileName=$(basename "$variantFile")
# newFile="${variantDir}${variantFileName}_no_driver_stat"

# Extract filename and extension
filename=$(basename "$variantFile")
extension="${filename##*.}"
basename="${filename%.*}"


newFile="${basename}_no_driver_stat.${extension}"
# Remove driver status (assumes last column is driver status)
awk '{print $1"\t"$2"\t"$3"\t"$4}' "$variantFile" > "$newFile"

# echo "Processed file saved to: $newFile"

cd /Users/edatkinson/Repos/canDrivr-Indel/canDrivr/CADD_Data/

python3 process_cadd.py $variantFile "process_data.bed" $outputDir
