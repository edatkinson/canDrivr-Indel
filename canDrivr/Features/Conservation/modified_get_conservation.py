
import pandas as pd 
import random
import numpy as np
import requests
import json
import pyBigWig
import time


def get_conservation_scores(variants_file, output_dir):
    '''
    
    - Obtains conservtion scores for the variants in the variants file.
    - The conservation scores are obtained from the UCSC genome browser.
    - Scores for Insertions and Deletions are calculated by taking the average of the scores for each position within the range of INDELs position.
    - 
    '''
    
    # Read in variant file
    variants = pd.read_csv(variants_file, sep="\t", header=None, names=["chrom", "start", "end", "ref_allele", "alt_allele", 'driver'])
    variants.drop(columns=['driver'], inplace=True)
    # variants = variants.drop_duplicates()
    #remove rows with duplicate start and stop positions:
    variants = variants.drop_duplicates(subset=['chrom', 'start', 'end'])
    # Remove variants with alt or ref allele length > 20
    variants = variants[variants["alt_allele"].apply(lambda x: pd.isna(x) or len(x) <= 20 )]
    variants = variants[variants["ref_allele"].apply(lambda x: pd.isna(x) or len(x) <= 20)]
    # done: 

    # variants = variants.sample(n=10, random_state=42)
    # Conservation file

    files = [
        "phastCons17way",
        "phastCons20way", "phastCons30way", "phastCons100way", "phastCons470way",
        "k24.Bismap.MultiTrackMappability", "k36.Umap.MultiTrackMappability",
        "k36.Bismap.MultiTrackMappability", "k24.Umap.MultiTrackMappability",
        "k100.Bismap.MultiTrackMappability", "k100.Umap.MultiTrackMappability",
        "k24.Umap.MultiTrackMappability","phyloP4way", "phyloP7way", "phyloP17way", 
        "phyloP20way", "phyloP30way", "phyloP100way", "phyloP470way", "phastCons4way", "phastCons7way"
    ]


    for file in files:
        variants2 = variants.copy()
        try:
            if "way" in file:
                # Path to BigWig file
                bw_path = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/{file}/hg38.{file}.bw'
            else:
                bw_path = f'http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/{file}.bw'

            # Open the BigWig file
            bw = pyBigWig.open(bw_path)

            # Initialise an empty list to store query results
            results = []

            #Insertions are always 1bp long, deletions can be longer, 0 length or 1bp long

            # Loop through the DataFrame rows
            for index, row in variants2.iterrows():
                #Insertions scores are calculated by taking the values of the adjacent positions and averaging them
                #Deletions scores are calculated by taking the values of the positions and averaging them
                if row['start'] != row['end']: # both insertions and deletions
                    values = bw.values(row['chrom'], row['start'], row['end'])
                    value = sum(values) / len(values)  # Average value over the range

                else: # if the start and end positions are the same, it is a deletion of length 1
                    value = bw.values(row['chrom'], row['start'], row['end']+1)[0]
                results.append((row['chrom'], row['start'], row['end'], value))
                print(results[-1])
                
            # Process the results
            variants2[file] = [result[3] for result in results]
            variants2 = variants2.drop("start", axis=1)
            # Save to CSV
            variants2.to_csv(f'{output_dir}hg38.{file}.bedGraph', header=None, sep="\t", index=None)
            print(f'Finished processing {file}')  
            bw.close()
            time.sleep(5)
        except Exception as e:
            print(f"Cannot access file {file}: {e}")  # Print the specific error message
            break

#For the driver variants

# variants_file = "filtered_data/cosmic_filtered.bed"
# output_dir = "conservation_scores/conservation_files/"

# get_conservation_scores(variants_file, output_dir) # This will take a while to run

#For the neutral variants

variants_file = "clinvar/filtered_clinvar_driver_stat.bed"
output_dir = "conservation_scores/conservation_files/"

get_conservation_scores(variants_file, output_dir) # This will take a while to run
