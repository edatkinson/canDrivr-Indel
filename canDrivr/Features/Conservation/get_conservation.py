import pyBigWig
import pandas as pd
import sys
import time

# if __name__ == "__main__":
#     # variants_file = sys.argv[1]
#     # output_dir = sys.argv[2]
#     variants_file = "/home/colin/canDrivr-Indel/canDrivr/CADD_Data/final_cadd_data.bed"
#     output_dir = "/home/colin/canDrivr-Indel/canDrivr/Features/Conservation/files"

#     # Read in variant file
#     variants = pd.read_csv(variants_file, sep = "\t")
#     # Conservation files
#     print(variants.head())
#     files = [
#         "phyloP4way", "phyloP7way", "phyloP17way", "phyloP20way", "phyloP30way",
#         "phyloP100way", "phyloP470way", "phastCons4way", "phastCons7way", "phastCons17way",
#         "phastCons20way", "phastCons30way", "phastCons100way", "phastCons470way",
#         "k24.Bismap.MultiTrackMappability", "k36.Umap.MultiTrackMappability",
#         "k36.Bismap.MultiTrackMappability", "k24.Umap.MultiTrackMappability",
#         "k50.Bismap.MultiTrackMappability", "k50.Umap.MultiTrackMappability",
#         "k100.Bismap.MultiTrackMappability", "k100.Umap.MultiTrackMappability"
#     ]

#     #Combined_Cadd_data only has Position. 
#     #Therefore we must calculate the end position judging by the length of the indel.

#     for file in files:
#         variants2 = variants.copy()  
#         try:
#             if "way" in file:
#                 # Path to BigWig file
#                 bw_path = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/{file}/hg38.{file}.bw'
#             else:
#                 bw_path = f'http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/{file}.bw'

#             # Open the BigWig file
#             bw = pyBigWig.open(bw_path)
#             # Initialise an empty list to store query results
#             results = []
#             chromosome_lengths = bw.chroms()

#             for _, row in variants2.iterrows():
#                 try:
#                     # Check if the chromosome exists in the BigWig file
#                     if row['chrom'] not in chromosome_lengths:
#                         print(f"Chromosome {row['chrom']} not found in BigWig file.")
#                         results.append(float('nan'))
#                         continue

#                     # Validate start and end positions against chromosome lengths
#                     chrom_length = chromosome_lengths[row['chrom']]
#                     if row['start'] >= chrom_length or row['end'] > chrom_length:
#                         print(f"Invalid interval {row['chrom']}:{row['start']}-{row['end']} (exceeds chromosome length {chrom_length}).")
#                         results.append(float('nan'))
#                         continue

#                     # Query the BigWig file based on indel type
#                     if len(row['ref']) > len(row['alt']):
#                         # Deletion: query the region from start to end
#                         if row["start"] == row['end']:
#                             end = row['end'] + 1
#                             values = bw.values(row['chrom'], int(row["start"]), end)
#                         else:
#                             end = row['end']
#                             values = bw.values(row['chrom'], int(row['start']), end)
#                     else:
#                         # Insertion: query a flanking region of +-5 bases
#                         start = max(row['start'] - 5, 0)  # Ensure start >= 0
#                         end = row['end'] + 5
#                         if end > chrom_length:
#                             print(f"Insertion region {row['chrom']}:{start}-{end} exceeds chromosome length {chrom_length}.")
#                             results.append(float('nan'))
#                             continue
#                         values = bw.values(row['chrom'], start, end)

#                     # Remove None values and compute the average
#                     values = [v for v in values if v is not None]
#                     value = sum(values) / len(values) if values else float('nan')
#                     results.append((row['chrom'], row['start'], end, value))
#                     print(results[-1])
#                 except Exception as e:
#                     print(f"Error processing {row['chrom']}:{row['start']}-{row['end']} ({e})")
#                     results.append(float('nan'))

#             # Process the results 
#             variants2[file] = [result[3] for result in results]
#             variants2 = variants2.drop("start", axis = 1)
#             # Save to CSV
#             variants2.to_csv(f'{output_dir}hg38.{file}.bedGraph', header = None, sep = "\t", index = None)
#             bw.close()
#         except Exception as e:
#            print("cannot access file", e) # 470-way not available at the moment


import pyBigWig
import pandas as pd
import sys
from collections import defaultdict

if __name__ == "__main__":
    # Input and output file paths
    # variants_file = "/home/colin/canDrivr-Indel/canDrivr/CADD_Data/final_cadd_data.bed"
    # output_dir = "/home/colin/canDrivr-Indel/canDrivr/Features/Conservation/output/"
    variants_file = "/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/CADD_Data/final_cadd_data.bed"
    output_dir = "/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/output/"

    # Read the variants file
    variants = pd.read_csv(variants_file, sep="\t", names=['chrom','start','end', 'ref','alt'])
    '''
     
    '''
    # Conservation files
    print(variants.head())
    files = ["phyloP7way", "phyloP17way", "phyloP20way", "phyloP30way",
        "phyloP100way","phyloP470way", "phastCons7way", "phastCons17way",
        "phastCons20way", "phastCons30way", "phastCons100way", "phastCons470way",
        "k24.Bismap.MultiTrackMappability", "k36.Umap.MultiTrackMappability",
        "k36.Bismap.MultiTrackMappability", "k24.Umap.MultiTrackMappability",
        "k50.Bismap.MultiTrackMappability", "k50.Umap.MultiTrackMappability",
        "k100.Bismap.MultiTrackMappability", "k100.Umap.MultiTrackMappability"
    ]

    for file in files:
        variants2 = variants.copy()
        try:
            # Set the BigWig file path
            if "way" in file:
                bw_path = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/{file}/hg38.{file}.bw'
            else:
                bw_path = f'http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/{file}.bw'

            # Open the BigWig file
            bw = pyBigWig.open(bw_path)
            chromosome_lengths = bw.chroms()

            # Group variants by chromosome for batch processing
            chrom_groups = defaultdict(list)
            for _, row in variants2.iterrows():
                chrom_groups[row['chrom']].append(row)

            results = []

            # Process each chromosome group
            for chrom, rows in chrom_groups.items():
                if chrom not in chromosome_lengths:
                    print(f"Chromosome {chrom} not found in BigWig file.")
                    results.extend([float('nan')] * len(rows))
                    continue

                chrom_length = chromosome_lengths[chrom]
                intervals = []

                # Collect intervals for batch query
                for row in rows:
                    if len(row['ref']) > len(row['alt']):  # Deletion
                        start, end = row['start'], row['end']
                        if start == end:
                            end += 1  # Ensure end > start
                    else:  # Insertion
                        start = max(row['start'] - 5, 0)
                        end = min(row['end'] + 5, chrom_length)

                    intervals.append((start, end))

                # Batch query intervals for the chromosome
                try:
                    batch_values = [bw.values(chrom, start, end) for start, end in intervals]

                    # Compute average for each interval
                    for row, values in zip(rows, batch_values):
                        values = [v for v in values if v is not None]
                        avg_value = sum(values) / len(values) if values else float('nan')
                        results.append(avg_value)
                        print(f"{chrom}:{row['start']}-{row['end']} -> {avg_value}")
                except Exception as e:
                    print(f"Error querying chromosome {chrom}: {e}")
                    results.extend([float('nan')] * len(rows))

            # Add results to the DataFrame
            variants2[file] = results
            variants2 = variants2.drop("end", axis=1)
            variants2 = variants2.rename(columns={"start":"pos", "alt":"alt_allele","ref":"ref_allele"})

            # Save the updated DataFrame to a file
            variants2.to_csv(f'{output_dir}hg38.{file}.bedGraph', header=True, sep="\t", index=None)
            bw.close()
        except Exception as e:
            print(f"Cannot access file {file}: {e}")
