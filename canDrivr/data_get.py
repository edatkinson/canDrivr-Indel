
import pandas as pd 
import numpy as np
import requests

def extract_data(file_path):
    # Load the data
        try:
            data = pd.read_csv(file_path, sep='\t', low_memory=False)
            df = pd.DataFrame(data)

            # Extract the relevant columns
            features = df.iloc[:, [14, 15, 16, 23, 24]]
            features.columns = ['chrom', 'start', 'stop', 'wt', 'mt']

            # Filter out chromosomes X, Y, and MT
            features = features[~features['chrom'].isin(['X', 'Y', 'MT'])]

            # Reset the index for consistent processing
            features.reset_index(drop=True, inplace=True)

            features.replace(np.nan, '-', regex=True, inplace=True)
            # Group by the specified columns and count duplicates
            try:
                features['count'] = features.groupby(['start', 'stop', 'wt', 'mt'])['start'].transform('size')
            except Exception as e:
                print("Error in groupby operation:", e)
                return None

            # Drop duplicate rows, keeping the first occurrence
            features = features.drop_duplicates(subset=['start', 'stop', 'wt', 'mt'], keep='first')
            return features
        except Exception as e:
            print("Error in reading the file:", e)
            return None


if __name__ == '__main__':
    # Test the function
    result = extract_data('data/cosmic_coding_indel_short.tsv')

    if result is not None:
        print(result.head())
    else:
        print("An error occurred during processing.")



"""

Interesting conservation scoring system:

import pyBigWig

def get_conservation_scores(chrom, pos, flank_size, bw_file):
    
    Extracts conservation scores from a bigWig file for the position and flanking region.
    
    Args:
        chrom (str): Chromosome (e.g., "chr1").
        pos (int): Position of the insertion.
        flank_size (int): Number of bases to consider around the position.
        bw_file (str): Path to the bigWig file.
    
    Returns:
        dict: Conservation scores for the region.
    
    bw = pyBigWig.open(bw_file)
    start = max(0, pos - flank_size)
    end = pos + flank_size
    scores = bw.values(chrom, start, end)
    bw.close()
    return {
        "region": f"{chrom}:{start}-{end}",
        "scores": scores,
        "average_score": sum(scores) / len(scores) if scores else None
    }

# Example usage
bigwig_path = "path_to_phastcons.bigwig"
chromosome = "chr1"
position = 100
flank = 10

scores = get_conservation_scores(chromosome, position, flank, bigwig_path)
print(scores)


"""