
import pandas as pd 
import numpy as np
import requests

def extract_data(file_path, type):
    # Load the data

    if type == 'pathogenic':
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

        
print(extract_data('data/cosmic_coding.tsv', 'pathogenic'))

cosmic_filtered = extract_data('data/cosmic_coding.tsv', 'pathogenic')
cosmic_filtered.to_csv('data/cosmic_filtered.tsv', sep='\t', index=False)


# Test the function
# result = extract_data('cosmic_coding.tsv')

# if result is not None:
#     print(result.head())
# else:
#     print("An error occurred during processing.")




