
from heapq import merge
import pandas as pd

def process_cadd(cadd_file):
    cadd = pd.read_csv(cadd_file, sep='\t', header=None)
    cadd = cadd.iloc[:,[0,1,2,3]]
    cadd.columns = ['chrom', 'pos', 'ref', 'alt']
    cadd['chrom'] = cadd['chrom'].apply(lambda x: "chr"+str(x))
    if 'sim' in cadd_file:
        cadd['driver_stat'] = [0 for i in range(len(cadd))]
    elif 'human' in cadd_file:
        cadd['driver_stat'] = [1 for i in range(len(cadd))]
    else:
        print("Please provide data which is human or simulated based. (Has the word 'Human' of 'simulated' in the file name)")
        return None

    return cadd

def merge_cadd(cadd_file, cadd_file2):
    cadd = process_cadd(cadd_file)
    cadd2 = process_cadd(cadd_file2)

    cadd = pd.concat([cadd, cadd2], axis=0, ignore_index=True)

    return pd.DataFrame(cadd)


def calc_length_indel(combined_file):
    # Read the input file
    data = pd.read_csv(combined_file, sep='\t')

    # Calculate lengths for ref and alt strings
    data['length_ref'] = data['ref'].str.len()
    data['length_alt'] = data['alt'].str.len()

    # Calculate the absolute difference in lengths
    data['lengths'] = (data['length_ref'] - data['length_alt']).abs()

    # Determine 'end' points based on insertion or deletion
    data['end'] = data['pos'] + data['lengths'] - 1  # Default case for deletions
    is_insertion = data['length_alt'] > data['length_ref']  # Identify insertions
    data.loc[is_insertion, 'end'] = data.loc[is_insertion, 'pos'] + 1  # Adjust for insertions

    # Drop intermediate columns and rename as needed
    data.drop(columns=['length_ref', 'length_alt', 'lengths'], inplace=True)
    data.rename(columns={'pos': 'start'}, inplace=True)

    # Reorder columns to the desired format
    data = data[['chrom', 'start', 'end', 'ref', 'alt', 'driver_stat']]

    return data


print(calc_length_indel('combined_cadd_data.tsv'))




