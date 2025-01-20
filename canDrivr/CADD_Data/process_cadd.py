
from heapq import merge
from itertools import count
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
    # Special cases:
        # Deletion of length 1 leads to start = end.
        # There are some cases where there are both insertions and deletions leading to start > end which is impossible

    data['end'] = data['pos'] + data['lengths'] - 1  # Default case for deletions 

    is_insertion = data['length_alt'] >= data['length_ref']  # Identify insertions
    data.loc[is_insertion, 'end'] = data.loc[is_insertion, 'pos'] + 1  # Adjust for insertions

    def adjust_positions(row):
        if row['pos'] == row['end']:
            row['end'] = row['pos'] + 1
        return row

    # Apply the function to each row
    data = data.apply(adjust_positions, axis=1)

    # Drop intermediate columns and rename as needed
    data.drop(columns=['length_ref', 'length_alt', 'lengths'], inplace=True)
    data.rename(columns={'pos': 'start'}, inplace=True)

    # Reorder columns to the desired format 
    data = data[['chrom', 'start', 'end', 'ref', 'alt', 'driver_stat']]

    return data


data = calc_length_indel('combined_cadd_data.tsv')


def select_equal_amounts(data):
    """
    Randomly select an equal number of simulated variants for each chromosome
    based on the counts of human-derived variants.
    
    Args:
        data (pd.DataFrame): Input DataFrame with columns 'driver_stat' and 'chrom'.
        
    Returns:
        pd.DataFrame: A DataFrame of simulated variants sampled to match human-derived counts.
    """
    # Filter for human-derived and simulated indels
    human_derived = data[data['driver_stat'] == 1]
    simulated_indels = data[data['driver_stat'] == 0]

    # Count the number of human-derived indels per chromosome
    count_chrom = human_derived.groupby('chrom').size()

    # Initialize an empty list to store sampled simulated indels
    sampled_simulated = []

    # Iterate through each chromosome and select the same number of simulated indels
    for chrom, count in count_chrom.items():
        # Filter simulated indels for the current chromosome
        sim_chrom = simulated_indels[simulated_indels['chrom'] == chrom]
        
        # Randomly sample 'count' indels, ensuring we don't exceed the available number
        if count <= len(sim_chrom):
            sampled_simulated.append(sim_chrom.sample(n=count, random_state=42))
        else:
            # If there are not enough simulated indels, sample all available
            sampled_simulated.append(sim_chrom)

    # Concatenate all sampled simulated indels into a single DataFrame
    result = pd.concat(sampled_simulated, ignore_index=True)
    result = pd.concat([human_derived, result], ignore_index=True)
    return result.iloc[:, [0,1,2,3,4]]

   
select_equal_amounts(data).to_csv('final_cadd_data.bed', sep='\t', header=None, index=None)


