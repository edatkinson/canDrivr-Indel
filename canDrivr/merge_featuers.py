import pandas as pd
import os
import numpy as np

"""
Sort the final_cadd_data_stat in the same way as the other data, this ensures that we don't need to do any matching and the merging process is swift.
- Done called /Users/edatkinson/Repos/canDrivr-Indel/canDrivr/CADD_Data/final_cadd_data_stat.sorted.bed

Need to locate all the files which contain features for my indels so far: done using locate_files()

Merge them together.

Files:
All files are within a folder called 'output' in there respective feature folder
"""



def locate_files(dir):
    output_paths = []
    for outputs in os.listdir(path=dir):
        if outputs == 'indtest':
            files = os.listdir(f'{dir}/{outputs}')
            for file in files:
                if file != '':
                    output_paths.append(f'{dir}/{outputs}/{file}')
        # elif outputs == 'pathogenic_out':
        #     files = os.listdir(f'{dir}/{outputs}')
        #     for file in files:
        #         if file != '':
        #             output_paths.append(f'{dir}/{outputs}/{file}')
        else:
            continue
    
    return output_paths


    
def merge_cons(cos_data, ben_data):
    file_paths = os.listdir('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/pathogenic_out')
    cos_data_formatted = cos_data[['chrom', 'start', 'ref_allele','alt_allele', 'driver_stat']]
    ben_data_formatted = ben_data[['chrom', 'start', 'ref_allele','alt_allele', 'driver_stat']]
    ben_data_formatted.rename(columns={'start':'pos'}, inplace=True)
    cos_data_formatted.rename(columns={'start':'pos'}, inplace=True)
    # print(file_paths)
    for file in file_paths:
        path  = os.path.join('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/pathogenic_out', file)
        data = pd.read_csv(path, sep='\t')
        cos_data_formatted = pd.merge(cos_data_formatted, data, on=['chrom','pos', 'ref_allele', 'alt_allele', 'driver_stat'])

    file_paths = os.listdir('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/output')
    # print(file_paths)
    for file in file_paths:
        path = os.path.join('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/output', file)
        data_ben = pd.read_csv(path, sep='\t')
        ben_data_formatted = pd.merge(ben_data_formatted, data_ben, on=['chrom','pos', 'ref_allele', 'alt_allele', 'driver_stat'])

    combined_data = pd.concat([cos_data_formatted, ben_data_formatted])
    combined_data = combined_data[combined_data['chrom'] != 'chrX']
    combined_data.to_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/merged_features/merged_cons.tsv', sep='\t', index=None)





def merge_files(file_paths, dataset_file):
    # Baseline data
    data = pd.read_csv(dataset_file, sep=',', names=['chrom', 'start','end','ref_allele', 'alt_allele', 'driver_stat'])
    data.drop(columns='end', inplace=True)


    # loop over file_paths and merge them on data

    for file in file_paths:
        splitted = file.split('/')
        features = pd.read_csv(file, sep='\t')
        # elif 'ProteinStructure' in splitted:
        #     features = pd.read_csv(file, sep='\t')
        # else:
        #     features = pd.read_csv(file, sep='\t')

        data = pd.merge(data, features, on=['chrom','pos', 'ref_allele', 'alt_allele'])
        # print(data.head())
    return data

def merge_data():
    list_of_features = ['/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/DNA_Shape/output','/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/FG5_gc_CpG/output','/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/FG3_dinucleotide_properties/output']

    for dir in list_of_features:
        merged = pd.DataFrame()
        files = os.listdir(dir)

        for file in files:
            filepath = os.path.join(dir, file)  # Use os.path.join for robust path construction
            try: # Try to read the file, handle potential errors
                data = pd.read_csv(filepath, sep='\t')
            except pd.errors.EmptyDataError: # Handle empty files
                print(f"Warning: File {filepath} is empty. Skipping.")
                continue
            except FileNotFoundError: # Handle if file not found (unlikely but good practice)
                print(f"Warning: File {filepath} not found. Skipping.")
                continue
            except Exception as e: # Catch any other exceptions during file reading
                print(f"Error reading file {filepath}: {e}. Skipping.")
                continue


            if 'driver_stat' not in data.columns:  # Check if column is absent
                if file.split('_')[-1] == 'cos.txt':
                    data['driver_stat'] = 1  # More efficient than insert for new columns
                    print(f"Added driver_stat=1 to {file}:") # Print the file where the column was added
                    print(data.head())
                else:
                    data['driver_stat'] = 0  # More efficient than insert for new columns
                    print(f"Added driver_stat=0 to {file}:") # Print the file where the column was added
                    print(data.head())
            else: # If driver_stat exists, just print the head
                print(f"driver_stat column found in {file}:")
                print(data.head())

            merged = pd.concat([merged, data], ignore_index=True) # ignore_index for clean concatenation

        # Save the merged dataframe for each directory if needed
        if not merged.empty: # Check if merged is not empty before saving
            output_filename = f"{dir}/merged_data.csv" # Or appropriate extension
            merged.to_csv(output_filename, sep='\t', index=False)
            print(f"Merged data saved to {output_filename}")
        else:
            print(f"No data to merge in {dir}")
    

def merge_all(dir, cons_data):
    list_dir = os.listdir(dir)
    for file in list_dir:
        if file != 'merged_cons.tsv':
            path = os.path.join(dir, file)
            data = pd.read_csv(path, sep='\t')
            cons_data = pd.merge(cons_data, data, on=['chrom', 'pos', 'alt_allele', 'ref_allele','driver_stat'], how='right')

    cons_data = cons_data.fillna(0)
    cons_data.to_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Anno_data.csv', index=None)

def merge_alpha(alpha, all):
    all_data = pd.read_csv(all, sep=',')
    print(all_data.head())
    alpha_data = pd.read_csv(alpha, sep='\t')
    print(alpha_data.head())

    merged = pd.merge(all_data, alpha_data, on=['chrom', 'pos','ref_allele', 'alt_allele', 'driver_stat'], how='outer', indicator=True)

    # Add a column indicating whether the row was merged or not
    merged['merged'] = merged['_merge'].apply(lambda x: True if x == 'both' else False)

    # Drop the _merge column as it is no longer needed
    merged = merged.drop(columns='_merge')

    merged.to_csv('merged_and_unmerged_data.csv', sep=',', index=None)

if __name__ == '__main__':

    file_paths = locate_files('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation')
    # merge_files(file_paths, '/Users/edatkinson/Repos/canDrivr-Indel/independentTesting/independent_test_set_combined.csv')
    # print(file_paths)
    # cos_data = pd.DataFrame(pd.read_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/gnomad/modified_2filtered_cosmic.tsv', sep='\t'))
    # ben_data = pd.DataFrame(pd.read_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/gnomad/benign_stat.tsv', sep='\t'))

    # data = merge_files(file_paths, dataset_file)

    # data.to_csv('Annotated_data.csv', sep=',', index=False)
    # cons_data = pd.read_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/merged_features/merged_cons.tsv', sep='\t')
    # merge_all('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/merged_features', cons_data)
    # merge_alpha('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/ProteinStructure/atomic_average_data.tsv','/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Anno_data.csv')

