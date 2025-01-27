import pandas as pd
import os

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
    for folder in os.listdir(path=dir):
        # Logic to access each folder in the path's directory
        out_path = f'{dir}/{folder}'
        if os.path.isdir(out_path):
            for outputs in os.listdir(path=out_path):
                if outputs == 'output':
                    files = os.listdir(f'{out_path}/{outputs}')
                    for file in files:
                        if file != '':
                            output_paths.append(f'{out_path}/{outputs}/{file}')
                else:
                    continue
        else:
            continue
    
    return output_paths



def merge_files(file_paths, dataset_file):
    # Baseline data
    data = pd.read_csv(dataset_file, sep='\t', names=['chrom', 'pos','end','ref_allele', 'alt_allele', 'driver_stat'])
    data.drop(columns='end', inplace=True)

    # loop over file_paths and merge them on data

    for file in file_paths:
        splitted = file.split('/')
        if 'Conservation' in splitted:
            features = pd.read_csv(file, sep='\t')
            features = features.rename(columns={"start":"pos", "alt":"alt_allele","ref":"ref_allele"})
            print(features.head())
        else:
            features = pd.read_csv(file, sep='\t')

        # Cons scores pos is actually the end instead of the start.
        data = pd.merge(data, features, on=['chrom','pos','ref_allele','alt_allele'])
        # print(data.head())
    return data
    


if __name__ == '__main__':
    file_paths = locate_files('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features')
    dataset_file = '/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/CADD_Data/final_cadd_data_stat.sorted.bed'
    data = merge_files(file_paths, dataset_file)

    data.to_csv('Annotated_data.csv', sep=',', index=False)
