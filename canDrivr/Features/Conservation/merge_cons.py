import pandas as pd
import os


def merge_on_common_columns(directory, base_file, output_file):
    # Define common columns
    common_cols = ['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat']
    
    # Read base file
    base_data = pd.read_csv(base_file, sep=',')  # Adjust separator if needed
    base_data.drop(columns=['stop'], errors='ignore', inplace=True)  # Ignore errors if 'stop' is missing
    base_data.rename(columns={'start': 'pos'}, inplace=True)
    base_data['pos'] = base_data['pos'].astype(str)  # Ensure position is a string for merging
    base_data.drop_duplicates(subset=common_cols, inplace=True)  # Remove duplicates before merge
    
    # Get list of CSV files in directory
    list_files = [f for f in os.listdir(directory) if f.endswith('.bedGraph')]
    
    for f in list_files:
        file_path = os.path.join(directory, f)
        df = pd.read_csv(file_path, sep='\t')  # Adjust separator if needed
        df['pos'] = df['pos'].astype(str)
        df.drop_duplicates(subset=common_cols, inplace=True)  # Remove duplicates in input file
        
        # Find common columns to merge on
        cols_to_merge = [col for col in common_cols if col in df.columns]
        
        if not cols_to_merge:
            print(f"Skipping {f}, no common columns found.")
            continue
        
        print(f"Merging {f} on {cols_to_merge}")
        base_data = pd.merge(base_data, df, on=cols_to_merge, how='outer', suffixes=('', f'_{f}'))
    
    # Drop duplicates after merging
    base_data.drop_duplicates(subset=common_cols, inplace=True)
    
    # Save merged file
    base_data.to_csv(output_file, index=False)
    print(f"Merged file saved at: {output_file}")

# Example usage
merge_on_common_columns(
    directory='/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/indtest',  # Change to your directory
    base_file='/Users/edatkinson/Repos/canDrivr-Indel/independentTesting/independent_test_set_combined.csv',  # Change to your base file
    output_file='/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/Conservation/merged_cons_indtest.csv'  # Change to your output path
)