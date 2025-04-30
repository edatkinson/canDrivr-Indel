
import pandas as pd
import requests
import os
import sys
import numpy as np
from unipressed import IdMappingClient
import time
from bioservices import UniProt
from Bio import SeqIO
from io import StringIO
import concurrent.futures
from config import *
from collections import Counter

def get_protein_sequence(uniprot_id):
    """Retrieves a single protein sequence from UniProt."""
    uniprot_service = UniProt()
    try:
        sequence = uniprot_service.retrieve(str(uniprot_id), frmt='fasta')
        if sequence:
            record = SeqIO.read(StringIO(sequence), "fasta")
            return str(record.seq)
        else:
            print(f"No sequence found for UniProt ID {uniprot_id}")
            return None
    except Exception as e:
        print(f"Error retrieving sequence for {uniprot_id}: {e}")
        return None

def get_protein_sequences_with_positions(vep_dataset, max_workers=8):  # Added max_workers
    """Retrieves UniProt sequences with positions from VEP data efficiently."""

    df = pd.read_csv(vep_dataset, sep='\t')

    # Group by UniProt ID to avoid redundant requests
    grouped = df.groupby("uniprot_res")

    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor: #Use multithreading
        futures = {uniprot_id: executor.submit(get_protein_sequence, uniprot_id) for uniprot_id in grouped.groups}

        for uniprot_id, future in futures.items():
            sequence = future.result()  # Get the result (this will block until the request is done)
            if sequence:
              group_df = grouped.get_group(uniprot_id) #Get the rows with this uniprot_id
              for index, row in group_df.iterrows(): #Iterate over the rows
                results.append({
                    'gene': row['gene'],
                    'uniprot_id': uniprot_id,
                    'sequence': sequence,
                    'protein_position': row['protein_position'],
                    'chrom': row['chrom'],
                    'pos': row['pos'],
                    'ref_allele': row['ref_allele'],
                    'alt_allele': row['alt_allele'],
                    'driver_stat': row['driver_stat'] if row['driver_stat'] == 1 else 0
                })

    if results:
        return pd.DataFrame(results)
    else:
        return None



def getUniprotIDs(vep_dataset):
    """
    Extract UniProt IDs and related information from a VEP (Variant Effect Predictor) dataset.

    Parameters:
    - vep_dataset (str): Path to the VEP dataset file.

    Returns:
    - Tuple: A tuple containing two DataFrames - the first with UniProt IDs and the second with processed information.
    """
    df3 = pd.read_csv(vep_dataset, sep='\t')

    # Extract unique genes for UniProt ID conversion
    genes = df3["gene"].unique().tolist()

    # Convert the list of genes to a string with curly braces
    genes_str = "{" + ", ".join(f'"{gene}"' for gene in genes) + "}"

    # Retrieve UniProt IDs using IdMappingClient
    try:
        request = IdMappingClient.submit(source="GeneCards", dest="UniProtKB", ids={genes_str})
        time.sleep(10)
        uniprotIDs = pd.DataFrame(list(request.each_result()))
        uniprotIDs.columns = ["gene", "uniprot_res"]
    except EncodingWarning as e:
        print(f"Error in IdMappingClient: {e}")
        uniprotIDs = pd.DataFrame(columns=["gene", "uniprot_res"])

    # Return a tuple with UniProt IDs and the processed DataFrame
    return uniprotIDs, df3


def get_alpha_fold_atom2(uniprot_id):
    """
    Get AlphaFold Atom Information for a UniProt ID

    This function takes a UniProt ID as input and retrieves relevant information from a CIF file associated with the AlphaFold protein structure prediction. It involves the following steps:

    Error Handling:
    - Handles FileNotFoundError if the CIF file is not found for the given UniProt ID.
    - Handles pd.errors.EmptyDataError if there is no data in the CIF file.
    - Catches other unexpected exceptions and prints an error message with details.

    Parameters:
    - uniprot_id (str): UniProt ID of the protein.

    Returns:
    - pd.DataFrame: DataFrame containing relevant information from the CIF file.
    """

    # Select rows from DataFrame where "uniprot_res" column matches the provided uniprot_id
    res = df[df["uniprot_res"] == uniprot_id].reset_index(drop=True)
    
    # Construct the path to the CIF file based on the uniprot_id
    cif_file_path = f"{alpha_fold_cif_location}/AF-{uniprot_id}-F1-model_v4.cif.gz"
    
    if os.path.isfile(cif_file_path):
        # Read the CIF file into a DataFrame with tab-separated values and no header
        cif_file = pd.read_csv(cif_file_path, sep="\t", header=None)

        # Initialise an empty DataFrame to store the results
        results = pd.DataFrame()

        if len(res) == 1:

            position_range_str = res.loc[0, "protein_position"]
            prot_pos = position_range_str

            # Check if it's a range
            if "-" in position_range_str:
                start_pos, end_pos = map(int, position_range_str.split("-")) #Extract Start and End positions
                mid_point = ((start_pos + end_pos) // 2)
                prot_pos = mid_point

            # If there is only one matching row in the DataFrame
            atom_mask = cif_file[0].str.contains("ATOM")  # Select rows with "ATOM"
            filtered_cif = cif_file[atom_mask][0]  # Extract column 0 from filtered rows
            position_mask = filtered_cif[filtered_cif.str.contains(str(prot_pos))]  # Select rows with protein position
            if not position_mask.empty:
                results = position_mask.str.split("?", expand=True)[1].to_string()  # Split and select the relevant column
                results = pd.DataFrame(" ".join(results.split()).split()).transpose().iloc[:, 1:6]
                results = pd.concat([res[["chrom", "pos", "ref_allele", "alt_allele", "driver_stat"]], results], axis=1)
            else:
                print("In (1), not found")
                print(f"Protein Position: {prot_pos} not found at {cif_file_path}")

        else:
            # If there are multiple rows with different protein positions
            results = []
            for position_range_str in res["protein_position"].tolist():
                # Check if it's a range
                if "-" in position_range_str:
                    split_pos = position_range_str.split("-")
                    if any(not s.isdigit() for s in split_pos):
                        print(f"Warning: Non-digit characters found in position range: {position_range_str}")
                        continue  # Skip to the next position_range_str
                    start_pos, end_pos = map(int, position_range_str.split("-")) #Extract Start and End positions
                    mid_point = ((start_pos + end_pos) // 2)
                    prot_pos = mid_point
                    # print(f'{start_pos}-{end_pos}: {mid_point}')
                else:
                    if not position_range_str.isdigit():
                        print(f"Warning: Non-digit character found in position: {position_range_str}")
                        continue  # Skip to the next position_range_str
                    prot_pos = position_range_str
                
                res_filtered = res[res["protein_position"] == position_range_str]
                atom_mask = cif_file[0].str.contains("ATOM")
                filtered_cif = cif_file[atom_mask][0]
                position_mask = filtered_cif[filtered_cif.str.contains(str(prot_pos))]
                if not position_mask.empty:
                    cif_res = position_mask.str.split("?", expand=True)[1].to_string()
                    cif_res = pd.DataFrame(" ".join(cif_res.split()).split()).transpose().iloc[:, 1:6]
                    cif_res = pd.concat([res_filtered[["chrom", "pos", "ref_allele", "alt_allele", "driver_stat"]], cif_res], axis=1)
                    results.append(cif_res)
                else:
                    print(f"Protein Position: {prot_pos}: not found at {uniprot_id}")
                    continue

            # Concatenate the results from multiple positions and drop any NaN values
            results = pd.concat(results)
        
        results = results.rename(columns = {1: "X_coordinate", 2: "Y_coordinate", 3: "Z_coordinate", 4: "atom_site_occupancy", 5: "isotropic_temperature"})
        results['driver_stat'].replace('.', 0, inplace=True)

        # Return the final DataFrame containing the results
        return results


def get_alpha_fold_atom(uniprot_id):
    """
    Get AlphaFold Atom Information for a UniProt ID

    This function takes a UniProt ID as input and retrieves relevant information from a CIF file associated with the AlphaFold protein structure prediction. It involves the following steps:

    Error Handling:
    - Handles FileNotFoundError if the CIF file is not found for the given UniProt ID.
    - Handles pd.errors.EmptyDataError if there is no data in the CIF file.
    - Catches other unexpected exceptions and prints an error message with details.

    Parameters:
    - uniprot_id (str): UniProt ID of the protein.

    Returns:
    - pd.DataFrame: DataFrame containing relevant information from the CIF file.
    """

    # # Select rows from DataFrame where "uniprot_res" column matches the provided uniprot_id
    res = df[df["uniprot_res"] == uniprot_id].reset_index(drop=True)

    # # Construct the path to the CIF file based on the uniprot_id
    cif_file_path = f"{alpha_fold_cif_location}/AF-{uniprot_id}-F1-model_v4.cif.gz"
    
    if os.path.isfile(cif_file_path):
        # Read the CIF file into a DataFrame with tab-separated values and no header
        cif_file = pd.read_csv(cif_file_path, sep="\t", header=None)

        # Initialise an empty DataFrame to store the results
        results = pd.DataFrame()

        if len(res) == 1:
            position_range_str = res.loc[0, "protein_position"]

            # Check if it's a range
            if "-" in position_range_str:
                start_pos, end_pos = map(int, position_range_str.split("-")) #Extract Start and End positions
                positions_in_range = list(range(start_pos, end_pos + 1)) # Create a list of positions in the range
            else:
                positions_in_range = [int(position_range_str)] #If not a range, make it a list with one element

            all_results_for_range = []

            for position in positions_in_range: #Iterate over the positions in the range
                atom_mask = cif_file[0].str.contains("ATOM")
                filtered_cif = cif_file[atom_mask][0]
                position_mask = filtered_cif[filtered_cif.str.contains(str(position))] # Convert position to string

                if not position_mask.empty:
                    cif_res_split = position_mask.str.split("?", expand=True)
                    if not cif_res_split.empty:
                        cif_res = cif_res_split[1].to_string()
                        cif_res = pd.DataFrame(" ".join(cif_res.split()).split()).transpose().iloc[:, 1:6]
                        cif_res.columns = ['X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature']
                        cif_res = pd.concat([res[["chrom", "pos", "ref_allele", "alt_allele","driver_stat"]], cif_res], axis=1)
                        all_results_for_range.append(cif_res)
                    else:
                        print(f"Warning: No data after split for uniprot_id: {uniprot_id}, position: {position}")
                        empty_df = pd.DataFrame(columns=['chrom', 'pos', 'ref_allele', 'alt_allele','driver_stat', 'X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature'])
                        all_results_for_range.append(empty_df)

                else:
                    print(f"Warning: No matching position data found for uniprot_id: {uniprot_id}, position: {position}")
                    empty_df = pd.DataFrame(columns=['chrom', 'pos', 'ref_allele', 'alt_allele','driver_stat', 'X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature'])
                    all_results_for_range.append(empty_df)


            if all_results_for_range:
                results = pd.concat(all_results_for_range, ignore_index=True)
            else:
                results = pd.DataFrame(columns=['chrom', 'pos', 'ref_allele', 'alt_allele','driver_stat', 'X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature'])

        else: #Multiple positions, handles ranges as well
            results = []
            for position_range_str in res["protein_position"].tolist():
                
                if "-" in position_range_str:
                    start_pos, end_pos = map(int, position_range_str.split("-"))
                    positions_in_range = list(range(start_pos, end_pos + 1))
                else:
                    positions_in_range = [int(position_range_str)]

                for position in positions_in_range:
                    res_filtered = res[res["protein_position"] == position_range_str] #res_filtered now contains the whole range row
                    atom_mask = cif_file[0].str.contains("ATOM")
                    filtered_cif = cif_file[atom_mask][0]
                    position_mask = filtered_cif[filtered_cif.str.contains(str(position))]

                    if not position_mask.empty:
                        cif_res_split = position_mask.str.split("?", expand=True)
                        if not cif_res_split.empty:
                            cif_res = cif_res_split[1].to_string()
                            cif_res = pd.DataFrame(" ".join(cif_res.split()).split()).transpose().iloc[:, 1:6]
                            
                            cif_res.columns = ['X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature']
        
                            res_filtered =res_filtered.reset_index(drop=True)
                            cif_res_new = pd.concat([res_filtered[["chrom", "pos", "ref_allele", "alt_allele", "driver_stat"]], cif_res], axis=1)
                            # print(cif_res)
                            results.append(cif_res_new)
                            # print(results[-1])
                            
                        else:
                            print(f"Warning: No data after split for uniprot_id: {uniprot_id}, position: {position}")
                            empty_df = pd.DataFrame(columns=['chrom', 'pos', 'ref_allele', 'alt_allele',"driver_stat", 'X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature'])
                            results.append(empty_df)
                    else:
                        print(f"Warning: No matching position data found for uniprot_id: {uniprot_id}, position: {position}")
                        empty_df = pd.DataFrame(columns=['chrom', 'pos', 'ref_allele', 'alt_allele',"driver_stat", 'X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature'])
                        results.append(empty_df)
                        

            if results:
                results = pd.concat(results, ignore_index=True)
            else:
                results = pd.DataFrame(columns=['chrom', 'pos', 'ref_allele', 'alt_allele','driver_stat', 'X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature'])

                # Concatenate the results from multiple positions and drop any NaN values
                # results = pd.concat(results).dropna()
            
            
            results = results.rename(columns = {1: "X_coordinate", 2: "Y_coordinate", 3: "Z_coordinate", 4: "atom_site_occupancy", 5: "isotropic_temperature"})

        # Return the final DataFrame containing the results
        return results


def average_vals(atom_file):
    data = pd.read_csv(atom_file, sep='\t', index_col=None)
    data.reset_index(drop=True, inplace=True)
    data.drop(columns=['Unnamed: 0'], inplace=True, errors='ignore')  # Remove the "Unnamed: 0" column. errors='ignore' will not raise an exception if the column is not found

    cols_to_avg = ['X_coordinate', 'Y_coordinate', 'Z_coordinate', 'atom_site_occupancy', 'isotropic_temperature']

    # 1. Identify and average duplicates:
    averaged_data = data.groupby('pos')[cols_to_avg].mean().reset_index()

    # 2. Identify rows to keep (either averaged or original if not duplicated):
    rows_to_keep = ~data.duplicated(['pos','ref_allele','alt_allele'], keep='first') #Mark the first occurence of a position as True, all others as False

    # 3. Filter the original data, then merge the averaged data
    data_unique = data[rows_to_keep].copy() #Create a copy so that the original dataframe is not altered
    data_unique = data_unique.drop(columns=cols_to_avg, errors='ignore') #Drop the original coordinate columns, ignore if they don't exist

    data_final = pd.merge(data_unique, averaged_data, on='pos', how='left')
    # data_final.dropna(inplace=True)
    data_final = data_final.reset_index(drop=True)
    data_final.replace('.', int(0), inplace=True)

    return data_final


def get_alpha_fold_struct(uniprot_id):
    """
    Get AlphaFold Atom Information for a UniProt ID

    This function takes a UniProt ID as input and retrieves relevant information from a CIF file associated with the AlphaFold protein structure prediction. It involves the following steps:

    Error Handling:
    - Handles FileNotFoundError if the CIF file is not found for the given UniProt ID.
    - Handles pd.errors.EmptyDataError if there is no data in the CIF file.
    - Catches other unexpected exceptions and prints an error message with details.

    Parameters:
    - uniprot_id (str): UniProt ID of the protein.

    Returns:
    - pd.DataFrame: DataFrame containing relevant information from the CIF file.
    """

    # Select rows from DataFrame where "uniprot_res" column matches the provided uniprot_id
    res = df[df["uniprot_res"] == uniprot_id].reset_index(drop=True)

    # Construct the path to the CIF file based on the uniprot_id
    cif_file_path = f"{alpha_fold_cif_location}/AF-{uniprot_id}-F1-model_v4.cif.gz"
    
    if os.path.isfile(cif_file_path):
        # Read the CIF file into a DataFrame with tab-separated values and no header
        cif_file = pd.read_csv(cif_file_path, sep="\t", header=None)

        # Initialise an empty DataFrame to store the results
        results = pd.DataFrame()

        if len(res) == 1:
            # If there is only one matching row in the DataFrame
            atom_mask = cif_file[0].str.contains("A ")  # Select rows with "ATOM"
            filtered_cif = cif_file[atom_mask][0]   # Extract column 0 from filtered rows

            atom_mask = filtered_cif.str.contains("\?") # Select rows with "ATOM"
            filtered_cif = filtered_cif[atom_mask]

            atom_mask = filtered_cif.str.contains("ATOM")
            filtered_cif = filtered_cif[~atom_mask]
            

            filtered_cif = filtered_cif[filtered_cif.str.contains(res.loc[0, "protein_position"])]


            filtered_cif = filtered_cif.str.split(expand = True).dropna()

            if len(filtered_cif) != 0:
                if 13 in filtered_cif.columns:
                    filtered_cif = filtered_cif.loc[filtered_cif[2] == res.loc[0, "protein_position"], [6, 13]].reset_index(drop = True)

                    filtered_cif.columns = ["struct_conf_type", "struct_conf_type_ID"]

                    results = pd.concat([res[["chrom", "pos", "ref_allele", "alt_allele"]].reset_index(drop = True), filtered_cif], axis=1)

        else:
            # If there are multiple rows with different protein positions
            results = []
            for position in res["protein_position"].tolist():
                res2 = res[res["protein_position"] == position]

                # If there is only one matching row in the DataFrame
                atom_mask = cif_file[0].str.contains("A ")  # Select rows with "A "
                filtered_cif = cif_file[atom_mask][0]   # Extract column 0 from filtered rows

                atom_mask = filtered_cif.str.contains("\?")
                filtered_cif = filtered_cif[atom_mask]
                atom_mask = filtered_cif.str.contains("ATOM")
                filtered_cif = filtered_cif[~atom_mask]
                filtered_cif = filtered_cif[filtered_cif.str.contains(position)]

                filtered_cif = filtered_cif.str.split(expand = True).dropna()

                if len(filtered_cif) != 0:
                    if 13 in filtered_cif.columns:
                        filtered_cif = filtered_cif.loc[filtered_cif[2] == position, [6, 13]].reset_index(drop = True)

                        filtered_cif.columns = ["struct_conf_type", "struct_conf_type_ID"]
                        res3 = pd.concat([res2[["chrom", "pos", "ref_allele", "alt_allele"]].reset_index(drop = True), filtered_cif], axis=1)
                        if len(res3) != 0:
                            results.append(res3)
            
            if len(results) > 1:
                # Concatenate the results from multiple positions and drop any NaN values
                results = pd.concat(results).dropna()
            elif len(results) == 1:
                results = results[0]
            else:
                results = results
        if type(results) is list:
            x="y"
        else:
            # Return the final DataFrame containing the results
            return results.dropna()



#Combine vepAA of cos and benign indels. Then merge these with sequence data. Then for the specified position of each mutation, make the correct changes to the sequence and update for full Mutated sequence.
#Then input these sequences into the Meta ESM api to get the results.

if __name__ == '__main__':
    alpha_fold_cif_location = alpha_fold_files

    # uniprotIDs, df3 = getUniprotIDs('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/ProteinStructure/vep_prot.tsv')
    # df = pd.merge(df3, uniprotIDs, on = "gene", how = "right")
    # df.to_csv('vep_prot.tsv',sep='\t', index=None)
    # results = get_protein_sequences_with_positions('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/ProteinStructure/vep_prot.tsv')
    # results.to_csv('vep_prot_2.tsv',sep='\t', index = None)
    df = pd.read_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/ProteinStructure/vep_indtest_prot.tsv', sep='\t')
    print(df.head())

    all_alphafold = [get_alpha_fold_atom2(uniprot_id) for uniprot_id in df["uniprot_res"]]
    all_alphafold = pd.concat(all_alphafold)
    # # print(all_alphafold)
    all_alphafold.to_csv('atomic_coordinates_indtest.tsv', sep='\t', index=None)
    # df = pd.read_csv('vep_prot.tsv', sep='\t')
    # struct_alphafold = [get_alpha_fold_struct(uniprot_id) for uniprot_id in df["uniprot_res"].unique()]
    # struct_alphafold = pd.concat(struct_alphafold)
    # # Get one hot encoding of columns B
    # struct_conf_type = pd.get_dummies(struct_alphafold['struct_conf_type'])
    
    # # Get one hot encoding of columns B
    # struct_conf_type_ID = pd.get_dummies(struct_alphafold['struct_conf_type_ID'])

    # # Drop column B as it is now encoded
    # struct_alphafold = struct_alphafold.drop(['struct_conf_type', 'struct_conf_type_ID'],axis = 1)

    # # Join the encoded df
    # struct_alphafold = pd.concat([struct_alphafold, struct_conf_type, struct_conf_type_ID], axis = 1)
    # struct_alphafold.to_csv("alphafold_structural.bed", sep = "\t", index = None)

    # data = average_vals('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/ProteinStructure/atomic_coordinates_indtest.tsv')
    # data.to_csv('atomic_average_data.tsv', sep='\t', index = None)
    

