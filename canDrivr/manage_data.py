import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter

def gnomad_to_tsv(file):
    data = pd.read_csv(file, sep='\t')
    data = data.rename(columns={'SYMBOL': 'gene', 'end':'stop', 'reference_name':'chrom', })
    
    data = data[['chrom', 'start', 'stop', 'ref_allele', 'alt_allele', 'gene']]

    data['driver_stat'] = [0 for _ in range(len(data))]
    data.to_csv('gnomad_all_chrom_filtered.tsv', sep='\t', index=None)


def filter_cosmic(file):
    cos = pd.read_csv(file, sep='\t')
    cos = cos.iloc[:, [0,14, 15, 16, 23, 24]]
    cos.columns = ['gene', 'chrom', 'start', 'stop', 'ref_allele','alt_allele']
    cos.replace(np.nan, '-', regex=True, inplace=True)
    try:
        cos['count'] = cos.groupby(['start', 'stop', 'ref_allele', 'alt_allele'])['start'].transform('size')
    except Exception as e:
        print("Error in groupby operation:", e)
        return None
    # Drop duplicate rows, keeping the first occurrence
    cos = cos.drop_duplicates(subset=['start', 'stop', 'ref_allele', 'alt_allele'], keep='first')
    cos = cos.reset_index(drop=True)

    gene = cos.iloc[:, [0]]  # Select the first column (gene names)

    # Count gene occurrences
    gene_counts = gene['gene'].value_counts()

    # Convert to DataFrame and reset index for sorting
    unique_genes_df = gene_counts.reset_index()
    unique_genes_df.columns = ['gene_name', 'count'] # Rename columns for clarity

    # Sort by count in descending order
    sorted_genes_df = unique_genes_df.sort_values(by='count', ascending=False)
    
    top_genes = sorted_genes_df[sorted_genes_df['count'] > 50]['gene_name'].tolist()
    
    
    # Top genes with mutaiton count > 50
    cos = cos[cos.iloc[:,0].isin(top_genes)].reset_index(drop=True)
    cos['driver_stat'] = [1 for _ in range(len(cos['chrom']))]
    print(Counter(cos['chrom']))
    # Reorder the columns
    cos = cos[['chrom', 'start', 'stop', 'ref_allele', 'alt_allele','gene', 'driver_stat']]
    return cos



def match_data(gnomad, cos):
    """
    No matches between datasets, this could be because gnomad isn't restricted to the genes with indel counts above 300 like cosmic.
    """
    gnomad_data = pd.read_csv(gnomad, sep='\t')
    cos_data = pd.read_csv(cos, sep='\t')

    cols = ['chrom','start', 'stop','ref_allele', 'alt_allele']

    # Merge df1 with df2 based on the specified columns
    inner_merged = gnomad_data.merge(cos_data, on=cols, how='inner')
    print("Rows present in both (using merge):\n", inner_merged)


    # Compare columms of gnomad start and stop and print matches


def alt_allele_maker(cos):
    data = pd.read_csv(cos, sep='\t')
    ref_fasta = "/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Features/FG5_gc_CpG/hg38.fa"  # Path to your reference FASTA file
    data['chrom'] = data['chrom'].apply(lambda x: 'chr'+str(x))
    vcf = data[['chrom', 'start', 'stop', 'ref_allele', 'alt_allele']]
    

    with open(ref_fasta) as handle:
        ref_sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta")) # Load the FASTA into a dictionary for faster lookup

    vcf['original_ref'] = ''  # Initialize a new column to store the original reference sequence

    for index, row in vcf.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['stop']
        ref_allele = row['ref_allele']
        alt_allele = row['alt_allele']

        if chrom in ref_sequences: # Check if chrom exists in the fasta file
            ref_record = ref_sequences[chrom]
            if alt_allele == '-':
                original_ref = ref_record.seq[start-2:end]  # Extract the original reference sequence (0-based indexing)
                if isinstance(original_ref, type(np.nan)):
                    vcf.drop(index, inplace=True)
                    continue
                else:
                    vcf.loc[index, 'original_ref'] = str(original_ref).upper() # store sequence in dataframe, convert from seq object to string
                    vcf.loc[index, 'alt_allele'] = str(original_ref[0]).upper()
                    vcf.loc[index, 'ref_allele'] = str(original_ref).upper()
                # vcf.loc[index, 'start'] = start -1
                if start == end:
                    vcf.loc[index, 'start'] = start -1 #Ask Amy about this.
            elif ref_allele == '-':
                if start > end:
                    vcf.drop(index, inplace=True)
                    continue
                else:
                    original_ref = ref_record.seq[start:end]
                    vcf.loc[index, 'original_ref'] = str(original_ref).upper()
                    vcf.loc[index, 'ref_allele'] = str(original_ref).upper()
                    vcf.loc[index, 'alt_allele'] = str(original_ref).upper()+alt_allele.upper()

        else:
            vcf.loc[index, 'original_ref'] = "Chromosome not found" # Handle cases where the chromosome is not found in FASTA

    vcf['driver_stat'] = [1 for _ in range(len(vcf))]
    vcf = vcf.drop(columns=['original_ref'])
    vcf = vcf.dropna()

    print(vcf)
    vcf.to_csv("modified_2" + cos, sep='\t', index=False) # Save the modified data to a new file




def tsv_to_bed(gnomad):
    data = pd.read_csv(gnomad, sep='\t')
    data = data[['chrom', 'start', 'stop', 'ref_allele', 'alt_allele']]
    data['chrom'] = data['chrom'].apply(lambda x: 'chr'+str(x))
    # data.drop(columns=['ref_allele', 'alt_allele'])
    data.to_csv('gnomad_formatted.bed', sep='\t', index=None)



def convert_to_vcf(input_file, output_vcf):
    df = pd.read_csv(input_file, sep='\t')

    # Add missing VCF columns (with placeholders)
    df['ID'] = '.'  # or generate IDs if you have a naming scheme
    df['QUAL'] = '.'  # or calculate quality scores if available
    df['FILTER'] = 'PASS'  # or add filter information
    df['INFO'] = '.'  # or add any additional info

    # Reorder columns to match VCF format
    vcf_columns = ['chrom', 'start', 'ID', 'ref_allele', 'alt_allele', 'QUAL', 'FILTER', 'INFO']
    df = df[vcf_columns]
    df.rename(columns={'chrom': '#CHROM', 'start': 'POS'}, inplace=True) # rename columns to #CHROM and POS

    # Write to VCF file (add header)
    with open(output_vcf, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')  # Or your VCF version
        f.write('#' + '\t'.join(vcf_columns) + '\n')
        df.to_csv(f, sep='\t', index=False, header=False)

# # Example usage:

def modify_gnomd_hg38(gnomad):
    data = pd.read_csv(gnomad, sep='\t', names=['chrom', 'start', 'stop', 'ref_allele', 'alt_allele'])

    data['driver_stat'] = [0 for _ in range(len(data))]

    data.to_csv('benign_stat.tsv', sep='\t', index=None, header=True)







if __name__ == '__main__':
    # cos = 'chr21.tsv'
    cos = ''
    alt_allele_maker(cos)
    
#     # gnomad_to_tsv('gnomad_all_chrom.tsv')
#     # match_data('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/gnomad/benign_stat.tsv', 'filtered_cosmic.tsv')


#     # filter_cosmic('new_cosmic.tsv').to_csv('filtered_cosmic.tsv', sep='\t', index=None)

#     # alt_allele_maker('filtered_cosmic.tsv')
#     # tsv_to_bed('gnomad_all_chrom_filtered.tsv')
#     # convert_to_vcf("gnomad_formatted.bed", "gnomad.vcf")
#     modify_gnomd_hg38('gnomad_hg38.bed')
