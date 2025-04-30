import pandas as pd
import numpy as np
from pubmed import get_pubmed_counts
import json    
from sklearn.preprocessing import MinMaxScaler


# def vep_to_df(vep_file, ml_data):
#     df = pd.read_csv(vep_file, sep='\t',low_memory=False)
#     df.sort_values(by=['POS'], ascending=True, inplace=True)
#     ml_df = pd.read_csv(ml_data, sep=',', low_memory=False)
#     ml_df.sort_values(by=['pos'], ascending=True, inplace=True)
#     ml_df.reset_index(drop=True)

#     genes = df['INFO'].str.split('|', expand=True)[3]
#     names = df['INFO'].str.split('|', expand=True)[1]
#     df['genes'] = genes
#     df['names'] = names


#     df.drop(columns=['QUAL', 'FILTER'], inplace=True)
#     df.rename(columns={'POS':'pos', 'ID':'end', 'REF':'ref_allele', 'ALT':'alt_allele', '#CHROM':'chrom'}, inplace=True)

#     sliced_df = df[['chrom','pos','ref_allele', 'alt_allele','genes','names']]

#     merged = ml_df.merge(sliced_df, on=['chrom','pos','ref_allele', 'alt_allele'])
#     allowed = ['frameshift_variant', 'inframe_insertion', 'inframe_deletion'] # and any variations

#     #Select data which include inframe and frameshift only
#     new = merged[merged['names'].isin(allowed)]
    
#     # split into inframe and frameshift ds
#     inframe = new[new['names'].isin(['inframe_insertion','inframe_deletion'])]
#     frameshift = new[~new['names'].isin(['inframe_insertion','inframe_deletion'])]
#     print(frameshift.shape, inframe.shape)
#     # other = ml_df[~ml_df['names'].isin(allowed)]
    
#     # inframe.drop(columns=['merged', 'types'], inplace=True)
#     # frameshift.drop(columns=['merged', 'types'], inplace=True)
#     # other.drop(columns=['merged', 'types'], inplace=True)

#     with open('new_approach/frameshift_genes.json','r') as f:
#         pmed_count = json.load(f)
    

#     frameshift['pmed_count'] = frameshift['genes'].map(pmed_count)


#     # print(len(inframe[inframe['driver_stat'] == 1]['chrom'].value_counts())) # missing chr10
#     # print(frameshift[frameshift['driver_stat'] == 1]['chrom'].value_counts()) # missing chr21 as always

#     # inframe.to_csv('data/inframe_data.csv', sep=',', index=None)
#     frameshift.to_csv('data/pmed_frameshift.csv', sep=',', index=None)
#     # other.to_csv('data/other_data.csv', sep=',', index=None)
#     return None
#     # return vep_frameshift


def frameshift_pmed(frameshift_ml, vep):
    df = pd.read_csv(frameshift_ml, sep=',')
    vep = pd.read_csv(vep, sep='\t')
    vep['genes'] = vep['INFO'].str.split('|', expand=True)[3]
    vep.rename(columns={'POS':'pos', 'ID':'end', 'REF':'ref_allele', 'ALT':'alt_allele', '#CHROM':'chrom'}, inplace=True)
    vep.drop(columns=['INFO', 'end', 'FILTER', 'QUAL'], inplace=True)
    merged = df.merge(vep, on = ['chrom', 'pos', 'ref_allele','alt_allele'])

    print(np.unique(merged['genes']).tolist())
    pmed = get_pubmed_counts(np.unique(merged['genes']).tolist())

    with open('/Users/edatkinson/Repos/Modelling/new_approach/inframe_genes_more.json', 'w') as f:
        json.dump(pmed, f)
        # pmed = json.load(f)

    pmed[''] = 0

    merged['pmed'] = merged['genes'].map(pmed)
    merged['pmed'].replace(0,1, inplace=True)
    merged['pmed'] = np.log10(merged['pmed'])

    merged.to_csv('/Users/edatkinson/Repos/Modelling/data/pmed_log_inframe.csv', sep=',', index=None)
    return merged

frameshift_pmed('/Users/edatkinson/Repos/Modelling/data/inframe_data.csv', '/Users/edatkinson/Repos/Modelling/all_variant_effect_output_all.txt')

# vep_to_df('/Users/edatkinson/Repos/Modelling/all_variant_effect_output_all.txt', '/Users/edatkinson/Repos/Modelling/merged_annotated_data.csv')

def combine_vep(file1, file2):
    df1 = pd.read_csv(file1, sep='\t', skiprows=4)
    df2 = pd.read_csv(file2, sep='\t', skiprows=4)
    
    both = pd.concat([df1, df2])

    names = both['INFO'].str.split('|', expand=True)[1]
    both['names'] = names

    frameshift = both[both['names'].isin(['frameshift_variant'])]
    inframe = both[both['names'].isin(['inframe_deletion','inframe_insertion'])]

    frameshift.drop(columns=['names'], inplace=True)
    inframe.drop(columns=['names'], inplace=True)
    # frameshift.to_csv('data/vep_frameshift.txt', sep='\t', index=None)
    # inframe.to_csv('data/vep_inframe.txt', sep='\t', index=None)
    genes = both['INFO'].str.split('|', expand=True)[3]

    # print(genes.value_counts())
    return np.unique(genes)



def combine_pubmed(file1, file2, frameshift_ml, inframe_ml):

    """
    File1: Frameshift
    File2: Inframe
    Genes: Gene Pubmed Counts
    Frameshift_ml : annotated ml data
    """

    df1 = pd.read_csv(file1, sep='\t', skiprows=4)
    df2 = pd.read_csv(file2, sep='\t', skiprows=4)


    df1['pmed'] = df1['INFO'].str.split('|', expand=True)[3] # genes

    df2['pmed'] = df2['INFO'].str.split('|', expand=True)[3] # genes

    # genes_df1 = get_pubmed_counts(np.unique(df1['pmed']))
    # with open('frameshift_genes.json', 'w') as f: 
    #     json.dump(genes_df1, f)
    #     print('Frameshift Json Genes made! \n')

    with open('new_approach/frameshift_genes.json', 'r') as f:
        loaded_dict = json.load(f)
    
    genes_df1 = loaded_dict
    


    # genes_df2 = get_pubmed_counts(np.unique(df2['pmed']))
    # with open('inframe_genes.json', 'w') as f: 
    #     json.dump(genes_df1, f)
    #     print('Inframe Json Genes made! \n')


    df1['pmed'] = df1['pmed'].map(genes_df1)
    # df2['pmed'].map(genes_df2)

    filtered_df1 = df1[['#CHROM', 'POS', 'REF','ALT', 'pmed']]
    filtered_df1.rename(columns={'#CHROM':'chrom', 'POS':'pos', 'REF':'ref_allele', 'ALT':'alt_allele'}, inplace=True)
    filtered_df1 = filtered_df1[~filtered_df1['chrom'].isin(['chrX','chrY'])]

    filtered_df1.sort_values(by=['chrom', 'pos'], inplace=True)
    


    # filtered_df2 = df2[['#CHROM', 'POS', 'REF','ALT', 'pmed']]
    # filtered_df2.rename(columns={'#CHROM':'chrom', 'POS':'pos', 'REF':'ref_allele', 'ALT':'alt_allele'}, inplace=True)
    # print(iltered_df1)

    frameshift_df = pd.read_csv(frameshift_ml, sep=',')
    frameshift_df.sort_values(by=['chrom', 'pos'], inplace=True)
    print(frameshift_df, '\n')
    print(filtered_df1, '\n')
    frameshift_df['pmed'] = df1['pmed']
    
    # inframe_df = pd.read_csv(inframe_ml, sep=',')

    # frameshift_df = frameshift_df.merge(filtered_df1, on=['chrom','pos', 'ref_allele', 'alt_allele'], how = 'right')
    # inframe_df.merge(filtered_df2, on=['chrom','pos','ref_allele', 'alt_allele'], how = 'right')
    
    # frameshift_df.to_csv('data/pmed_frameshift.csv', sep='\t', index=None)
    # inframe_df.to_csv('data/pmed_inframe.csv', sep='\t', index=None)

    


#combine_pubmed('/Users/edatkinson/Repos/Modelling/data/vep_frameshift.txt', '/Users/edatkinson/Repos/Modelling/data/vep_inframe.txt', '/Users/edatkinson/Repos/Modelling/data/frameshift_data.csv', '/Users/edatkinson/Repos/Modelling/data/inframe_data.csv')


# def new_seperator(data, get='frameshift_variant'):

#     ### ------ ### get frameshift

#     data['var_type'] = data['INFO'].str.split('|', expand=True)[1]
    
#     data['gene'] = data['INFO'].str.split('|', expand=True)[3]

#     select_data = data[data['var_type'].isin([get])]

#     print(select_data)

#     rename = {'#CHROM':'chrom', 'POS':'pos', 'REF':'ref_allele', 'ALT':'alt_allele'}
#     select_benign.rename(columns=rename, inplace=True) 
#     select_pathogenic.rename(columns=rename, inplace=True)


#     # genes = [select_pathogenic['gene'], select_benign['gene']]
#     # unique_genes = np.unique(genes)

#     ml = pd.read_csv(ml, sep=',')
#     matcher = ml[['chrom','pos', 'ref_allele', 'alt_allele']]

#     benign_match = matcher.merge(select_benign[['chrom','pos', 'ref_allele', 'alt_allele']], on=['chrom','pos', 'ref_allele', 'alt_allele'])
#     pathogenic_match = matcher.merge(select_pathogenic[['chrom','pos', 'ref_allele', 'alt_allele']], on=['chrom','pos', 'ref_allele', 'alt_allele'])

#     print(benign_match.shape)
#     print(pathogenic_match.shape)



# new_seperator('modified_2filtered_cosmic.tsv_variant_effect_output_all.txt', '/Users/edatkinson/Repos/Modelling/benign_stat.tsv_variant_effect_output_all.txt', '/Users/edatkinson/Repos/Modelling/merged_annotated_data.csv')


# df = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_frameshift.csv', sep=',')

# print(df[df['driver_stat']==0])

# print(df[df['driver_stat']==1])

# df = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_inframe.csv')


