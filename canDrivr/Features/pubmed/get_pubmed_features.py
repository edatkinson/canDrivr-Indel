import sys
sys.path.append('/Users/edatkinson/Repos/Modelling')
from pubmed import get_pubmed_counts
import pandas as pd
import json
import numpy as np


def get_gene_counts(vep):
    """File = vep_indtest_prot"""
    # df = pd.read_csv(file, sep='\t') 
    # genes = df['gene']
    vep = pd.read_csv(vep, sep='\t', skiprows=4)
    vep['genes'] = vep['INFO'].str.split('|', expand=True)[3]
    genes = vep['genes']
    vep.rename(columns={'POS':'pos', 'ID':'end', 'REF':'ref_allele', 'ALT':'alt_allele', '#CHROM':'chrom', 'QUAL':'driver_stat'}, inplace=True)
    vep.drop(columns=['INFO', 'end', 'FILTER'], inplace=True)
    df = vep[['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat']]
    df['genes'] = genes
    # counts = get_pubmed_counts(np.unique(genes).tolist())

    with open('indTestPubmedCounts.json', 'r') as f:
        # json.dump(counts, f)
        counts = json.load(f)

    df['pmed'] = df['genes'].map(counts)

    df['pmed'].replace(0,1, inplace=True)
    df['pmed'] = np.log10(df['pmed'])

    df.to_csv('pmed_features.csv', sep=',', index=None)
    return df

df = get_gene_counts('/Users/edatkinson/Repos/canDrivr-Indel/independentTesting/independent_test_set_combined_copy.csv_variant_effect_output_all.txt')
