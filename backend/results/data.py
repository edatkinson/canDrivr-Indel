
import pandas as pd
import numpy as np
import random
import sys
import requests
import json
#Row 550 missing chm1 and chm5 features up until row 1628, check drivr base for the missing features
#Interesting features: Rows 1635 onwards -
#VEP: Splice_acceptor_variant - Sequence_variant


'''
Need to get all rows which have a mutation decription = 'frameshift_variant', 'inframe_insertion', 'inframe_deletion'
From https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences, found this link in the comsmic data read me files under: 'Mutation Description'

TODO: Get the full dataset and filter out the rows with the above mutation descriptions.
TODO: Filter out all repeat rows
TODO: Change the format of the data to be in the format of SNV. For now I am using Missense Variants

All of the mutations are INDELs, so they are all insertions or deletions of bases in the genome.
This means that they are not in the format of SNV, so how do I use drivr base to get the labels for the dataset?

CHROMOSOME, GENOME_START, GENOME_STOP, GENOMIC_WT_ALLELE, GENOMIC_MUT_ALLELE
chr1, 1234, 1235, ATAGTA, NaN

We need to change this to the format of SNV:
CHROMOSOME, POSITION, REF, ALT
chr1, 1234, A, T


'''


#Now Define a function which does all of this to a file:

def get_labels(file_path, output_path):
    cosmic = pd.read_csv(file_path, delimiter='\t', low_memory=False)
    # Extract random 100000 rows from the dataset for testing
    df = pd.DataFrame(cosmic, [i for i in random.sample(range(len(cosmic)), 10000)]) 
    df = df[~df['CHROMOSOME'].isin(['X','Y','MT'])] #Disregard rows with chromosome = X, Y 
    filtered_df = df[df['MUTATION_DESCRIPTION'].isin(['missense_variant'])]

    #Genomic_WT_ALLELE = WT allele, Genomic_MUT_ALLELE = Mutant allele resulting from the mutation
    
    features = filtered_df[['CHROMOSOME','GENOME_START','GENOME_STOP','GENOMIC_WT_ALLELE','GENOMIC_MUT_ALLELE']]
    features_copy = features.copy()
    
    #GENOME_START AND STOP NEED TO BE TYPE INT
    features_copy['GENOME_START'] = features_copy['GENOME_START'].astype(int)
    features_copy['GENOME_STOP'] = features_copy['GENOME_STOP'].astype(int)
    #Fill in blank spaces with NaN - doesn't work for some reason
    features_copy = features_copy.replace('', np.nan) 
    
    
    #Add driver status column
    # features_copy['DRIVER'] = np.ones(len(features_copy))

    #Change the feature names:
    chromosome = features.loc[:,'CHROMOSOME']
    for chr in chromosome:
        chromosome = chromosome.replace(chr, 'chr'+str(chr))

    features_copy['CHROMOSOME'] = chromosome
    features_copy.sort_values(by=['GENOME_START'], ascending=True, inplace=True)

    #Remove repeated rows
    features_copy_no_repeats = features_copy.drop_duplicates(keep='first')
    
    features_copy_no_repeats.to_csv(output_path, sep='\t', index=False, header=False)

data_path = 'Cosmic/Cosmic_MutantCensus_Tsv_v100_GRCh38/Cosmic_MutantCensus_v100_GRCh38.tsv'

output_path = 'Cosmic/Cosmic_data/missense.tsv'

# data_example = 'Cosmic/Cosmic_INDEL_100/Cosmic_MutantCensus_v100_GRCh38 copy.tsv'
# example_output_path = 'Cosmic/Cosmic_INDEL_100/new.bed'
get_labels(data_path, output_path)

'''

To illustrate how single–nucleotide features are integrated into an indel score, 
let us consider a three– nucleotide reference sequence AGC. 
A two–nucleotide insertion TT starting at the second position 
will yield the mutant sequence ATTGC. 
The mutation’s overall conservation score includes the 
score for the first nucleotide, A, along with the scores from the inserted nucleotides TT. 
We compute the effect of a transition from the three–nucleotide 
wildtype sequence AGC to the three–nucleotide mutant sequence ATT by averaging 
the conversion scores for three transitions: A → A, G → T, and C → T. 
The approach for deletions is similar: we simply average the scores for mutations along 
the length of a deletion.

Insertions in the cosmic data are rows with a NaN value in the Wild-Type allele column
Deletions in the cosmic data are rows with a NaN value in the Mutant allele column - and delete the reference allele

I need to find the adjacent Alleles to the indel and then calculate the conservation score for the indel
How to find adjacent alleles? - Use the reference genome to find the adjacent alleles

Insertions positions are represented as a space of length 1 between the positions.


Deletions where l, r are the left and right adjacent alleles to the deletion:

Deletion:
chr8    1903380     1903380     T         NaN

If length of the deletion is 1, keep the SNV as it is, else separate into discrete SNVs:

chr11	1020723	1020727	TGGTG  NaN

becomes: 

chr11	1020723	1020723	T	NaN
chr11	1020724	1020724	G	NaN
chr11	1020725	1020725	G	NaN
chr11	1020726	1020726	T	NaN
chr11	1020727	1020727	G	NaN


Insertion: 
chr20	1915303	    1915304	    NaN	            GT


becomes: 
chr20   1915303     1915304     ref (C)         ref (C) C in this case
chr20   1915304     1915305     r1_adj (C)     G 
chr20   1915305     1915306     r2_adj (T)     T
chr20   1915306     1915307     C              C     
chr20   1915307     1915308     T              T



'''


def fetch_ucsc_variants_insertions(indel): 
    '''Gets the first 2 right adjacent alleles to the indel usng the UCSC API'''

    chromosome, start, end, wt, mut = indel[0], indel[1], indel[2], indel[3], indel[4]
    
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chromosome};start={start};end={end+len(mut)}"
    
    response = requests.get(url)
    
    results = response.json()

    
    return results['dna'].upper()



def snv_parametrisation_of_indel(indel):
    '''Converts the indel to discrete SNVs'''
    #Indel = [chromosome, start, end, ref, alt]
    #Get the adjacent alleles to the indel
    
    #Insertion
    if indel[3] == 'NaN':
        dna = fetch_ucsc_variants_insertions(indel) #retrieve the right adjacent alleles
        alt = indel[4]
        #Initialise the dataframe with the insertion position
        df = pd.DataFrame({'CHROMOSOME': indel[0], 'GENOME_START': indel[1],'GENOME_END':indel[1]+1, 'GENOMIC_WT_ALLELE': dna[0], 'GENOMIC_MUT_ALLELE': dna[0]}, index=[0])

        trimmed_dna = dna[1:]
        for i in range(len(alt)):
            df = pd.concat([df, pd.DataFrame({'CHROMOSOME': indel[0], 'GENOME_START': indel[1]+i+1,'GENOME_END':indel[1]+i+2, 'GENOMIC_WT_ALLELE': trimmed_dna[i], 'GENOMIC_MUT_ALLELE': alt[i]}, index=[0])])
        return df
    #Deletion: 
    elif indel[4] == 'NaN':
        ref = indel[3]
        #Initialise the dataframe with the deletion position
        df = pd.DataFrame()
        trimmed_ref = ref[1:] 
        for i in range(len(ref)):
            df = pd.concat([df, pd.DataFrame({'CHROMOSOME': indel[0], 'GENOME_START': indel[1]+i,'GENOME_END':indel[1]+i+1, 'GENOMIC_WT_ALLELE': ref[i], 'GENOMIC_MUT_ALLELE': np.nan}, index=[0])])
        return df

#Need to reformat the location of the indel insertion so that it covers the range of the insertion
#At the moment it is just the location of the insertion

# indel = ['chr11',1020723, 1020727, 'NaN', 'TGGTG']
# indel = ['chr20',1915303, 1915304, 'NaN', 'GTGT']
# df = snv_parametrisation_of_indel(indel)
# print(df)


def snv_converted_data_set(file_path, output_path):
    """Convert the INDELS from get_labels function to a dataset of discrete SNVs"""
    frameshift_data = pd.read_csv(file_path, delimiter='\t',low_memory=False, header=None)
    frameshift_data = frameshift_data.replace(np.nan, 'NaN', regex=True)
    df = pd.DataFrame(frameshift_data, [i for i in range(0,100)]) #if frameshift_data.iloc[i,3] == 'NaN' and len(frameshift_data.iloc[i,3]) > 1])

    snv = pd.DataFrame(columns=['CHROMOSOME', 'GENOME_START', 'GENOME_END', 'GENOMIC_WT_ALLELE', 'GENOMIC_MUT_ALLELE'])
    # count = 0
    for i,row in df.iterrows():
        row_list = list(row)
        #Just checking the number of single snvs match with the length of the indel
        # length = len(row_list[4]) + 1
        # count+=length 
        snv_df = snv_parametrisation_of_indel(row_list)
        snv = pd.concat([snv, snv_df])

    snv.to_csv(output_path, sep='\t', index=False, header=False)
    

input_file = 'Cosmic/Cosmic_data/frameshift.tsv'
output = 'Cosmic/Cosmic_data/frameshift_snv.tsv'
# snv_converted_data_set(input_file, output)


'''

Convert the discrete SNVs to a bed file format, where it is in the format compatible with DrivR-Base:
[chromosome, position, ref, alt]

Could do this concurrently with the snv_converted_data_set function, but for now I will do it separately

'''

def snv_to_bed(file_path, output_path):
    '''
    Converts the discrete SNVs to a bed file format, in the format compatible with DrivR-Base
    Also adds a driver stat column which are all 1's becuase all of the mutations are drivers
    '''
    snv = pd.read_csv(file_path, delimiter='\t', low_memory=False, header=None)
    snv = snv.replace(np.nan, '-', regex=True)
    snv = snv.drop_duplicates(keep='first')
    #drop the 3rd column
    snv = snv.drop(columns=[2])
    snv[5] = np.ones(len(snv))
    snv.to_csv(output_path, sep='\t', index=False, header=False)

input_file = 'Cosmic/Cosmic_data/frameshift_snv.tsv'
output = 'Cosmic/Cosmic_data/frameshift_snv.bed'
# snv_to_bed(input_file, output)

def missense_to_bed(file_path, output_path):
    '''
    Converts the missense mutations to a bed file format, in the format compatible with DrivR-Base
    Also adds a driver stat column which are all 1's becuase all of the mutations are drivers
    '''
    missense = pd.read_csv(file_path, delimiter='\t', low_memory=False, header=None)
    missense = missense.drop_duplicates(keep='first')
    missense = missense.drop(columns=[2])
    missense[5] = np.ones(len(missense))
    #remove rows where the alt or ref allele is > length 1:
    missense = missense[missense[3].map(len) == 1]
    missense = missense[missense[4].map(len) == 1]
    missense.to_csv(output_path, sep='\t', index=False, header=False)
input_file = 'Cosmic/Cosmic_data/missense.tsv'
output = 'Cosmic/Cosmic_data/missense.bed'
missense_to_bed(input_file, output)



'''
TODO: Add command line arguments to the functions so that they can be run from the command line


TODO: Test the newly created dataset (frameshift_snv.bed) with DrivR-Base and see if it works... 
Testing: 
    - The variants_with_driver_stat.bed file is the file that is used to test the DrivR-Base tool.
    - This then changes to variant.bed, which for some reason moves driver status values to the empty NaN values in the 
        genomic_mut_allele column. Not sure why.
    - Not sure what to do about this? Maybe instead of using NaN values use something else? 
    - Look into the vcf2bed tools on how to convert Indels to SNVs?
    - FATHMM-INDEL creates scores for the INDELs and doesn't convert them to SNVs, I need to convert to SNVs 
        in order to use DrivR-Base



- Will need to check with Colin to see if what I've done is a valid approach.
    - For the insertions, I could adapt the same approach as deletions, if the NaN values work with DrivR-Base. But I'm not sure if they do. 
    - So far I've just followed the FATHMM-INDEL approach 
- If it does, then we need to collate some neutral data to balance the dataset
- If it doesn't, then we need to find out why and fix it
- Then upon testing the neutral data with DrivR-Base, we can start to build an ML model.

'''

# if __name__ == '__main__':
#     bigdata_path = sys.argv[1]
#     frameshift_path = sys.argv[2]
#     snv_path = sys.argv[3]
#     get_labels(bigdata_path, frameshift_path) #Get the frameshift data
#     snv_converted_data_set(frameshift_path, snv_path) #Convert the frameshift data to SNVs
#     snv_to_bed(snv_path, snv_path.replace('.tsv', '.bed')) #Convert the SNVs to bed file format




