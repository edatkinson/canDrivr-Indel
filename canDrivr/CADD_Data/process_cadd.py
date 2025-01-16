
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

combined_data = merge_cadd('cadd_simulated_data.tsv', 'cadd_human_data.tsv')

combined_data.to_csv('combined_cadd_data.tsv', sep='\t', index=False)


