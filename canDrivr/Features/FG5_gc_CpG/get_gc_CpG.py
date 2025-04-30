from config import *

import os
import pandas as pd
from strkernel.mismatch_kernel import MismatchKernel
from strkernel.mismatch_kernel import preprocess
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import classification_report  # classification summary
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
import sys
from multiprocessing import Pool
from functools import reduce

if __name__ == "__main__":
    variantType = sys.argv[1]
    variants = sys.argv[1] + sys.argv[2]
    outputDir = sys.argv[3]

    # Reads in the human GRCh38 genome in fasta format
    record_dict = SeqIO.to_dict(SeqIO.parse("hg38.fa", "fasta"))

    # Reading in the variant file
    variants = pd.read_csv(variants, sep=",")
    variants.rename(columns={"start": "pos"}, inplace=True)
    
    # Drops variants on the sex chromosomes
    variants = variants[(variants['chrom'] != "chrX") & (variants['chrom'] != "chrY")]
    variants = variants.reset_index(drop=True)

    # Function to extract wild-type sequence and compute GC content and CpG count
    def getGCContent(variantIndex, windowSize):
        chrom = variants.loc[variantIndex, "chrom"]
        pos = int(variants.loc[variantIndex, "pos"])
        
        try:
            wildType = str(record_dict[chrom].seq[pos - 1 - windowSize:pos - 1 + windowSize]).upper()
        except KeyError:
            # If chromosome is not found in record_dict
            return 0, 0, 0
        
        if not wildType or len(wildType) == 0:
            return 0, 0, 0  # Handle empty sequences safely
        
        GCContent = round((wildType.count("G") + wildType.count("C")) / len(wildType) * 100, 2)
        CpGCount = wildType.count("CG")
        
        if wildType.count("G") != 0 and wildType.count("C") != 0:
            CpG_obs_exp = round((CpGCount * len(wildType)) / (wildType.count("G") * wildType.count("C")), 2)
        else:
            CpG_obs_exp = 0
        
        return GCContent, CpGCount, CpG_obs_exp

    # Applies the above function to all variants for a given window size
    def getGC(windowSize):
        GCDf = [getGCContent(x, windowSize) for x in range(len(variants))]
        GCdataframe = pd.concat([variants, pd.DataFrame(GCDf)], axis=1)
        GCdataframe["chrom"] = GCdataframe["chrom"].str.replace("chr", "").astype(str)
        GCdataframe["pos"] = GCdataframe["pos"].astype(int)
        GCdataframe = GCdataframe.rename(columns={
            0: f"{windowSize*2}GCContent",
            1: f"{windowSize*2}CpGCount",
            2: f"{windowSize*2}CpG_obs_exp"
        })
        
        if "stop" in GCdataframe.columns:
            GCdataframe = GCdataframe.drop("stop", axis=1)
        
        GCdataframe['chrom'] = 'chr' + GCdataframe['chrom'].astype(str)
        return GCdataframe

    # Applies function to different window sizes
    dataframeList = [getGC(i) for i in [10, 20, 30, 40, 50, 100, 250, 500, 1000]]

    # Merge DataFrames
    data_merge = reduce(lambda left, right: 
                         pd.merge(left, right,
                                  on=["chrom", "pos", "ref_allele", "alt_allele", "driver_stat"],
                                  how="outer"),
                         dataframeList)
    
    # Write the merged dataframe to a CSV file
    output_file = os.path.join(outputDir, "GC_indTest.csv")
    data_merge = data_merge.drop_duplicates(keep="first")
    data_merge.to_csv(output_file, index=False, sep=",")

