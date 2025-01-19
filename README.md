

Plan for the architecture of this project
---------------------------------------------------------------------------------------

- Get Data 

    - Filter indels from pathogenic dataset (cosmic)
        - `get_COSMIC.sh`
            - gets coding indels from cosmic from `data/Cosmic_GenomeScreensMutant_v100_GRCh38.tsv.gz`
            - Outputs `cosmic_coding_indels_short.tsv` and `cosmic_coding_indels_short_filtered.tsv`
        - data_get.py 
            - `extract_data()` Extracts the data for pathogenic indels from cosmic from `cosmic_coding_indels_short.tsv`
                - Count number of identical variants and add coloumn for the count
                - Results in `cosmic_filtered.tsv`, filtered for ~150k pathogenic indels
    - Filter for benign indels from clinvar
        - `get_clinvar.sh` Processes `variant_summary.txt` from clinvar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
            - Produces `clinvar_filtered`, filtered for ~13k benign indels.
    - Verify indels are within the coding region:
        - Hopefully reduces the size of the dataset.
    
    //??// Many more pathogenic indels than benign. How do I handle this/ how do I reduce the number of pathogenic indels. Do I select the indels from cosmic with a count above a certain threshold? 

    - TODO;
        - Assign driver status.
        - Combine datasets.
        - Save final dataset.


New Objective:

- Obtain Cadd Training data: from https://krishna.gs.washington.edu/download/CADD-development/v1.4/training_data/GRCh38/ 
    - Done: `cadd_human_data.tsv` and `cadd_simulated_data.tsv`. These are derived from `humanDerived_InDels.tsv` and `simulation_InDels.tsv` respectively.
    - Need to:
        - Randomly select a set of pathogenic indels with equal chromosomal group sizes to match the lower number of neutral coding indels.
        - Combine Dataset: Done in `process_cadd.py` and outputs `combined_cadd_data.tsv`
        - Annotate

---------------------------------------------------------------------------------------
- Feature Allocation
    
    - Conservation Scores:
        - Batch Queries to UCSC genome browser to obtain cons scores for each indel.
        - Takes a very, very long time. Try in UNI, or find a way to do it locally.
        - How? 
            - Run `get_conservation.py` - this uses `final_cadd_data.bed` which is obtained from running `process_cadd.py` with `combined_cadd_data.tsv`
        - Output?
            - In files/ will appear all the cons files separately
        - Needed? Merging of all features, will do this at the end when I have all of the features
    
    - VEP Scores:
        - Outputs one-hot coded scores for each indel in `final_cadd_data.bed` in /output/
        - In order to work, you must have vep downloaded.
            - clone ensemble vep, perl install.pl, then download the cache from ensemble documentation.
            - For now it is located in `/home/colin./.vep/`.
            - Python files from DrivR-Base and ./get_vep.sh do the one-hot coding of the output from vep.
            - I have run `./get_vep.sh /home/colin/canDrivr-Indel/canDrivr/CADD_Data/ final_cadd_data.bed /home/colin/canDrivr-Indel/canDrivr/Features/FG2_vep/output/ ` which is (variant directory) (filename) (output directory)
        
        - Look into using values from other predictors in `VEP.md` for extra features - opportunity to expand report.
    
    Next:

    - Dinucleotide Properties:

    - DNA Shape: - half done.

    - gc_CpG:

    - Kernal:

    - Amino Acid Substitution Matrices:

    - Amino Acid Properties

    - Encode

    - Alpha Fold

    When finished with annotation:
        - Create a pipeline which does it all automatically so I can be more efficient with annotating my data.
        - Provide options to omit some features (like the slow ones) to streamline testing and modelling.

---------------------------------------------------------------------------------------

- Modelling 

    - Comparative Analysis:

        - CADD

        - Fathmm-INDEl

        - Pred-CID

        - Sift
    
    - xgBoost

    - SVM

    - other models


---------------------------------------------------------------------------------------


