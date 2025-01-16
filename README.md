

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
        - Will annotate features to the whole dataset once it is done.
            - During testing use a diluted dataset, i expect the whole annotation will take a long time for 100,000+ variants.
    - DNA_shape:
        - This should hopefully work easily like before with integration using a shell script.
        - Look out for NaN bias here.
    - New features:
        - Once I have my combined dataset (TODOs are done), I'll start using the linux machine to add more features. 
        - May be able to copy and paste the drivr base code and install the libraries locally which are used instead of setting up the venv listed on the github
    - Feature combination:
        - Should features be individually found, then combined to a final big file at the end?
        - Explore options for this.
    
    - PubMed Features

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


