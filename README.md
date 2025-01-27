

# Plan for the architecture of this project
---------------------------------------------------------------------------------------

<!-- - Get Data 

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
        - Save final dataset. -->


## New Objective:

- Obtain Cadd Training data: from https://krishna.gs.washington.edu/download/CADD-development/v1.4/training_data/GRCh38/ 
    - Done: `cadd_human_data.tsv` and `cadd_simulated_data.tsv`. These are derived from `humanDerived_InDels.tsv` and `simulation_InDels.tsv` respectively.
    - Need to:
        - Randomly select a set of pathogenic indels with equal chromosomal group sizes to match the lower number of neutral coding indels.
        - Combine Dataset: Done in `process_cadd.py` and outputs `combined_cadd_data.tsv`
        - Annotate

---------------------------------------------------------------------------------------
### Feature Allocation
    
    - Conservation Scores:
        - Done on Linux although can be also done on mac.
        - Batch Queries to UCSC genome browser to obtain cons scores for each indel.
        - Takes a very, very long time. Try in UNI, or find a way to do it locally.
        - How? 
            - Run `get_conservation.py` - this uses `final_cadd_data.bed` which is obtained from running `process_cadd.py` with `combined_cadd_data.tsv`
        - Output?
            - In output/ will appear all the cons files separately
        - Remember to sort the individual files by position, so they match up with other features.
    
    - VEP Scores:
        - Done on Linux
        - Outputs one-hot coded scores for each indel in `final_cadd_data.bed` in /output/
        - In order to work, you must have vep downloaded.
            - clone ensemble vep, perl install.pl, then download the cache from ensemble documentation.
            - For now it is located in `/home/colin./.vep/`.
            - Python files from DrivR-Base and ./get_vep.sh do the one-hot coding of the output from vep.
            - I have run `./get_vep.sh /home/colin/canDrivr-Indel/canDrivr/CADD_Data/ final_cadd_data.bed /home/colin/canDrivr-Indel/canDrivr/Features/FG2_vep/output/ ` which is (variant directory) (filename) (output directory)
        
        - Look into using values from other predictors in `VEP.md` for extra features - opportunity to expand report.
    

    - Dinucleotide Properties:
        - Done on mac.
        - Output located in features/FG3_dinucleotide_properites/output
        - Need to check if it is producing the correct results.

    - DNA Shape:
        - Done on mac.
        - Output located in Features/DNA_shape/output
        - Need to check if it is producing the correct results.

    - gc_CpG:
        - Done on mac.
        - Output located in Features/FG5_gc_CpG/output
        - Need to check if it is producing the correct results.

    - Kernal:
        - TODO
        - Need to make changes to the logic to include the following:
            - Handle Indel lengths: sequence extraction logic needs to account for the length of the indel.
            - Adjust sequence extraction for insertions and deletions.
                - For insertions:
                    - Insert the alternate allele into the reference sequence
                - Deletions:
                    - Remove the reference allele from the sequence.

    - Amino Acid Substitution Matrices:
        - TODO

    - Amino Acid Properties:
        - TODO

    - Encode:
        - TODO

    - Alpha Fold
        - TODO


    When finished with annotation:
        - Create a pipeline which does it all automatically so I can be more efficient with annotating my data.
        - Gonna be hard because I've done it on both mac and linux
        - Provide options to omit some features (like the slow ones) to streamline testing and modelling.


Merged Features in `merge_featues.py` which produces `Annotated_data.csv` which may need a reduction in some features as it has 1700+ columns. But this can be done in preprocessing for the ML models. 

---------------------------------------------------------------------------------------

### Plan for `get_features.sh`:
- user inputs variant file with driver stat
- create file with no driver stat
- make output folder to store all features.
- run scripts.
- merge outputs.
- delete prev feature files to free space.
- 



---------------------------------------------------------------------------------------

### Modelling 

    Will start off by playing around in a jupyter notebook as this is more efficient for testing in long processing time scripts.
    
    - `xgboost_model.ipynb` is where I will be first implementing models with the set of features I've found so far (DNA_shape, VEP, gc CpG, dinucleotide properties)
        - Still need to add conservation scores which I will try to do in uni with faster internet connection. This will probably produce best indicators and so are very important.
    
    - Tested a model using the above features. Get quite a low average CV score of: 73%, but get an overall test accuracy of 83%. This is without the conservation features, so hopefully they can help out with the accuracy by a few percentage points. 

    - Research methods in ML which can boost performance which are relevant to what I am doing, try something which my supervisor hasn't explicitly mentioned, such as LOCO-CV etc.

    - I've added the conservation features, bumps up accuracy to ~85-86%. Feature reduction improves accuracy by 1%.
    
    - Plan for improving model:
        - Feature importance selection, include with CV and find the best features on average and select those for final model.
        - I don't think LOCO CV is working atm, look into other applications of it and do that.
        - Hyperparameter tuning using GridSearch or other optimisation methods.
        - Class imbalance resampling, as I have more of some chromosomes than others.
        - Look into permutation importance and decide if it is worth implementing compared to feature importance.
        - 


---------------------------------------------------------------------------------------
### Comparative Analysis

    - The real test comes on an independent test set and seeing how it compares to the models below. 

    - Comparative Analysis:

        - CADD

        - Fathmm-INDEl

        - Pred-CID

        - Sift
    
    - xgBoost

    - SVM

    - other models


---------------------------------------------------------------------------------------