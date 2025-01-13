

Plan for the new architecture of this project
---------------------------------------------------------------------------------------

- Get Data 

    - Filter indels from pathogenic dataset (cosmic)
        - get_COSMIC.sh
            - gets coding indels from cosmic
        - data_get.py 
            - `extract_data` Extracts the data
            - Count number of identical variants and add coloumn for the count
        - TODO;
            - Repeat for Neutral Data and combine the functionality so it works automatically
        - TODO;
            - Assign driver status.
            - Combine datasets.
            - Save final dataset.

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



##gff-version 3
##gvf-version 1.07
##file-date 2024-08-30
##genome-build ensembl GRCh38
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.283
##data-source Source=ensembl;version=113;url=http://vertebrates.ensembl.org/Homo_sapiens
##file-version 113
##sequence-region 19 1 58617616
##sequence-region 21 1 46709983
##sequence-region 13 1 114364328
##sequence-region 6 1 170805979
##sequence-region 1 1 248956422
##sequence-region MT 1 16569
##sequence-region Y 1 57227415
##sequence-region 4 1 190214555
##sequence-region 2 1 242193529
##sequence-region 12 1 133275309
##sequence-region 22 1 50818468
##sequence-region 10 1 133797422
##sequence-region 8 1 145138636
##sequence-region 15 1 101991189
##sequence-region X 1 156040895
##sequence-region 7 1 159345973
##sequence-region 16 1 90338345
##sequence-region 17 1 83257441
##sequence-region 3 1 198295559
##sequence-region 5 1 181538259
##sequence-region 11 1 135086622
##sequence-region 18 1 80373285
##sequence-region 14 1 107043718
##sequence-region 9 1 138394717
##sequence-region 20 1 64444167