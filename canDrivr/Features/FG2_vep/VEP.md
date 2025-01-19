
Plan for VEP annotations;
  - Don't delete the original files which are made in the ./get_vep.sh and use them to test the python scripts and make sure they are not specific to
    SNVs, i.e. they work for indels and it all works. I think it does, but further investigation is needed.
  - Add in some predictions from down below which give good indications on deleteriousness of indels / difference between simulated and human indels like in CADD
  -  


Other prediction tools (plugins) for vep:
# Pathogenicity predictions
   1: dbNSFP                   - dbNSFP provides pathogenicity predictions for missense variants from various algorithms
   2: AlphaMissense            - Annotates missense variants with the pre-computed AlphaMissense pathogenicity scores. AlphaMissense is a deep learning model developed by Google DeepMind that predicts the pathogenicity of single nucleotide missense variants.
   3: CADD                     - Combined Annotation Dependent Depletion (CADD) is a tool for scoring the deleteriousness of single nucleotide variants and insertion/deletion variants in the human genome. CADD integrates multiple annotations into one metric by contrasting variants that survived natural selection with simulated mutations. CADD is only available here for non-commercial use. See CADD website for more information.
   4: CAPICE                   - Consequence-Agnostic Pathogenicity Interpretation of Clinical Exome variations (CAPICE) is a tool for scoring the pathogenicity of single nucleotide variants and insertion/deletion variants in the human genome. CAPICE uses a gradient boosting tree model trained with multiple genomic annotations used by CADD score and based on clinical significance.
   5: FATHMM                   - 
   6: FATHMM_MKL               - FATHMM-MKL predicts functional consequences of variants, both coding and non-coding.
   7: Gwava                    - Retrieves precomputed Genome Wide Annotation of VAriants (GWAVA) scores for any  variant that overlaps a known variant from the Ensembl variation database
   8: Carol                    - Calculates the Combined Annotation scoRing toOL (CAROL) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
   9: Condel                   - Calculates the Consensus Deleteriousness (Condel) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
  10: LoF                      - LOFTEE identifies LoF (loss-of-function) variation
  11: LoFtool                  - Provides a per-gene rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations in ExAC data
  12: MPC                      - MPC is a missense deleteriousness metric based on the analysis of genic regions depleted of missense mutations in ExAC
  13: MTR                      - MTR scores quantify the amount of purifying selection acting specifically on missense variants in a given window of protein-coding sequence
  14: PrimateAI                - PrimateAI scores can be used to assess pathogenicity of variants in humans
  15: REVEL                    - An ensemble method for predicting the pathogenicity of rare missense variants
  16: PON_P2                   - 
  17: ClinPred                 - A prediction tool for the identification of disease-relevant nonsynonymous single nucleotide variants
  18: EVE                      - EVE (evolutionary model of variant effect) is a computational method for the classification of human genetic variants trained solely on evolutionary sequences.
  19: BayesDel                 - BayesDel adds a deleteriousness meta-score combining multiple deleteriousness predictors.
  20: PolyPhen_SIFT            - 
  21: VARITY                   - VARITY is a tool to predict pathogenicity of rare human missense variants. VARITY uses a gradient boosting tree model trained with rare variants from multiple genomic annotations and based on clinical significance.
# Phenotype data and citations
  22: AVADA                    - 
  23: Phenotypes               - Retrieves overlapping phenotype annotations
  24: GO                       - Retrieves Gene Ontology terms associated with transcripts/translations
  25: G2P                      - Assesses variants using G2P allelic requirements for potential phenotype involvement.
  26: Geno2MP                  - Geno2MP is a web-accessible database of rare variant genotypes linked to individual-level phenotypic profiles defined by human phenotype ontology (HPO) terms
  27: PostGAP                  - 
  28: satMutMPRA               - Reports saturation mutagenesis MPRA measures of variant effects on gene RNA expression for 21 enhancers/promoters
  29: DisGeNET                 - Variant-Disease-PMID associations from the DisGeNET database (https://www.disgenet.org/)
  30: Mastermind               - Variants with clinical evidence cited in the medical literature reported by Mastermind Genomic Search Engine (https://www.genomenon.com/mastermind)
  31: GWAS                     - Genome-wide associated study data from the NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/)
  32: PhenotypeOrthologous     - Phenotypes associated with orthologous genes in model organisms (mouse and rat)
# Pathogenicity predictions
  33: dbNSFP                   - dbNSFP provides pathogenicity predictions for missense variants from various algorithms
  34: AlphaMissense            - Annotates missense variants with the pre-computed AlphaMissense pathogenicity scores. AlphaMissense is a deep learning model developed by Google DeepMind that predicts the pathogenicity of single nucleotide missense variants.
  35: CADD                     - Combined Annotation Dependent Depletion (CADD) is a tool for scoring the deleteriousness of single nucleotide variants and insertion/deletion variants in the human genome. CADD integrates multiple annotations into one metric by contrasting variants that survived natural selection with simulated mutations. CADD is only available here for non-commercial use. See CADD website for more information.
  36: CAPICE                   - Consequence-Agnostic Pathogenicity Interpretation of Clinical Exome variations (CAPICE) is a tool for scoring the pathogenicity of single nucleotide variants and insertion/deletion variants in the human genome. CAPICE uses a gradient boosting tree model trained with multiple genomic annotations used by CADD score and based on clinical significance.
  37: FATHMM                   - 
  38: FATHMM_MKL               - FATHMM-MKL predicts functional consequences of variants, both coding and non-coding.
  39: Gwava                    - Retrieves precomputed Genome Wide Annotation of VAriants (GWAVA) scores for any  variant that overlaps a known variant from the Ensembl variation database
  40: Carol                    - Calculates the Combined Annotation scoRing toOL (CAROL) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
  41: Condel                   - Calculates the Consensus Deleteriousness (Condel) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
  42: LoF                      - LOFTEE identifies LoF (loss-of-function) variation
  43: LoFtool                  - Provides a per-gene rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations in ExAC data
  44: MPC                      - MPC is a missense deleteriousness metric based on the analysis of genic regions depleted of missense mutations in ExAC
  45: MTR                      - MTR scores quantify the amount of purifying selection acting specifically on missense variants in a given window of protein-coding sequence
  46: PrimateAI                - PrimateAI scores can be used to assess pathogenicity of variants in humans
  47: REVEL                    - An ensemble method for predicting the pathogenicity of rare missense variants
  48: PON_P2                   - 
  49: ClinPred                 - A prediction tool for the identification of disease-relevant nonsynonymous single nucleotide variants
  50: EVE                      - EVE (evolutionary model of variant effect) is a computational method for the classification of human genetic variants trained solely on evolutionary sequences.
  51: BayesDel                 - BayesDel adds a deleteriousness meta-score combining multiple deleteriousness predictors.
  52: PolyPhen_SIFT            - 
  53: VARITY                   - VARITY is a tool to predict pathogenicity of rare human missense variants. VARITY uses a gradient boosting tree model trained with rare variants from multiple genomic annotations and based on clinical significance.
# Variant data
  54: Paralogues               - Retrieves ClinVar variants that overlap genomic coordinates corresponding to aligned amino acid positions in paralogous proteins
  55: DeNovo                   - Identifies de novo variants in a VCF file
  56: OpenTargets              - Returns locus-to-gene (L2G) scores to predict causal genes at GWAS loci from Open Targets Genetics
  57: LD                       - Finds variants in linkage disequilibrium with any overlapping existing variants from the Ensembl variation databases
  58: SameCodon                - Reports existing variants that fall in the same codon
  59: LOVD                     - Retrieves LOVD variation data
  60: SubsetVCF                - 
  61: GeneBe                   - Retrieves automatic ACMG variant classification data from GeneBe
# Pathogenicity predictions
  62: dbNSFP                   - dbNSFP provides pathogenicity predictions for missense variants from various algorithms
  63: AlphaMissense            - Annotates missense variants with the pre-computed AlphaMissense pathogenicity scores. AlphaMissense is a deep learning model developed by Google DeepMind that predicts the pathogenicity of single nucleotide missense variants.
  64: CADD                     - Combined Annotation Dependent Depletion (CADD) is a tool for scoring the deleteriousness of single nucleotide variants and insertion/deletion variants in the human genome. CADD integrates multiple annotations into one metric by contrasting variants that survived natural selection with simulated mutations. CADD is only available here for non-commercial use. See CADD website for more information.
  65: CAPICE                   - Consequence-Agnostic Pathogenicity Interpretation of Clinical Exome variations (CAPICE) is a tool for scoring the pathogenicity of single nucleotide variants and insertion/deletion variants in the human genome. CAPICE uses a gradient boosting tree model trained with multiple genomic annotations used by CADD score and based on clinical significance.
  66: FATHMM                   - 
  67: FATHMM_MKL               - FATHMM-MKL predicts functional consequences of variants, both coding and non-coding.
  68: Gwava                    - Retrieves precomputed Genome Wide Annotation of VAriants (GWAVA) scores for any  variant that overlaps a known variant from the Ensembl variation database
  69: Carol                    - Calculates the Combined Annotation scoRing toOL (CAROL) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
  70: Condel                   - Calculates the Consensus Deleteriousness (Condel) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
  71: LoF                      - LOFTEE identifies LoF (loss-of-function) variation
  72: LoFtool                  - Provides a per-gene rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations in ExAC data
  73: MPC                      - MPC is a missense deleteriousness metric based on the analysis of genic regions depleted of missense mutations in ExAC
  74: MTR                      - MTR scores quantify the amount of purifying selection acting specifically on missense variants in a given window of protein-coding sequence
  75: PrimateAI                - PrimateAI scores can be used to assess pathogenicity of variants in humans
  76: REVEL                    - An ensemble method for predicting the pathogenicity of rare missense variants
  77: PON_P2                   - 
  78: ClinPred                 - A prediction tool for the identification of disease-relevant nonsynonymous single nucleotide variants
  79: EVE                      - EVE (evolutionary model of variant effect) is a computational method for the classification of human genetic variants trained solely on evolutionary sequences.
  80: BayesDel                 - BayesDel adds a deleteriousness meta-score combining multiple deleteriousness predictors.
  81: PolyPhen_SIFT            - 
  82: VARITY                   - VARITY is a tool to predict pathogenicity of rare human missense variants. VARITY uses a gradient boosting tree model trained with rare variants from multiple genomic annotations and based on clinical significance.
# Splicing predictions
  83: dbscSNV                  - Retrieves data for splicing variants from a tabix-indexed dbscSNV file
  84: GeneSplicer              - Detects splice sites in genomic DNA
  85: MaxEntScan               - Sequence motif and maximum entropy based splice site consensus predictions
  86: SpliceRegion             - More granular predictions of splicing effects
  87: SpliceAI                 - Pre-calculated annotations from SpliceAI a deep neural network, developed by Illumina, Inc that predicts splice junctions from an arbitrary pre-mRNA transcript sequence. Used for non-commercial purposes (https://github.com/Illumina/SpliceAI)
  88: SpliceVault              - Retrieves the most common variant-associated mis-splicing events for variants that overlap a near-splice-site region
# Conservation
  89: Blosum62                 - BLOSUM62 amino acid conservation score
  90: Conservation             - Retrieves a conservation score from the Ensembl Compara databases for variant positions
  91: AncestralAllele          - Retrieves the ancestral allele for variants inferred from the Ensembl Compara Enredo-Pecan-Ortheus (EPO) pipeline
# Frequency data
  92: gnomADc                  - Reports coverage data from the gmomAD coverage files
# Phenotype data and citations
  93: AVADA                    - 
  94: Phenotypes               - Retrieves overlapping phenotype annotations
  95: GO                       - Retrieves Gene Ontology terms associated with transcripts/translations
  96: G2P                      - Assesses variants using G2P allelic requirements for potential phenotype involvement.
  97: Geno2MP                  - Geno2MP is a web-accessible database of rare variant genotypes linked to individual-level phenotypic profiles defined by human phenotype ontology (HPO) terms
  98: PostGAP                  - 
  99: satMutMPRA               - Reports saturation mutagenesis MPRA measures of variant effects on gene RNA expression for 21 enhancers/promoters
 100: DisGeNET                 - Variant-Disease-PMID associations from the DisGeNET database (https://www.disgenet.org/)
 101: Mastermind               - Variants with clinical evidence cited in the medical literature reported by Mastermind Genomic Search Engine (https://www.genomenon.com/mastermind)
 102: GWAS                     - Genome-wide associated study data from the NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/)
 103: PhenotypeOrthologous     - Phenotypes associated with orthologous genes in model organisms (mouse and rat)
# Variant data
 104: Paralogues               - Retrieves ClinVar variants that overlap genomic coordinates corresponding to aligned amino acid positions in paralogous proteins
 105: DeNovo                   - Identifies de novo variants in a VCF file
 106: OpenTargets              - Returns locus-to-gene (L2G) scores to predict causal genes at GWAS loci from Open Targets Genetics
 107: LD                       - Finds variants in linkage disequilibrium with any overlapping existing variants from the Ensembl variation databases
 108: SameCodon                - Reports existing variants that fall in the same codon
 109: LOVD                     - Retrieves LOVD variation data
 110: SubsetVCF                - 
 111: GeneBe                   - Retrieves automatic ACMG variant classification data from GeneBe
# Nearby features
 112: NearestGene              - Finds the nearest gene to non-genic variants
 113: NearestExonJB            - 
 114: Downstream               - Predicts the downstream effects of a frameshift variant on the protein sequence of a transcript
 115: TSSDistance              - Calculates the distance from the transcription start site for upstream variants   
# Sequence
 116: ProteinSeqs              - Prints out the reference and mutated protein sequences of any proteins found with non-synonymous mutations
 117: ReferenceQuality         - 
# Visualisation
 118: Draw                     - Creates images of the transcript model showing variant location
# Look up
 119: LocalID                  - Allows you to use variant IDs as VEP input without making a database connection.
# External ID
 120: FlagLRG                  - 
# Motif
 121: FunMotifs                - 
# HGVS
 122: HGVSIntronOffset         - 
 123: SingleLetterAA           - 
# Structural variant data
 124: StructuralVariantOverlap - 
# Gene tolerance to change
 125: DosageSensitivity        - DosageSensitivity determines the likelihood of a gene being haploinsufficient or triplosensitive
 126: LOEUF                    - LOEUF stands for the 'loss-of-function observed/expected upper bound fraction'. This plugin adds constraint scores derived from gnomAD to VEP
 127: pLI                      - Provides a per-gene or per transcript probability of being loss-of-function intolerant (pLI)
# Transcript annotation
 128: NMD                      - Nonsense-mediated mRNA decay escaping variants prediction
 129: UTRAnnotator             - A VEP plugin that annotates the effect of 5' UTR variant especially for variant creating/disrupting upstream ORFs.
 130: TranscriptAnnotator      - 
 131: RiboseqORFs              - A VEP plugin that calculates consequences for variants overlapping Ribo-seq ORFs
# Protein data
 132: neXtProt                 - Retrieves protein-related data from neXtProt
# Functional effect
 133: IntAct                   - IntAct provides molecular interaction data for variants as reported by IntAct database
 134: MaveDB                   - MaveDB holds experimentally determined measures of variant effect
# Pathogenicity predictions
 135: dbNSFP                   - dbNSFP provides pathogenicity predictions for missense variants from various algorithms
 136: AlphaMissense            - Annotates missense variants with the pre-computed AlphaMissense pathogenicity scores. AlphaMissense is a deep learning model developed by Google DeepMind that predicts the pathogenicity of single nucleotide missense variants.
 137: CADD                     - Combined Annotation Dependent Depletion (CADD) is a tool for scoring the deleteriousness of single nucleotide variants and insertion/deletion variants in the human genome. CADD integrates multiple annotations into one metric by contrasting variants that survived natural selection with simulated mutations. CADD is only available here for non-commercial use. See CADD website for more information.
 138: CAPICE                   - Consequence-Agnostic Pathogenicity Interpretation of Clinical Exome variations (CAPICE) is a tool for scoring the pathogenicity of single nucleotide variants and insertion/deletion variants in the human genome. CAPICE uses a gradient boosting tree model trained with multiple genomic annotations used by CADD score and based on clinical significance.
 139: FATHMM                   - 
 140: FATHMM_MKL               - FATHMM-MKL predicts functional consequences of variants, both coding and non-coding.
 141: Gwava                    - Retrieves precomputed Genome Wide Annotation of VAriants (GWAVA) scores for any  variant that overlaps a known variant from the Ensembl variation database
 142: Carol                    - Calculates the Combined Annotation scoRing toOL (CAROL) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
 143: Condel                   - Calculates the Consensus Deleteriousness (Condel) score for a missense mutation based on the pre-calculated SIFT and PolyPhen scores
 144: LoF                      - LOFTEE identifies LoF (loss-of-function) variation
 145: LoFtool                  - Provides a per-gene rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations in ExAC data
 146: MPC                      - MPC is a missense deleteriousness metric based on the analysis of genic regions depleted of missense mutations in ExAC
 147: MTR                      - MTR scores quantify the amount of purifying selection acting specifically on missense variants in a given window of protein-coding sequence
 148: PrimateAI                - PrimateAI scores can be used to assess pathogenicity of variants in humans
 149: REVEL                    - An ensemble method for predicting the pathogenicity of rare missense variants
 150: PON_P2                   - 
 151: ClinPred                 - A prediction tool for the identification of disease-relevant nonsynonymous single nucleotide variants
 152: EVE                      - EVE (evolutionary model of variant effect) is a computational method for the classification of human genetic variants trained solely on evolutionary sequences.
 153: BayesDel                 - BayesDel adds a deleteriousness meta-score combining multiple deleteriousness predictors.
 154: PolyPhen_SIFT            - 
 155: VARITY                   - VARITY is a tool to predict pathogenicity of rare human missense variants. VARITY uses a gradient boosting tree model trained with rare variants from multiple genomic annotations and based on clinical significance.
# Protein annotation
 156: mutfunc                  - mutfunc predicts destabilization effect of protein structure, interaction, regulatory region etc. caused by a variant
# Regulatory impact
 157: Enformer                 - Predictions of variant impact on gene expression
# Other plugins
 158: CSN    


