import pandas as pd
from collections import Counter
import numpy as np
file = '/Users/edatkinson/Repos/Modelling/new_cosmic.tsv'
cons_features = ['phyloP20way','k36.Bismap.MultiTrackMappability','k50.Bismap.MultiTrackMappability',
                 'phyloP30way','phyloP17way','phyloP100way','phastCons470way','k24.Umap.MultiTrackMappability',
                 'phyloP7way','phastCons20way','k100.Umap.MultiTrackMappability','phyloP470way','k50.Umap.MultiTrackMappability',
                 'phastCons100way','phastCons17way','phastCons30way','k36.Umap.MultiTrackMappability','phastCons7way',
                 'k24.Bismap.MultiTrackMappability','k100.Bismap.MultiTrackMappability','phyloP4way','phastCons4way'
            ]

dna_shape = [
    '3_MGW', '4_MGW', '5_MGW', '6_MGW', '7_MGW', '8_MGW', '9_MGW', '10_MGW', '11_MGW', '12_MGW', 
    '13_MGW', '14_MGW', '15_MGW', '16_MGW', '17_MGW', '18_MGW', '2_HelT', '3_HelT', '4_HelT', '5_HelT', 
    '6_HelT', '7_HelT', '8_HelT', '9_HelT', '10_HelT', '11_HelT', '12_HelT', '13_HelT', '14_HelT', 
    '15_HelT', '16_HelT', '17_HelT', '18_HelT', '3_ProT', '4_ProT', '5_ProT', '6_ProT', '7_ProT', 
    '8_ProT', '9_ProT', '10_ProT', '11_ProT', '12_ProT', '13_ProT', '14_ProT', '15_ProT', '16_ProT', 
    '17_ProT', '18_ProT', '2_Roll_x', '3_Roll_x', '4_Roll_x', '5_Roll', '6_Roll', '7_Roll', '8_Roll', 
    '9_Roll', '10_Roll', '11_Roll', '12_Roll', '13_Roll', '14_Roll', '15_Roll', '16_Roll', '17_Roll', 
    '18_Roll', '3_EP', '4_EP', '5_EP', '6_EP', '7_EP', '8_EP', '9_EP', '10_EP', '11_EP', '12_EP', 
    '13_EP', '14_EP', '15_EP', '16_EP', '17_EP', '18_EP'
]


alphaFold = ['X_coordinate','Y_coordinate','Z_coordinate','atom_site_occupancy','isotropic_temperature']

gcContent = [
    '20GCContent', '20CpGCount', '20CpG_obs_exp',
    '40GCContent', '40CpGCount', '40CpG_obs_exp',
    '60GCContent', '60CpGCount', '60CpG_obs_exp',
    '80GCContent', '80CpGCount', '80CpG_obs_exp',
    '100GCContent', '100CpGCount', '100CpG_obs_exp',
    '200GCContent', '200CpGCount', '200CpG_obs_exp',
    '500GCContent', '500CpGCount', '500CpG_obs_exp',
    '1000GCContent', '1000CpGCount', '1000CpG_obs_exp',
    '2000GCContent', '2000CpGCount', '2000CpG_obs_exp'
]





dna_props = [
    'WTtrinuc', 'mutTrinuc', '1_Bend', '1_Tip', '1_Inclination', '1_Major_Groove_Width', '1_Major_Groove_Depth',
    '1_Major_Groove_Size', '1_Major_Groove_Distance', '1_Minor_Groove_Width', '1_Minor_Groove_Depth',
    '1_Minor_Groove_Size', '1_Minor_Groove_Distance', '1_Persistance_Length',
    '1_Probability_contacting_nucleosome_core', '1_Mobility_to_bend_towards_major_groove',
    '1_Mobility_to_bend_towards_minor_groove', '1_Propeller_Twist', '1_Clash_Strength', '1_Shift_(RNA)',
    '1_Hydrophilicity_(RNA)', '1_Twist_(DNA-protein_complex)', '1_Twist_twist', '1_Tilt_tilt', '1_Roll_roll',
    '1_Twist_tilt', '1_Twist_roll', '1_Tilt_roll', '1_Shift_shift', '1_Slide_slide', '1_Rise_rise',
    '1_Shift_slide', '1_Shift_rise', '1_Slide_rise', '1_Twist_shift', '1_Twist_slide', '1_Twist_rise',
    '1_Tilt_shift', '1_Tilt_slide', '1_Tilt_rise', '1_Roll_shift', '1_Roll_slide', '1_Roll_rise',
    '1_Slide_stiffness', '1_Shift_stiffness', '1_Roll_stiffness', '1_Tilt_stiffness', '1_Twist_stiffness',
    '1_GC_content', '1_Purine_(AG)_content', '1_Keto_(GT)_content', '1_Adenine_content', '1_Guanine_content',
    '1_Cytosine_content', '1_Thymine_content', '1_Tilt_(DNA-protein_complex)', '1_Roll_(DNA-protein_complex)',
    '1_Shift_(DNA-protein_complex)', '1_Slide_(DNA-protein_complex)', '1_Rise_(DNA-protein_complex)', '1_Shift',
    '1_Slide', '1_Rise', '1_Wedge', '1_Direction', '1_Slide_(RNA)', '1_Rise_(RNA)', '1_Tilt_(RNA)',
    '1_Roll_(RNA)', '1_Twist_(RNA)', '1_Stacking_energy_(RNA)', '1_Rise_stiffness', '1_Melting_Temperature',
    '1_Stacking_energy', '1_Free_energy_(RNA)', '1_Enthalpy_(RNA)', '1_Entropy_(RNA)', '1_Tilt', '1_Roll',
    '1_Twist', '1_Flexibility_slide', '1_Flexibility_shift', '1_Enthalpy', '1_Entropy', '1_Free_energy',
    '2_Bend', '2_Tip', '2_Inclination', '2_Major_Groove_Width', '2_Major_Groove_Depth', '2_Major_Groove_Size',
    '2_Major_Groove_Distance', '2_Minor_Groove_Width', '2_Minor_Groove_Depth', '2_Minor_Groove_Size',
    '2_Minor_Groove_Distance', '2_Persistance_Length',
    '2_Probability_contacting_nucleosome_core', '2_Mobility_to_bend_towards_major_groove',
    '2_Mobility_to_bend_towards_minor_groove', '2_Propeller_Twist', '2_Clash_Strength', '2_Shift_(RNA)',
    '2_Hydrophilicity_(RNA)', '2_Twist_(DNA-protein_complex)', '2_Twist_twist', '2_Tilt_tilt', '2_Roll_roll',
    '2_Twist_tilt', '2_Twist_roll', '2_Tilt_roll', '2_Shift_shift', '2_Slide_slide', '2_Rise_rise',
    '2_Shift_slide', '2_Shift_rise', '2_Slide_rise', '2_Twist_shift', '2_Twist_slide', '2_Twist_rise',
    '2_Tilt_shift', '2_Tilt_slide', '2_Tilt_rise', '2_Roll_shift', '2_Roll_slide', '2_Roll_rise',
    '2_Slide_stiffness', '2_Shift_stiffness', '2_Roll_stiffness', '2_Tilt_stiffness', '2_Twist_stiffness',
    '2_GC_content', '2_Purine_(AG)_content', '2_Keto_(GT)_content', '2_Adenine_content', '2_Guanine_content',
    '2_Cytosine_content', '2_Thymine_content', '2_Tilt_(DNA-protein_complex)', '2_Roll_(DNA-protein_complex)',
    '2_Shift_(DNA-protein_complex)', '2_Slide_(DNA-protein_complex)', '2_Rise_(DNA-protein_complex)', '2_Shift',
    '2_Slide', '2_Rise', '2_Wedge', '2_Direction', '2_Slide_(RNA)', '2_Rise_(RNA)', '2_Tilt_(RNA)',
    '2_Roll_(RNA)', '2_Twist_(RNA)', '2_Stacking_energy_(RNA)', '2_Rise_stiffness', '2_Melting_Temperature',
    '2_Stacking_energy', '2_Free_energy_(RNA)', '2_Enthalpy_(RNA)', '2_Entropy_(RNA)', '2_Tilt', '2_Roll_y',
    '2_Twist', '2_Flexibility_slide', '2_Flexibility_shift', '2_Enthalpy', '2_Entropy', '2_Free_energy',
    '3_Bend', '3_Tip', '3_Inclination', '3_Major_Groove_Width', '3_Major_Groove_Depth', '3_Major_Groove_Size',
    '3_Major_Groove_Distance', '3_Minor_Groove_Width', '3_Minor_Groove_Depth', '3_Minor_Groove_Size',
    '3_Minor_Groove_Distance', '3_Persistance_Length', '3_Probability_contacting_nucleosome_core',
    '3_Mobility_to_bend_towards_major_groove', '3_Mobility_to_bend_towards_minor_groove', '3_Propeller_Twist',
    '3_Clash_Strength', '3_Shift_(RNA)', '3_Hydrophilicity_(RNA)', '3_Twist_(DNA-protein_complex)',
    '3_Twist_twist', '3_Tilt_tilt', '3_Roll_roll', '3_Twist_tilt', '3_Twist_roll', '3_Tilt_roll',
    '3_Shift_shift', '3_Slide_slide', '3_Rise_rise', '3_Shift_slide', '3_Shift_rise', '3_Slide_rise',
    '3_Twist_shift', '3_Twist_slide', '3_Twist_rise', '3_Tilt_shift', '3_Tilt_slide', '3_Tilt_rise',
    '3_Roll_shift', '3_Roll_slide', '3_Roll_rise', '3_Slide_stiffness', '3_Shift_stiffness', '3_Roll_stiffness',
    '3_Tilt_stiffness', '3_Twist_stiffness', '3_GC_content', '3_Purine_(AG)_content', '3_Keto_(GT)_content',
    '3_Adenine_content', '3_Guanine_content', '3_Cytosine_content', '3_Thymine_content',
    '3_Tilt_(DNA-protein_complex)', '3_Roll_(DNA-protein_complex)', '3_Shift_(DNA-protein_complex)',
    '3_Slide_(DNA-protein_complex)', '3_Rise_(DNA-protein_complex)', '3_Shift', '3_Slide', '3_Rise',
    '3_Wedge', '3_Direction', '3_Slide_(RNA)', '3_Rise_(RNA)', '3_Tilt_(RNA)', '3_Roll_(RNA)', '3_Twist_(RNA)',
    '3_Stacking_energy_(RNA)', '3_Rise_stiffness', '3_Melting_Temperature', '3_Stacking_energy',
    '3_Free_energy_(RNA)', '3_Enthalpy_(RNA)', '3_Entropy_(RNA)', '3_Tilt', '3_Roll_y', '3_Twist',
    '3_Flexibility_slide', '3_Flexibility_shift', '3_Enthalpy', '3_Entropy', '3_Free_energy',
    '4_Bend', '4_Tip', '4_Inclination', '4_Major_Groove_Width', '4_Major_Groove_Depth', '4_Major_Groove_Size',
    '4_Major_Groove_Distance', '4_Minor_Groove_Width', '4_Minor_Groove_Depth', '4_Minor_Groove_Size',
    '4_Minor_Groove_Distance', '4_Persistance_Length', '4_Probability_contacting_nucleosome_core',
    '4_Mobility_to_bend_towards_major_groove', '4_Mobility_to_bend_towards_minor_groove', '4_Propeller_Twist',
    '4_Clash_Strength', '4_Shift_(RNA)', '4_Hydrophilicity_(RNA)', '4_Twist_(DNA-protein_complex)',
    '4_Twist_twist', '4_Tilt_tilt', '4_Roll_roll', '4_Twist_tilt', '4_Twist_roll', '4_Tilt_roll',
    '4_Shift_shift', '4_Slide_slide', '4_Rise_rise', '4_Shift_slide', '4_Shift_rise', '4_Slide_rise',
    '4_Twist_shift', '4_Twist_slide', '4_Twist_rise', '4_Tilt_shift', '4_Tilt_slide', '4_Tilt_rise',
    '4_Roll_shift', '4_Roll_slide', '4_Roll_rise', '4_Slide_stiffness', '4_Shift_stiffness', '4_Roll_stiffness',
    '4_Tilt_stiffness', '4_Twist_stiffness', '4_GC_content', '4_Purine_(AG)_content', '4_Keto_(GT)_content',
    '4_Adenine_content', '4_Guanine_content', '4_Cytosine_content', '4_Thymine_content',
    '4_Tilt_(DNA-protein_complex)', '4_Roll_(DNA-protein_complex)', '4_Shift_(DNA-protein_complex)',
    '4_Slide_(DNA-protein_complex)', '4_Rise_(DNA-protein_complex)', '4_Shift', '4_Slide', '4_Rise', '4_Wedge',
    '4_Direction', '4_Slide_(RNA)', '4_Rise_(RNA)', '4_Tilt_(RNA)', '4_Roll_(RNA)', '4_Twist_(RNA)',
    '4_Stacking_energy_(RNA)', '4_Rise_stiffness', '4_Melting_Temperature', '4_Stacking_energy',
    '4_Free_energy_(RNA)', '4_Enthalpy_(RNA)', '4_Entropy_(RNA)', '4_Tilt', '4_Roll_y', '4_Twist',
    '4_Flexibility_slide', '4_Flexibility_shift', '4_Enthalpy', '4_Entropy', '4_Free_energy'
]
features = {'Conservation': len(cons_features),
            'DNA Shape': len(dna_shape),
            'DNA Properties': len(dna_props),
            'GC Content': len(gcContent),
            'AlphaFold2': len(alphaFold)
            }

print(features)

# Could do some data vis to see which genes have the highest count.
def filter_cosmic(file):
    cos = pd.read_csv(file, sep='\t')
    cos = cos.iloc[:, [0,14, 15, 16, 23, 24]]
    cos.columns = ['gene', 'chrom', 'start', 'stop', 'ref_allele','alt_allele']
    cos.replace(np.nan, '-', regex=True, inplace=True)
    try:
        cos['count'] = cos.groupby(['start', 'stop', 'ref_allele', 'alt_allele'])['start'].transform('size')
    except Exception as e:
        print("Error in groupby operation:", e)
        return None
    # Drop duplicate rows, keeping the first occurrence
    cos = cos.drop_duplicates(subset=['start', 'stop', 'ref_allele', 'alt_allele'], keep='first')
    cos = cos.reset_index(drop=True)

    gene = cos.iloc[:, [0]]  # Select the first column (gene names)

    # Count gene occurrences
    gene_counts = gene['gene'].value_counts()

    # Convert to DataFrame and reset index for sorting
    unique_genes_df = gene_counts.reset_index()
    unique_genes_df.columns = ['gene_name', 'count'] # Rename columns for clarity

    # Sort by count in descending order
    sorted_genes_df = unique_genes_df.sort_values(by='count', ascending=False)
    
    top_genes = sorted_genes_df[sorted_genes_df['count'] > 10]['gene_name'].tolist()
    
    
    # Top genes with mutaiton count > 50
    cos = cos[cos.iloc[:,0].isin(top_genes)].reset_index(drop=True)
    cos['driver_stat'] = [1 for _ in range(len(cos['chrom']))]
    print(Counter(cos['chrom']))
    # Reorder the columns
    # cos = cos[cos['chrom'] == '21']
    cos = cos[['chrom', 'start', 'stop', 'ref_allele', 'alt_allele','gene', 'driver_stat']]
    return cos




