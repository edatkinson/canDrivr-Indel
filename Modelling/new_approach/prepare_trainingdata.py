from sklearn.model_selection import train_test_split, LeaveOneGroupOut
from imblearn.under_sampling import RandomUnderSampler
import numpy as np
import time
import random
import pandas as pd
from collections import Counter



def prepare_data(df, features=None):
    """
    Prepare the data for training.

    Parameters:
        df (pd.DataFrame): The input dataframe containing the dataset.
        features (list, optional): List of feature columns. If None, default columns are used. Default is None.

    Returns:
        tuple: A tuple containing X_train_val, y_train_val, groups_train_val, X_test, y_test, groups_test.
    """
    if features is None:
        X = df.drop(["chrom", "pos", "ref_allele", "alt_allele", "driver_stat"], axis=1)
    elif isinstance(features, list):
        if len(features) > 1:
            X = df[features]
        elif len(features) == 1:
            X = df[features]

    y = df["driver_stat"]
    groups = df["grouping"].values

    # Get unique groups
    unique_groups = np.unique(groups)
    # print(unique_groups)

    # Define the test group
    test_group = unique_groups[:3] # first group, where each group comprises of 2 chromosomes

    # Define train_val groups
    train_val_groups = unique_groups[3:] #remaining 18 groups

    # Filter test and train_val groups
    # test_indices = groups == test_group
    test_indices = np.isin(groups, test_group)

    train_val_indices = np.isin(groups, train_val_groups)

    # X = X.drop("grouping", axis = 1)

    # Define test and train_val data
    X_test, y_test = X[test_indices], y[test_indices]
    X_train_val, y_train_val = X[train_val_indices], y[train_val_indices]
    groups_test, groups_train_val = groups[test_indices], groups[train_val_indices]

    return X_train_val, y_train_val, groups_train_val, X_test, y_test, groups_test

# groups = list(range(1, 12)) * 2  # Duplicating group range to cover 22 assignments
# groups = list(range(1,23))
# print(groups)
# random.shuffle(groups)  # Shuffle to randomize distribution
# chromosomes = [
#     'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
#     'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
#     'chr18', 'chr19', 'chr20','chr21', 'chr22'
# ]


# chromosome_to_group = dict(zip(chromosomes, groups))
# print(chromosome_to_group)



# df = pd.read_csv('/Users/edatkinson/Repos/Modelling/merged_annotated_data.csv', sep=',')
# df['grouping'] = df['chrom'].map(chromosome_to_group)
# # print(np.unique(df['grouping']))
# print(np.unique(df['chrom']))

# df.drop(columns=['WTtrinuc', 'mutTrinuc'], inplace=True)
# df.replace(np.nan, 0, inplace=True)



# X_train_val, y_train_val, groups_train_val, X_test, y_test, groups_test = prepare_data(df)

# print(Counter(groups_train_val))
# print(Counter(groups_test))
