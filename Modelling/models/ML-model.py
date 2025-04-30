'''
File for my machine learning model which:
- Processes and splits data
- Trains a model
- Performs LOCO-CV
- Obtains a Confidence measure
- Performs Feature Selection
- Performs Visualisation of key metrics & data
- Think I'll do it all in one file.
'''
import pandas as pd
import xgboost as xgb
from xgboost import XGBClassifier
import numpy as np
from sklearn.model_selection import train_test_split, KFold, LeaveOneGroupOut, cross_val_score
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt
from imblearn.under_sampling import RandomUnderSampler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.metrics import balanced_accuracy_score



def drop_bias(data: pd.DataFrame):
    """ 'WTtrinuc', 'mutTrinuc' are also objects """
    bias_cols = ['merged']
    object_cols = data.select_dtypes(include=['object']).columns

    if not object_cols.empty:
        print(f"Columns with dtype('O'): {object_cols[-2:]}")
        data = data.drop(columns=object_cols[-2:])
        print("Object columns removed.")

    data = data.drop(columns=bias_cols)
    return data
    
def select_equal_amounts(data: pd.DataFrame):
    """
    Randomly select an equal number of simulated variants for each chromosome
    based on the counts of human-derived variants.
    
    Args:
        data (pd.DataFrame): Input DataFrame with columns 'driver_stat' and 'chrom'.
        
    Returns:
        pd.DataFrame: A DataFrame of simulated variants sampled to match human-derived counts.
    """
    human_derived = data[data['driver_stat'] == 0]
    simulated_indels = data[data['driver_stat'] == 1]

    count_chrom = human_derived.groupby('chrom').size()
    sampled_simulated = []

    for chrom, count in count_chrom.items():
        sim_chrom = simulated_indels[simulated_indels['chrom'] == chrom]
        
        if count <= len(sim_chrom):
            sampled_simulated.append(sim_chrom.sample(n=count, random_state=42))
        else:
            sampled_simulated.append(sim_chrom)

    result = pd.concat(sampled_simulated, ignore_index=True)
    result = pd.concat([human_derived, result], ignore_index=True)
    return result


def train_and_evaluate_model(classifier, X_train_val, y_train_val, groups_train_val, X_test, y_test):
    """
    Train the model using Leave-One-Group-Out cross-validation and evaluate its performance.

    Parameters:
        classifier: The classifier model object.
        X_train_val (pd.DataFrame): Features of the training/validation set.
        y_train_val (pd.Series): Target variable of the training/validation set.
        groups_train_val (pd.Series): Groups corresponding to the training/validation set.
        X_test (pd.DataFrame): Features of the test set.
        y_test (pd.Series): Target variable of the test set.

    Returns:
        tuple: A tuple containing metrics, actual_targets, predicted_targets, final_model, test_accuracy, y_test_pred.
    """

    metrics = {
        "accuracy": [],
        "precision": [],
        "recall": [],
        "f1": [],
        "roc_auc": [],
    }

    # Leave-One-Group-Out cross-validation
    logo = LeaveOneGroupOut()
    actual_targets = np.array([])
    predicted_targets = np.array([])
    for train_index, val_index in logo.split(X_train_val, y_train_val, groups_train_val):
        model = classifier
        X_train, X_val = X_train_val.iloc[train_index], X_train_val.iloc[val_index]
        y_train, y_val = y_train_val.iloc[train_index], y_train_val.iloc[val_index]
        # Random undersampling to balance the training data
        sampler = RandomUnderSampler(random_state=42)
        X_train_resampled, y_train_resampled = sampler.fit_resample(X_train, y_train)
        print(X_train_resampled.shape)
        # Model fitting
        model.fit(X_train_resampled, y_train_resampled)

        # Predict on validation set
        y_val_pred = model.predict(X_val)
        actual_targets = np.append(actual_targets, y_val)
        predicted_targets = np.append(predicted_targets, y_val_pred)

        # Evaluate and store metrics
        y_val_pred_proba = model.predict_proba(X_val)[:, 1]
        metrics["accuracy"].append(accuracy_score(y_val, y_val_pred))
        metrics["precision"].append(precision_score(y_val, y_val_pred))
        metrics["recall"].append(recall_score(y_val, y_val_pred))
        metrics["f1"].append(f1_score(y_val, y_val_pred))
        metrics["roc_auc"].append(roc_auc_score(y_val, y_val_pred_proba))

    # Fit the final model on the entire training/validation set
    final_model = classifier
    # final_sampler = RandomUnderSampler(random_state=42)
    # X_train_val_resampled, y_train_val_resampled = final_sampler.fit_resample(X_train_val, y_train_val)
    final_model.fit(X_train_val, y_train_val)

    # Predict and evaluate on the test set
    y_test_pred = final_model.predict(X_test[X_train_val.columns.tolist()])
    test_accuracy = balanced_accuracy_score(y_test, y_test_pred)

    test_results = pd.DataFrame(
        {
            "balanced_accuracy": [test_accuracy],
            "precision": [precision_score(y_test, y_test_pred)],
            "recall": [recall_score(y_test, y_test_pred)],
            "f1": [f1_score(y_test, y_test_pred)],
            "roc_auc": [roc_auc_score(y_test, y_test_pred)],
        }
    )

    actual_predicted_targets_cross_val = (actual_targets, predicted_targets)
    actual_predicted_targets_test = (y_test, y_test_pred)
    actual_predicted = (actual_predicted_targets_cross_val, actual_predicted_targets_test)

    return metrics, final_model, test_results, actual_predicted


class Data:
    def __init__(self, data_file: str):
        self.data_file = data_file
        self.data = pd.read_csv(data_file, sep=',')
        self.clean_data = drop_bias(self.data)
        self.sampled_data = select_equal_amounts(self.clean_data)
    def __str__(self):
         return str(self.sampled_data[self.sampled_data['driver_stat'] == 1].shape,self.sampled_data[self.sampled_data['driver_stat'] == 0].shape)
    


class Model(Data):
    def __init__(self,classifier: XGBClassifier, data_file: str):
        self.data = Data(data_file)

        self.classifier = classifier
        
    def target_and_label(self, data: pd.DataFrame):
        """Split data into train, val, test"""
        print(data)
        X = data.drop(columns=['chrom','pos', 'ref_allele', 'alt_allele', 'driver_stat'])

        y = data["driver_stat"]

        return X, y
    
    def get_idx(self, data, X, y):
        """Define chromosomal groups"""
        groups = data['chrom'].values
        unique_groups = np.unique(groups)
        np.random.shuffle(unique_groups)

        test_groups = []
        test_size = 0
        target_size = int(0.2 * len(data))

        for chrom in unique_groups:
            test_groups.append(chrom)
            test_size += (groups == chrom).sum()
            if test_size >= target_size:
                break
        
        test_indices = np.isin(groups, test_groups)
        train_val_indices = ~test_indices
        test_indices, train_val_indices, groups

        X_train_val, y_train_val = X[train_val_indices], y[train_val_indices]
        X_test, y_test = X[test_indices], y[test_indices]
        groups_train_val = groups[train_val_indices]

        return X_train_val, y_train_val, X_test, y_test, groups_train_val

    def train_evaluate(self, data: pd.DataFrame):
        X, y = self.target_and_label(data)
        X_train_val, y_train_val, X_test, y_test, groups_train_val = self.get_idx(data, X, y)
        metrics, final_model, test_results, actual_predicted = train_and_evaluate_model(self.classifier, X_train_val, y_train_val, groups_train_val, X_test, y_test)

        return metrics, final_model, test_results, actual_predicted
    

data_file = 'merged_annotated_data.csv'

candrivr = Model(XGBClassifier(), data_file)

metrics, final_model, test_results, actual_predicted = candrivr.train_evaluate(candrivr.data.sampled_data)

