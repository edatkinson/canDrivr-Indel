from sklearn.model_selection import train_test_split, LeaveOneGroupOut
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, classification_report
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from imblearn.under_sampling import RandomUnderSampler
from sklearn.inspection import permutation_importance
import pandas as pd
import numpy as np
import time
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from collections import Counter, defaultdict

def train_and_evaluate_model(classifier, X_train_val, y_train_val, groups_train_val, X_test, y_test, seed=20):
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
        'pos_confidence': [],
        'neg_confidence': []
    }



    # Leave-One-Group-Out cross-validation
    logo = LeaveOneGroupOut()
    actual_targets = np.array([])
    predicted_targets = np.array([])
    chrom_performances = {}
    for train_index, val_index in logo.split(X_train_val, y_train_val, groups_train_val):
        group = groups_train_val[val_index[0]]
        model = classifier
        X_train, X_val = X_train_val.iloc[train_index], X_train_val.iloc[val_index]
        y_train, y_val = y_train_val.iloc[train_index], y_train_val.iloc[val_index]
        # Random undersampling to balance the training data
        sampler = RandomUnderSampler(random_state=42)
        X_train_resampled, y_train_resampled = sampler.fit_resample(X_train, y_train)
        # print(y_train_resampled.value_counts())
        # Model fitting
        model.fit(X_train_resampled, y_train_resampled)

        # Predict on validation set
        y_val_pred = model.predict(X_val)
        # print("y_val_pred", Counter(y_val_pred))
        # print("y_val", Counter(y_val))
        actual_targets = np.append(actual_targets, y_val)
        predicted_targets = np.append(predicted_targets, y_val_pred)
        

        # Evaluate and store metrics
        y_val_pred_proba = model.predict_proba(X_val)[:, 1]
        y_val_pred_proba_neg = model.predict_proba(X_val)[:,0]


        metrics['pos_confidence'].append([*zip(y_val_pred_proba, val_index, y_val, y_val_pred)])
        metrics['neg_confidence'].append([*zip(y_val_pred_proba_neg, val_index, y_val, y_val_pred)])
        metrics["accuracy"].append(accuracy_score(y_val, y_val_pred))
        metrics["precision"].append(precision_score(y_val, y_val_pred))
        metrics["recall"].append(recall_score(y_val, y_val_pred, zero_division=0))
        metrics["f1"].append(f1_score(y_val, y_val_pred))
        metrics["roc_auc"].append([roc_auc_score(y_val, y_val_pred_proba) if len(np.unique(y_val)) == 2 else np.nan])

        chrom_performances[group] = {
            "accuracy": accuracy_score(y_val, y_val_pred),
            "precision": precision_score(y_val, y_val_pred),
            "recall": recall_score(y_val, y_val_pred, zero_division=0),
            "f1": f1_score(y_val, y_val_pred),
            "roc_auc": roc_auc_score(y_val, y_val_pred_proba) if len(np.unique(y_val)) == 2 else np.nan,
            "count": len(y_val),
            "actual": y_val,
            "pred": y_val_pred,
            'conf': [*zip(y_val_pred_proba, val_index, y_val, y_val_pred)]
        }

    # Fit the final model on the entire training/validation set
    final_model = classifier
    #Under Sampling
    final_sampler = RandomUnderSampler(random_state=42)
    X_train_val_resampled, y_train_val_resampled = final_sampler.fit_resample(X_train_val, y_train_val)
    # print(f'ytrain: {y_train_resampled.value_counts()}')
    #Fit the model using the whole train+val set
    final_model.fit(X_train_val_resampled, y_train_val_resampled)
    # print(y_test.value_counts())
    # Predict and evaluate on the test set: 
    y_test_pred = final_model.predict(X_test[X_train_val.columns.tolist()])
    y_test_pred_prob = final_model.predict_proba(X_test[X_train_val.columns.tolist()])[:,1]
    y_test_pred_prob_neg = final_model.predict_proba(X_test[X_train_val.columns.tolist()])[:,0]
    # print(Counter(y_test_pred), Counter(y_test))
    test_accuracy = balanced_accuracy_score(y_test, y_test_pred)


    test_results = pd.DataFrame(
        {
            "balanced_accuracy": [test_accuracy],
            "precision": [precision_score(y_test, y_test_pred)],
            "recall": [recall_score(y_test, y_test_pred)],
            "f1": [f1_score(y_test, y_test_pred)],
            "roc_auc": [roc_auc_score(y_test, y_test_pred) if len(np.unique(y_test)) == 2 else np.nan],
            "actual": [Counter(y_test)],
            "pred": [Counter(y_test_pred)],
            "neg_conf": [[*zip(y_test_pred_prob_neg, X_test.index, y_test, y_test_pred)]],
            "pos_conf": [[*zip(y_test_pred_prob, X_test.index, y_test, y_test_pred)]]

        }
    )
    
    chrom_performances[1] = {
            "accuracy": accuracy_score(y_test, y_test_pred),
            "precision": precision_score(y_test, y_test_pred),
            "recall": recall_score(y_test, y_test_pred, zero_division=0),
            "f1": f1_score(y_test, y_test_pred),
            "roc_auc": roc_auc_score(y_test, y_test_pred_prob) if len(np.unique(y_val)) == 2 else np.nan,
            "count": len(y_test),
            "actual": y_test,
            "pred": y_test_pred,
            "conf": [*zip(y_test_pred_prob, X_test.index, y_test, y_test_pred)]
        }

    actual_predicted_targets_cross_val = (actual_targets, predicted_targets)
    actual_predicted_targets_test = (y_test, y_test_pred)
    actual_predicted = (actual_predicted_targets_cross_val, actual_predicted_targets_test)


    return metrics, final_model, test_results, actual_predicted, chrom_performances
# %%






