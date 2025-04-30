# import pandas as pd
# import numpy as np
# from sklearn.model_selection import LeaveOneGroupOut, train_test_split, GridSearchCV
# from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, balanced_accuracy_score, matthews_corrcoef
# from imblearn.under_sampling import RandomUnderSampler
# from xgboost import XGBClassifier
# import warnings

# warnings.filterwarnings("ignore")

# def drop_bias(data: pd.DataFrame):
#     """Removes features causing bias and unnecessary object columns."""
#     bias_cols = ['merged']
#     object_cols = data.select_dtypes(include=['object']).columns
    
#     if len(object_cols) > 1:
#         data = data.drop(columns=object_cols[-2:])
    
#     return data.drop(columns=bias_cols, errors='ignore')

# def balance_data(data: pd.DataFrame):
#     """Balances simulated and real variants for each chromosome."""
#     human_derived = data[data['driver_stat'] == 0]
#     simulated_indels = data[data['driver_stat'] == 1]
#     count_chrom = human_derived.groupby('chrom').size()
#     sampled_simulated = []

#     for chrom, count in count_chrom.items():
#         sim_chrom = simulated_indels[simulated_indels['chrom'] == chrom]
#         sampled_simulated.append(sim_chrom.sample(n=min(count, len(sim_chrom)), replace=False))
    
#     return pd.concat([human_derived, pd.concat(sampled_simulated)], ignore_index=True)

# def train_and_evaluate_model(classifier, X_train, y_train, X_test, y_test, groups):
#     """Trains and evaluates the model with Leave-One-Group-Out CV and hyperparameter tuning."""
#     logo = LeaveOneGroupOut()
    
#     param_grid = {
#         'learning_rate': [0.01, 0.1, 0.2],
#         'n_estimators': [100, 200, 300],
#         'max_depth': [3, 5, 7]
#     }
    
#     grid_search = GridSearchCV(classifier, param_grid, cv=3, scoring='roc_auc', n_jobs=-1)
    
#     best_model = None
#     best_auc = -1
    
#     for train_idx, val_idx in logo.split(X_train, y_train, groups):
#         X_tr, X_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
#         y_tr, y_val = y_train.iloc[train_idx], y_train.iloc[val_idx]

#         if len(np.unique(y_val)) < 2:
#             print(f"Skipping LOGO fold with chromosome {groups[val_idx][0]} (only one class present).")
#             continue  # Skip this fold
        
#         sampler = RandomUnderSampler()
#         X_tr_res, y_tr_res = sampler.fit_resample(X_tr, y_tr)
        
#         grid_search.fit(X_tr_res, y_tr_res)
#         model = grid_search.best_estimator_
        
#         y_val_pred = model.predict_proba(X_val)[:, 1]
#         auc = roc_auc_score(y_val, y_val_pred)
        
#         if auc > best_auc:
#             best_auc = auc
#             best_model = model
    
#     y_test_pred = best_model.predict(X_test)
#     test_results = {
#         "balanced_accuracy": balanced_accuracy_score(y_test, y_test_pred),
#         "precision": precision_score(y_test, y_test_pred),
#         "recall": recall_score(y_test, y_test_pred),
#         "f1": f1_score(y_test, y_test_pred),
#         "roc_auc": roc_auc_score(y_test, y_test_pred),
#         "mcc": matthews_corrcoef(y_test, y_test_pred)
#     }
    
#     return best_model, test_results

# class CandrivrModel:
#     def __init__(self, data_file):
#         self.data = pd.read_csv(data_file, sep=',')
#         self.cleaned_data = drop_bias(self.data)
#         self.sampled_data = balance_data(self.cleaned_data)
#         self.classifier = XGBClassifier(eval_metric='logloss')

#     def prepare_data(self):
#         X = self.sampled_data.drop(columns=['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat'])
#         y = self.sampled_data['driver_stat']
#         groups = self.sampled_data['chrom'].values
#         unique_groups = np.unique(groups)
#         np.random.shuffle(unique_groups)

#         test_groups, test_size, target_size = [], 0, int(0.2 * len(y))
#         for chrom in unique_groups:
#             test_groups.append(chrom)
#             test_size += (groups == chrom).sum()
#             if test_size >= target_size:
#                 break
        
#         test_indices = np.isin(groups, test_groups)
#         X_train, X_test = X[~test_indices], X[test_indices]
#         y_train, y_test = y[~test_indices], y[test_indices]
#         groups_train = groups[~test_indices]
        
#         return X_train, y_train, X_test, y_test, groups_train
    
#     def run(self):
#         X_train, y_train, X_test, y_test, groups_train = self.prepare_data()
#         best_model, test_results = train_and_evaluate_model(self.classifier, X_train, y_train, X_test, y_test, groups_train)
#         return best_model, test_results

# # Execute
# candrivr = CandrivrModel('merged_annotated_data.csv')
# best_model, test_results = candrivr.run()
# print("Final Model Performance:", test_results)


import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneGroupOut, train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, balanced_accuracy_score, matthews_corrcoef
from imblearn.under_sampling import RandomUnderSampler
from xgboost import XGBClassifier
from boruta import BorutaPy
import warnings

warnings.filterwarnings("ignore")

def drop_bias(data: pd.DataFrame):
    """Removes features causing bias and unnecessary object columns."""
    bias_cols = ['merged', 'atom_site_occupancy']
    object_cols = data.select_dtypes(include=['object']).columns
    
    if len(object_cols) > 1:
        data = data.drop(columns=object_cols[-2:])
    
    return data.drop(columns=bias_cols, errors='ignore')

def balance_data(data: pd.DataFrame):
    """Balances simulated and real variants for each chromosome.
    Was made for CADD data, but can be adapted for other datasets.
    """
    human_derived = data[data['driver_stat'] == 0]
    simulated_indels = data[data['driver_stat'] == 1]
    count_chrom = human_derived.groupby('chrom').size()
    sampled_simulated = []

    for chrom, count in count_chrom.items():
        sim_chrom = simulated_indels[simulated_indels['chrom'] == chrom]
        sampled_simulated.append(sim_chrom.sample(n=min(count, len(sim_chrom)), replace=False))
    
    return pd.concat([human_derived, pd.concat(sampled_simulated)], ignore_index=True)

def select_features_with_boruta(X_train, y_train):
    """Performs feature selection using Boruta."""
    xgb = XGBClassifier(n_jobs=-1, eval_metric='logloss')
    boruta = BorutaPy(xgb, n_estimators='auto', verbose=2, random_state=42)
    boruta.fit(X_train.values, y_train.values)
    
    selected_features = X_train.columns[boruta.support_]
    return X_train[selected_features], selected_features

def train_and_evaluate_model(classifier, X_train, y_train, X_test, y_test, groups):
    """Trains and evaluates the model with Leave-One-Group-Out CV."""
    logo = LeaveOneGroupOut()
    best_model = None
    best_auc = -1
    
    for train_idx, val_idx in logo.split(X_train, y_train, groups):
        X_tr, X_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_tr, y_val = y_train.iloc[train_idx], y_train.iloc[val_idx]
        print(X_tr.shape, X_val.shape)
        if len(np.unique(y_val)) < 2:
            print(f"Skipping LOGO fold with chromosome {groups[val_idx][0]} (only one class present).")
            continue  # Skip this fold
        
        sampler = RandomUnderSampler(random_state=42)
        X_tr_res, y_tr_res = sampler.fit_resample(X_tr, y_tr)
        
        classifier.fit(X_tr_res, y_tr_res)
        y_val_pred = classifier.predict_proba(X_val)[:, 1]
        auc = roc_auc_score(y_val, y_val_pred)
        
        if auc > best_auc:
            best_auc = auc
            best_model = classifier
    
    y_test_pred = best_model.predict(X_test)
    test_results = {
        "balanced_accuracy": balanced_accuracy_score(y_test, y_test_pred),
        "precision": precision_score(y_test, y_test_pred),
        "recall": recall_score(y_test, y_test_pred),
        "f1": f1_score(y_test, y_test_pred),
        "roc_auc": roc_auc_score(y_test, y_test_pred),
        "mcc": matthews_corrcoef(y_test, y_test_pred)
    }
    
    return best_model, test_results

class CandrivrModel:
    def __init__(self, data_file):
        self.data = pd.read_csv(data_file, sep=',')
        self.cleaned_data = drop_bias(self.data)
        self.sampled_data = balance_data(self.cleaned_data)
        self.classifier = XGBClassifier(eval_metric='logloss')

    def prepare_data(self):
        X = self.sampled_data.drop(columns=['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat'])
        y = self.sampled_data['driver_stat']
        X.replace(np.nan, 0, inplace=True)
        groups = self.sampled_data['chrom'].values
        unique_groups = np.unique(groups)
        np.random.shuffle(unique_groups)

        test_groups, test_size, target_size = [], 0, int(0.2 * len(y))
        for chrom in unique_groups:
            test_groups.append(chrom)
            test_size += (groups == chrom).sum()
            if test_size >= target_size:
                break
        
        test_indices = np.isin(groups, test_groups)
        X_train, X_test = X[~test_indices], X[test_indices]
        y_train, y_test = y[~test_indices], y[test_indices]
        groups_train = groups[~test_indices]
        
        # Perform feature selection
        X_train_selected, selected_features = select_features_with_boruta(X_train, y_train)
        X_test_selected = X_test[selected_features]

        return X_train_selected, y_train, X_test_selected, y_test, groups_train, selected_features
    
    def run(self):
        X_train, y_train, X_test, y_test, groups_train, selected_features = self.prepare_data()
        print(f"Selected Features: {list(selected_features)}")
        best_model, test_results = train_and_evaluate_model(self.classifier, X_train, y_train, X_test, y_test, groups_train)
        return best_model, test_results

# Execute
# candrivr = CandrivrModel('/Users/edatkinson/Repos/Modelling/merged_annotated_data.csv')
# best_model, test_results = candrivr.run() # uses full dataset, not equal samples
# print("Final Model Performance:", test_results)

# # Now I have my best model, I need to get an independent test set to find it's actual accuracy.
# best_model.save_model('xgboost_feature_selected_model.json')
