# %%
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/edatkinson/Repos/Modelling/new_approach')
from prepare_trainingdata import prepare_data
from train_classifier import train_classifier, train_and_evaluate, process_df, confidence_analysis, plot_metrics, train_classifier_calibrated,process_results, get_confidence, alphaFold, dna_shape,gcContent,dna_props,cons_features, print_res
import itertools
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix

def predict_with_trained_model(independent_test_file: str, trained_model, features_ind: list, pmed=False):
    """
    Loads an independent test dataset, preprocesses it, and makes predictions using a provided trained model.

    Args:
        independent_test_file (str): Path to the independent test dataset CSV file.
        trained_model: A trained model object with a 'predict' method.
        features_ind (list): List of feature names to be used for prediction.

    Returns:
        tuple: A tuple containing the true labels (y_test) and the predicted labels (y_pred).
    """
    ind = pd.read_csv(independent_test_file, sep=',')
    ind = process_df(ind).dropna()
    ind = ind.drop(columns=['Unnamed: 0','genes'])
    if pmed == True:
        ind = ind[['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'grouping']+features_ind+['pmed']]
    else:
        ind = ind[['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'grouping']+features_ind]
    
    y_test = ind['driver_stat']
    y_pred = trained_model.predict(ind.drop(columns=['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'grouping']))

    return y_test, y_pred

def combine(inframe, frameshift):
    infr = pd.read_csv(inframe, sep=',')
    frame = pd.read_csv(frameshift, sep=',')
    combined = pd.concat([infr,frame])
    print(combined.shape)
    combined.to_csv('/Users/edatkinson/Repos/Modelling/data/combined_data.csv', index=None)
    return combined

if __name__ == '__main__':
    pmed = False
    verbose = True
    features = [gcContent,cons_features]

    cadd_file = '/Users/edatkinson/Repos/Modelling/CADD_data.csv'
    inframe = '/Users/edatkinson/Repos/Modelling/data/pmed_log_inframe.csv'
    frameshift = '/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv'
    ind_test  = '/Users/edatkinson/Repos/Modelling/Merged_Ind_test.csv'
    combined = '/Users/edatkinson/Repos/Modelling/data/combined_data.csv'

    actual_predicted, feature_importance, final_model_tuple, results_tuple = train_and_evaluate(frameshift, features, verbose=verbose, pmed=pmed)
    features_ind = list(itertools.chain.from_iterable(features))
    trained_model = final_model_tuple[0]
    y_test, y_pred = predict_with_trained_model(ind_test, trained_model, features_ind, pmed=pmed)
    message = f'Pmed: {pmed}'
    print_res(classification_report(y_test, y_pred), message)
# %%
