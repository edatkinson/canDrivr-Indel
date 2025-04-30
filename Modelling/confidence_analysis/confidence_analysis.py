# %%
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/edatkinson/Repos/Modelling/new_approach')
from prepare_trainingdata import prepare_data
from train_classifier import train_classifier, process_df, confidence_analysis, plot_metrics, train_classifier_calibrated,process_results, get_confidence, alphaFold, dna_shape,gcContent,dna_props,cons_features
import itertools
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix
import random
random.seed(100)


# %%



# ## CADD TESTING
# df = pd.read_csv('/Users/edatkinson/Repos/Modelling/CADD_data.csv', sep=',')
# df = process_df(df)


# def handle_CADD(cadd):
#     cadd.drop(columns=['WTtrinuc', 'mutTrinuc'], inplace=True)
#     #Combine features together
#     features = list(itertools.chain.from_iterable([alphaFold, dna_shape,gcContent,dna_props,cons_features, ['frameshift_variant', 'inframe_deletion', 'inframe_insertion']]))
#     features.remove('phyloP4way')
#     features.remove('phastCons4way')
#     cadd = cadd[['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'grouping']+features]
#     cadd_inframe = cadd[(cadd['inframe_insertion'] == 1) | (cadd['inframe_deletion'] == 1)].drop(columns=['inframe_insertion', 'inframe_deletion','frameshift_variant'])
#     cadd_frameshift = cadd[cadd['frameshift_variant'] == 1].drop(columns=['frameshift_variant','inframe_insertion', 'inframe_deletion'])
#     print(cadd_frameshift.shape, cadd_inframe.shape)
#     return (cadd_frameshift.dropna(), cadd_inframe)
# cdf_f, cdf_i = handle_CADD(df)

# # %%
# trained_model = final_model_tuple[0]
# y_test = cdf_f['driver_stat']
# y_pred = trained_model.predict(cdf_f.drop(columns=['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'grouping']))

# print(classification_report(y_test, y_pred))

# %%
inframe = '/Users/edatkinson/Repos/Modelling/data/pmed_log_inframe.csv'
frameshift = '/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv'

df = pd.read_csv(frameshift, sep=',')
df, dic = process_df(df, return_chroms=True)
# df = df.drop(columns=['WTtrinuc', 'mutTrinuc', 'genes'])
# %%

params = {
        "learning_rate": 0.2,
        "max_depth": 5,
        "n_estimators": 200,
        "reg_alpha": 0,
        "subsample": 0.5,
        "random_state": 30,
        'colsample_bytree': 0.9
    }

#
actual_predicted, feature_importance, final_model_tuple, results_tuple, data, chrom_performances = train_classifier(
            df, features=None, classifier=XGBClassifier(**params), feature_importance=True
        )


# %%

test_rep, cv_report = process_results(results_tuple, actual_predicted)


# %%

def remap_chroms(dic, chrom_performances):
    inv_dic = {v: k for k, v in dic.items()}
    test_case = inv_dic[1], inv_dic[2], inv_dic[3] # always first 3

    print('random test group =', test_case)

    # Remap chrom_performances keys
    remapped_chrom = {
        inv_dic.get(k, f'chr{k}') if k >3 else 'test': v for k, v in chrom_performances.items()
    }

    return remapped_chrom, test_case, inv_dic

remapped_chrom, test_case, inv_dic = remap_chroms(dic, chrom_performances)


# %%

probabilities_cv = results_tuple[0]['pos_confidence'].to_dict()
probabilities_test = results_tuple[1]['pos_conf'].to_dict()



# combined_probs = {**probabilities_cv, **probabilities_test}

def handle_probs(probs):
    # Ensure the input is a dictionary
    if not isinstance(probs, dict):
        raise TypeError("Expected input to be a dictionary where values are lists of tuples (confidence, index, actual_y, pred_y).")

    # Flatten the data into a list of (confidence, index, actual_y, pred_y) tuples
    flattened_data = []
    for sublist in probs.values():
        if isinstance(sublist, list):  # Ensure it's a list before iterating
            flattened_data.extend(sublist)
    
    # Sort by index
    sorted_data = sorted(flattened_data, key=lambda x: x[1])

    # Convert to a DataFrame
    df = pd.DataFrame(sorted_data, columns=["probability", "Index", "actual", "predicted"])
    df['confidence'] = get_confidence(df['probability'])[0]
    return df



probs_test = handle_probs(probabilities_test)
probs_cv = handle_probs(probabilities_cv)


# %%
test_rep, cv_report = process_results(results_tuple, actual_predicted)


# %%

confidence_df = actual_predicted[1]

confidence_metrics = confidence_analysis(confidence_df)
plot_metrics(confidence_metrics)


# %%
# Define the confidence level
confidence_level = 0.95

# Bootstrapping function to compute confidence intervals
def bootstrap_confidence_interval(data, num_bootstrap_samples=1000, confidence=0.95):
    boot_means = [np.mean(np.random.choice(data, size=len(data), replace=True)) for _ in range(num_bootstrap_samples)]
    lower_bound = np.percentile(boot_means, (1 - confidence) / 2 * 100)
    upper_bound = np.percentile(boot_means, (1 + confidence) / 2 * 100)
    return lower_bound, upper_bound


def calc_prob_CI(confidence_df, confidence_level=0.95):
    # Separate probabilities for positive (1) and negative (0) classes based on actual labels
    probabilities_0 = 1 - confidence_df[confidence_df["actual"] == 0]["probability"].values
    probabilities_1 = confidence_df[confidence_df["actual"] == 1]["probability"].values

    # Compute confidence intervals for each class
    lower_ci_0, upper_ci_0 = bootstrap_confidence_interval(probabilities_0, num_bootstrap_samples=1000, confidence=confidence_level)
    lower_ci_1, upper_ci_1 = bootstrap_confidence_interval(probabilities_1, num_bootstrap_samples=1000, confidence=confidence_level)

    # Display results
    return {
        "Negative Class (0) CI": (lower_ci_0, upper_ci_0),
        "Positive Class (1) CI": (lower_ci_1, upper_ci_1)
    }

# %%
print(calc_prob_CI(probs_test, confidence_level=0.95))
print(calc_prob_CI(probs_cv, confidence_level=0.95))


# %%

def calc_conf_CI(confidence_df, confidence_level=0.95):
    # Extract confidence values for each actual class
    confidence_0 = confidence_df[confidence_df["actual"] == 0]["confidence"].values
    confidence_1 = confidence_df[confidence_df["actual"] == 1]["confidence"].values

    # Compute confidence intervals for the confidence scores
    lower_ci_conf_0, upper_ci_conf_0 = bootstrap_confidence_interval(confidence_0, num_bootstrap_samples=1000, confidence=confidence_level)
    lower_ci_conf_1, upper_ci_conf_1 = bootstrap_confidence_interval(confidence_1, num_bootstrap_samples=1000, confidence=confidence_level)

    # Display results
    return {
        "Neg Confidence CI": (lower_ci_conf_0, upper_ci_conf_0),
        "Pos Confidence CI": (lower_ci_conf_1, upper_ci_conf_1)
    }

print(calc_conf_CI(probs_test, confidence_level=0.95))
print(calc_conf_CI(probs_cv, confidence_level=0.95))

# %%
#Analysis of predictions and their confidence:
def prediction_conf(confidence_df):
    # Plot probability distributions for correct vs. incorrect predictions
    confidence_df["correct"] = confidence_df["actual"] == confidence_df["predicted"]

    # Split into correct and incorrect predictions
    correct_probs = confidence_df[confidence_df["correct"]]["confidence"]
    incorrect_probs = confidence_df[~confidence_df["correct"]]["confidence"]
    
    ci = calc_conf_CI(confidence_df, confidence_level=0.95)['Pos Confidence CI']



    # Plot probability distributions
    plt.figure(figsize=(10, 5))
    sns.histplot(correct_probs, bins=50, kde=True, color="green", label="Correct Predictions", alpha=0.6)
    sns.histplot(incorrect_probs, bins=50, kde=True, color="red", label="Incorrect Predictions", alpha=0.6)
    plt.axvline(x=0.5, color="black", linestyle="--", label="Threshold (0.5)")
    plt.xlabel("Predicted Probability of Class 1")
    plt.ylabel("Count")
    plt.xlim((0.8,1))
    plt.title("Distribution of Prediction Probabilities for Correct vs. Incorrect Predictions")
    plt.legend()
    plt.show()

    # Analyze the confidence levels for incorrect predictions
    incorrect_confidences = confidence_df[~confidence_df["correct"]]["confidence"]
    correct_confidences = confidence_df[confidence_df["correct"]]["confidence"]

    # Compute mean confidence for correct and incorrect predictions
    mean_confidence_correct = correct_confidences.mean()
    mean_confidence_incorrect = incorrect_confidences.mean()

    # Compute proportion of incorrect predictions with high confidence (>0.9)
    high_confidence_errors = (incorrect_confidences > 0.9).mean()

    return {
        "Mean Confidence for Correct Predictions": mean_confidence_correct,
        "Mean Confidence for Incorrect Predictions": mean_confidence_incorrect,
        "Proportion of High Confidence Errors (>0.9)": high_confidence_errors
    }

print(prediction_conf(probs_cv))
print(prediction_conf(probs_test))


# %%
#Adjusted prediction threshold:
#Vary this to minimise False Negatives as this is whats important in medecine. Then maybe plot the number of each in relation to threshold
def confusion_matrix_maker(confidence_df, threshold=0.09):

    # Apply new threshold to generate adjusted predictions
    confidence_df["adjusted_prediction"] = (confidence_df["probability"] >= threshold).astype(int)

    # Compute classification metrics
    report = classification_report(confidence_df["actual"], confidence_df["adjusted_prediction"], output_dict=True)
    conf_matrix = confusion_matrix(confidence_df["actual"], confidence_df["adjusted_prediction"])

    # Convert classification report to DataFrame for better visualization
    report_df = pd.DataFrame(report).transpose()

    print(conf_matrix)
    return report_df

# confusion_matrix_maker(probs_cv, threshold=0.5)
confusion_matrix_maker(probs_test, threshold=0.5)
confusion_matrix_maker(probs_test, threshold=0.1)

# %%



