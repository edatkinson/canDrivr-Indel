

# %%
import sys
# sys.path.append('/Users/edatkinson/Repos/Modelling/new_approach')
from prepare_trainingdata import prepare_data
from train_evaluate import train_and_evaluate_model
from train_classifier import process_df
import random
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from sklearn.utils import resample
from sklearn.model_selection import train_test_split
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_curve, roc_auc_score
from xgboost import XGBClassifier
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
import matplotlib.patches as mpatches
from train_classifier import train_classifier, process_df, alphaFold, dna_shape,gcContent,dna_props,cons_features
random.seed(100)


# %%

# Load dataset
df = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv', sep=',')
df, remapped_chroms = process_df(df, return_chroms=True)
# print(remapped_chroms)
# Split data into training and testing
selected_features_boruta = [ # Add/Remove pmed for changes in prediction
    'phyloP17way',
    'phyloP100way',
    'phastCons20way',
    'phyloP470way',
    'phastCons17way',
    'phastCons7way',
    'phyloP4way',
    'k24.Umap.MultiTrackMappability',
    'k50.Umap.MultiTrackMappability',
    '200GCContent',
    '200CpGCount',
    '1000CpGCount',
    '9_MGW',
    '17_HelT',
    '8_ProT',
    '9_ProT',
    '1_Shift_stiffness',
    '1_Stacking_energy_(RNA)',
    '2_Shift_shift',
    '3_Shift_slide',
    '4_Twist_slide',
    '4_Stacking_energy_(RNA)',
    'pmed'
]


X_train_val, y_train_val, groups_train_val, X_test, y_test, groups_test = prepare_data(df, selected_features_boruta)

def get_confidence(probs):
    confidence = np.maximum(probs, 1 - probs)
    probs_safe = np.clip(probs, 1e-10, 1 - 1e-10)
    h_p = -probs_safe * np.log2(probs_safe) - (1 - probs_safe) * np.log2(1 - probs_safe)
    entropy = 1 - h_p  # Normalized entropy
    return confidence, entropy

def train_classifier_with_bootstrap(X_train, y_train, X_test, y_test, n_bootstraps=100):
    proba_predictions = np.zeros((n_bootstraps, len(X_test)))

    for i in range(n_bootstraps):

        X_resampled, y_resampled = resample(X_train.drop(columns=['pmed']), y_train, replace=True, random_state=i)
        boot_model = XGBClassifier(n_estimators=200, learning_rate=0.2, max_depth=5, reg_alpha=0, subsample=1, random_state=i)
        boot_model.fit(X_resampled, y_resampled)

        test_probs = boot_model.predict_proba(X_test.drop(columns=['pmed']))[:, 1]

        proba_predictions[i] = test_probs
        
    mean_proba = np.mean(proba_predictions, axis=0)
    std_proba = np.std(proba_predictions, axis=0)
    lower_bound = np.percentile(proba_predictions, 2.5, axis=0)
    upper_bound = np.percentile(proba_predictions, 97.5, axis=0)
    
    return mean_proba, std_proba, lower_bound, upper_bound

def train_classifier_with_bootstrap_iso(X_train, y_train, X_test, y_test, n_bootstraps=100):
    proba_predictions = np.zeros((n_bootstraps, len(X_test)))

    for i in range(n_bootstraps):
        X_train, X_calib, y_train, y_calib = train_test_split(X_train_val.drop(columns=['pmed']), y_train_val, test_size=0.2, random_state=42)

        X_resampled, y_resampled = resample(X_train, y_train, replace=True, random_state=i)
        boot_model = XGBClassifier(n_estimators=200, learning_rate=0.2, max_depth=5, reg_alpha=0, subsample=1, random_state=i)
        boot_model.fit(X_resampled, y_resampled)

        calib_probs = boot_model.predict_proba(X_calib)[:, 1]

        iso_reg = IsotonicRegression(out_of_bounds='clip')
        iso_reg.fit(calib_probs, y_calib)
        test_probs = boot_model.predict_proba(X_test.drop(columns=['pmed']))[:, 1]

        proba_predictions[i] = iso_reg.predict(test_probs)
        
    mean_proba = np.mean(proba_predictions, axis=0)
    std_proba = np.std(proba_predictions, axis=0)
    lower_bound = np.percentile(proba_predictions, 2.5, axis=0)
    upper_bound = np.percentile(proba_predictions, 97.5, axis=0)
    
    return mean_proba, std_proba, lower_bound, upper_bound



def train_classifier_with_bootstrap_platt(X_train_val, y_train_val, X_test, y_test, n_bootstraps=100):
    proba_predictions = np.zeros((n_bootstraps, len(X_test)))

    for i in range(n_bootstraps):
        # Split for calibration (you had this *inside* the loop but overwriting X_train — moved outside)
        X_train_split, X_calib_split, y_train_split, y_calib_split = train_test_split(
            X_train_val.drop(columns=['pmed']), y_train_val, test_size=0.2, random_state=i
        )

        # Bootstrap from the calibration training set
        X_resampled, y_resampled = resample(X_train_split, y_train_split, replace=True, random_state=i)

        # Train base model
        boot_model = XGBClassifier(
            n_estimators=200,
            learning_rate=0.2,
            max_depth=5,
            reg_alpha=0,
            subsample=0.6,
            random_state=i
        )
        boot_model.fit(X_resampled, y_resampled)

        # Wrap in CalibratedClassifierCV (Platt scaling with sigmoid)
        calibrated_clf = CalibratedClassifierCV(boot_model, method='sigmoid', cv='prefit')  # Important: prefit
        calibrated_clf.fit(X_calib_split, y_calib_split)

        # Predict calibrated probabilities on the test set
        probs_calibrated = calibrated_clf.predict_proba(X_test.drop(columns=['pmed']))[:, 1]
        proba_predictions[i] = probs_calibrated

    # Aggregate bootstrap results
    mean_proba = np.mean(proba_predictions, axis=0)
    std_proba = np.std(proba_predictions, axis=0)
    lower_bound = np.percentile(proba_predictions, 2.5, axis=0)
    upper_bound = np.percentile(proba_predictions, 97.5, axis=0)

    return mean_proba, std_proba, lower_bound, upper_bound


def plot_calibration_curve(y_true, y_proba, n_bins=30, label='Platt-Scaled'):
    prob_true, prob_pred = calibration_curve(y_true, y_proba, n_bins=n_bins, strategy='quantile')
    
    return prob_true, prob_pred
# %%
def plot_confidence_intervals(y_test, mean_proba, lower_ci, upper_ci):
    plt.figure(figsize=(8, 6))

    # Error bar calculation (asymmetric)
    lower_error = mean_proba - lower_ci
    upper_error = upper_ci - mean_proba
    asymmetric_error = [lower_error, upper_error]

    # Set colors based on class labels
    colors = ['red' if y == 1 else 'blue' for y in y_test]

    # Plot each point with its individual color and error bar
    for i in range(len(mean_proba)):
        plt.errorbar(i, mean_proba[i],
                     yerr=[[lower_error[i]], [upper_error[i]]],
                     fmt='o', color=colors[i], ecolor='grey', elinewidth=1.5, capsize=3, alpha=0.7)

    # plt.hlines([0.8, 0.2], 0, len(mean_proba), colors='grey', linestyles='dashed', alpha=0.7, label='CI width')
    plt.hlines([0.5], 0, len(mean_proba), colors='green', linestyles='dashed', alpha=0.7, label='Decision Threshold')
    # Legend
    red_patch = mpatches.Patch(color='red', label='True Driver')
    blue_patch = mpatches.Patch(color='blue', label='True Passenger')
    ci_patch = mpatches.Patch(color='grey', alpha=0.5, label='95% CI (error bars)')

    plt.legend(handles=[red_patch, blue_patch, ci_patch], loc='lower right', fontsize=10)
    plt.xlabel('Sample Index', fontsize=14)
    plt.ylabel('Predicted Probability', fontsize=14)
    # plt.title('Calibrated Probabilities with Error Bars (95% CI)', fontsize=15)
    plt.tight_layout()
    plt.show()


# %%

# Train using bootstrap aggregation
start_time = time.time()
mean_proba, std_proba, lower_ci, upper_ci = train_classifier_with_bootstrap_iso(X_train_val, y_train_val, X_test, y_test, n_bootstraps=100)
mean_proba_iso = mean_proba
iso_lower_ci = lower_ci
iso_upper_ci = upper_ci
iso_std_proba = std_proba
end_time = time.time()


start_time = time.time()
mean_proba, std_proba, lower_ci, upper_ci = train_classifier_with_bootstrap_platt(X_train_val, y_train_val, X_test, y_test, n_bootstraps=100)
mean_proba_platt = mean_proba
platt_lower_ci = lower_ci
platt_upper_ci = upper_ci
platt_std_proba = std_proba
end_time = time.time()


start_time = time.time()
mean_proba, std_proba, lower_ci, upper_ci = train_classifier_with_bootstrap(X_train_val, y_train_val, X_test, y_test, n_bootstraps=100)
mean_proba_no_cali = mean_proba
raw_lower_ci = lower_ci
raw_upper_ci = upper_ci
raw_std_proba = std_proba
end_time = time.time()

# %%
# Plott the calibration curves for each method
prob_true_platt, prob_pred_platt = plot_calibration_curve(y_test, mean_proba_platt)
prob_true_iso, prob_pred_iso = plot_calibration_curve(y_test, mean_proba_iso)
prob_true_none, prob_pred_none = plot_calibration_curve(y_test, mean_proba_no_cali)
plt.figure(figsize=(8, 6))
plt.plot(prob_pred_platt, prob_true_platt, marker='o', label='Platt Scaling')
plt.plot(prob_pred_iso, prob_true_iso, marker='o', label='Isotonic Regression')
plt.plot(prob_pred_none, prob_true_none, marker='o', label='Raw')
plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Perfect Calibration')

plt.xlabel('Bagged Mean Predicted Probability', fontsize=14)
plt.ylabel('True Fraction of Positives', fontsize=14)
plt.legend(fontsize=14)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
# %%

decision_threshold = 0.5
#choosing platt as mean_proba and cis
mean_proba = mean_proba_platt
lower_ci = platt_lower_ci
upper_ci = platt_upper_ci
std_proba = platt_std_proba
y_pred = (mean_proba > decision_threshold).astype(int)
confidence, entropy = get_confidence(mean_proba)

# Store results
actual_predicted_df = pd.DataFrame({
    'actual': y_test,
    'predicted': y_pred,
    'probability': mean_proba,
    'confidence': confidence,
    'entropy': entropy,
    'std_proba': std_proba,
    'lower_ci': lower_ci,
    'upper_ci': upper_ci
})

print(f"Training completed in {end_time - start_time:.2f} seconds")
print(actual_predicted_df)
import seaborn as sns

def prediction_conf(confidence_df):
    # Plot probability distributions for correct vs. incorrect predictions
    confidence_df["correct"] = confidence_df["actual"] == confidence_df["predicted"]

    # Split into correct and incorrect predictions
    correct_probs = confidence_df[confidence_df["correct"]]["confidence"]
    incorrect_probs = confidence_df[~confidence_df["correct"]]["confidence"]
    
    # ci = calc_conf_CI(confidence_df, confidence_level=0.95)['Pos Confidence CI']
    ci = confidence_df['upper_ci'] - confidence_df['lower_ci']



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

# %%
print(prediction_conf(actual_predicted_df))
# print(prediction_conf(probs_test))
from sklearn.metrics import classification_report, confusion_matrix
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
confident_predictions = actual_predicted_df[actual_predicted_df["confidence"] >= 0.9]

confusion_matrix_maker(confident_predictions, threshold=0.5)
# confusion_matrix_maker(confident_predictions, threshold=0.1)



# %%

# Get uncertain values
uncertain_values = actual_predicted_df[(actual_predicted_df['upper_ci'] > 0.55) & (actual_predicted_df['lower_ci'] < 0.45)]
uncertain_values = uncertain_values.reset_index(drop=True)

# Plot them
# plot_confidence_intervals(uncertain_values['actual'], uncertain_values['probability'], uncertain_values['lower_ci'], uncertain_values['upper_ci'])
# plot_confidence_intervals(y_test, mean_proba, lower_ci, upper_ci)

# Predicted class from mean probability
actual_predicted_df["predicted"] = (actual_predicted_df["probability"] >= 0.5).astype(int)

# Identify misclassified predictions
misclassified = actual_predicted_df["predicted"] != actual_predicted_df["actual"]

# Identify uncertain confidence intervals
uncertain = (actual_predicted_df['upper_ci'] > 0.8) & (actual_predicted_df['lower_ci'] < 0.2)

# Combine both
bad_predictions = actual_predicted_df[misclassified | uncertain].reset_index(drop=True)
plot_confidence_intervals(bad_predictions['actual'], bad_predictions['probability'], bad_predictions['lower_ci'], bad_predictions['upper_ci'])

bad_pred_index = list(bad_predictions['actual'].index)


# %%




# %%

# actual_predicted_df = uncertain_values # focus on uncertain values
uncertain_indices = list(uncertain_values.index)
uncertain_features = X_test.loc[uncertain_indices] # 139 rows
bad_pred_features = X_test.loc[bad_pred_index] # 27 rows


#Annotate with real and prob, remember to remove before shap
uncertain_features['actual'] = y_test.loc[uncertain_indices]
uncertain_features['prob'] = actual_predicted_df['probability'].loc[uncertain_indices]
bad_pred_features['actual'] = y_test.loc[bad_pred_index]
bad_pred_features['prob'] = actual_predicted_df['probability'].loc[bad_pred_index]

# # Now select the features which are taken from feature selection
selected_features_boruta = ['phyloP17way',
 'phyloP100way',
 'k24.Umap.MultiTrackMappability',
 'phastCons20way',
 'phyloP470way',
 'k50.Umap.MultiTrackMappability',
 'phastCons17way',
 'phastCons7way',
 '200GCContent',
 '200CpGCount',
 '1000CpGCount',
 '9_MGW',
 '17_HelT',
 '8_ProT',
 '9_ProT',
 '1_Shift_stiffness',
 '1_Stacking_energy_(RNA)',
 '2_Shift_shift',
 '3_Shift_slide',
 '4_Twist_slide',
 '4_Stacking_energy_(RNA)',
 'phyloP4way']

uncertain_features_selected = uncertain_features[selected_features_boruta]
bad_features_selected = bad_pred_features[selected_features_boruta]


# %%
# DO SHAP ANALYSIS OF THESE INDIVIDUAL PREDICTIONS
# First get selected features, the select those features from X_test, then look at which ones screw it up
# We need to see the differences between models, or do some sort of correlation test to see similarities between each variant to deduce why these are struggling.

import shap

#Train a new model using the selected features:

params = {
        "learning_rate": 0.2,
        "max_depth": 5,
        "n_estimators": 200,
        "reg_alpha": 0,
        "subsample": 0.5,
        "random_state": 30,
        'colsample_bytree': 0.9
    }
actual_predicted, feature_importance, final_model_tuple, results_tuple, data, chrom_performances = train_classifier(
            df, features=selected_features_boruta, classifier=XGBClassifier(**params), feature_importance=True
        )


xgb_model = final_model_tuple[0]



# %%
import shap
#uncertain
explainer = shap.Explainer(xgb_model)
shap_values_uncertain = explainer(uncertain_features_selected)
shap_values_bad = explainer(bad_features_selected)
shap_values_all = explainer(X_test[selected_features_boruta])

# shap.summary_plot(shap_values, uncertain_features_selected)
# shap.summary_plot(shap_values, bad_features_selected)
# shap.summary_plot(shap_values, X_test[selected_features_boruta])
# Local explanation for one prediction
shap.plots.waterfall(shap_values_bad[22]) # USE THIS FOR LOOKING AT LOW CONFIDENCE PREDICTIONS
# shap.plots.beeswarm(shap_values)

# %%

def summarise_shap_values(shap_values, feature_names):
    abs_vals = np.abs(shap_values.values)
    summary_df = pd.DataFrame({
        'feature': feature_names,
        'mean_abs_shap': abs_vals.mean(axis=0),
        'std_shap': abs_vals.std(axis=0)
    })
    return summary_df.sort_values(by='mean_abs_shap', ascending=False)

summary_uncertain = summarise_shap_values(shap_values_uncertain, selected_features_boruta)
summary_bad = summarise_shap_values(shap_values_bad, selected_features_boruta)
summary_all = summarise_shap_values(shap_values_all, selected_features_boruta)

shap_vals_all = np.abs(shap_values_all.values)
shap_vals_uncertain = np.abs(shap_values_uncertain.values)

from scipy.stats import mannwhitneyu

p_mannwhitney = []
for i in range(len(selected_features_boruta)):
    vals_all = np.abs(shap_values_all.values[:, i])
    vals_uncertain = np.abs(shap_values_uncertain.values[:, i])
    
    try:
        u_stat, p = mannwhitneyu(vals_all, vals_uncertain, alternative='two-sided')
    except:
        p = np.nan
    
    p_mannwhitney.append(p)

# Store results
stat_df = pd.DataFrame({
    'feature': selected_features_boruta,
    'p_mannwhitney_uncertain_vs_all': p_mannwhitney
}).sort_values(by='p_mannwhitney_uncertain_vs_all', na_position='last')

# %%
# --- Compute means and log2 fold-change ---
mean_shap_all = np.abs(shap_values_all.values).mean(axis=0)
mean_shap_uncertain = np.abs(shap_values_uncertain.values).mean(axis=0)

epsilon = 1e-8  # To avoid division by zero
log2_fc = np.log2((mean_shap_uncertain + epsilon) / (mean_shap_all + epsilon))

# --- Add to existing stat_df ---
stat_df['mean_shap_all'] = mean_shap_all
stat_df['mean_shap_uncertain'] = mean_shap_uncertain
stat_df['log2_fold_change'] = log2_fc



# %%

# Add class labels
pred_classes = (actual_predicted_df['probability'] > 0.5).astype(int)
drivers = X_test[pred_classes == 1][selected_features_boruta]
passengers = X_test[pred_classes == 0][selected_features_boruta]

# Compute SHAP values by class
shap_values_drivers = explainer(drivers)
shap_values_passengers = explainer(passengers)

# Compare average SHAP values by class
mean_driver_shap = np.abs(shap_values_drivers.values).mean(axis=0)
mean_passenger_shap = np.abs(shap_values_passengers.values).mean(axis=0)

shap_class_df = pd.DataFrame({
    'feature': selected_features_boruta,
    'mean_driver': mean_driver_shap,
    'mean_passenger': mean_passenger_shap,
    'driver > passenger': mean_driver_shap > mean_passenger_shap
})
import matplotlib.pyplot as plt

def plot_shap_comparison(summary1, summary2, label1="Uncertain", label2="All"):
    df = summary1.merge(summary2, on="feature", suffixes=(f"_{label1}", f"_{label2}"))
    df = df.sort_values(by=f"mean_abs_shap_{label1}", ascending=False)

    plt.figure(figsize=(10, 6))
    bar_width = 0.35
    index = np.arange(len(df))
    plt.bar(index, df[f"mean_abs_shap_{label1}"], bar_width, label=label1)
    plt.bar(index + bar_width, df[f"mean_abs_shap_{label2}"], bar_width, label=label2)
    plt.xticks(index + bar_width / 2, df['feature'], rotation=90)
    plt.ylabel("Mean |SHAP| Value")
    plt.title(f"SHAP Comparison: {label1} vs {label2}")
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_shap_comparison(summary_uncertain, summary_all, label1="Uncertain", label2="All")

# %%
# --------------------------------------------------
# Assume your data has been split, and you have:
# X_train, y_train, X_test, y_test, and a list of selected features.
selected_features_boruta = [ # Add/Remove pmed for changes in prediction
    'phyloP17way',
    'phyloP100way',
    'phastCons20way',
    'phyloP470way',
    'phastCons17way',
    'phastCons7way',
    'phyloP4way',
    'k24.Umap.MultiTrackMappability',
    'k50.Umap.MultiTrackMappability',
    '200GCContent',
    '200CpGCount',
    '1000CpGCount',
    '9_MGW',
    '17_HelT',
    '8_ProT',
    '9_ProT',
    '1_Shift_stiffness',
    '1_Stacking_energy_(RNA)',
    '2_Shift_shift',
    '3_Shift_slide',
    '4_Twist_slide',
    '4_Stacking_energy_(RNA)' 
]

# Restrict the data to the selected features
X_train_sel = X_train_val[selected_features_boruta]
X_test_sel = X_test[selected_features_boruta]

# %%



# %%

# --------------------------------------------------
# 1. Train the model using the selected features
params = {
    "learning_rate": 0.2,
    "max_depth": 5,
    "n_estimators": 200,
    "reg_alpha": 0,
    "subsample": 0.5,
    "random_state": 30,
    "colsample_bytree": 0.9
}
model = XGBClassifier(**params)
model.fit(X_train_sel, y_train_val)

# --------------------------------------------------
# 2. Compute SHAP Interaction Values Using TreeExplainer
# This returns an array of shape (n_samples, n_features, n_features)
tree_explainer = shap.TreeExplainer(model)
shap_interaction_values = tree_explainer.shap_interaction_values(bad_features_selected)

# --------------------------------------------------
# 3. Aggregate the Interaction Values
# Here we compute the mean absolute interaction value for each feature pair across all test samples.
avg_interaction = np.mean(np.abs(shap_interaction_values), axis=0)

# Create a DataFrame from the aggregated interaction values
interaction_df = pd.DataFrame(avg_interaction, 
                              index=bad_features_selected.columns, 
                              columns=bad_features_selected.columns)
print(interaction_df)

# --------------------------------------------------
# 4. Visualize the Aggregated Interaction Values with a Heatmap
import seaborn as sns
plt.figure(figsize=(17, 13))
sns.heatmap(interaction_df, annot=True, fmt=".3f", cmap="coolwarm")
# plt.title("Aggregated SHAP Interaction Values (Mean Absolute)")
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.show()

# %%

#  1. Split X_test_sel by class (assumes y_test is aligned with X_test_sel)
X_test_benign = X_test_sel[y_test == 0]
X_test_pathogenic = X_test_sel[y_test == 1]

# 2. Compute SHAP interaction values
explainer = shap.TreeExplainer(model)
shap_vals_benign = explainer.shap_interaction_values(X_test_benign)
shap_vals_pathogenic = explainer.shap_interaction_values(X_test_pathogenic)

# 3. Compute mean absolute interaction values
mean_interact_benign = np.mean(np.abs(shap_vals_benign), axis=0)
mean_interact_pathogenic = np.mean(np.abs(shap_vals_pathogenic), axis=0)

# 4. Compute the difference matrix (Pathogenic - Benign)
interaction_diff = mean_interact_pathogenic - mean_interact_benign
interaction_diff_df = pd.DataFrame(interaction_diff, 
                                   index=X_test_sel.columns, 
                                   columns=X_test_sel.columns)

# 5. Visualize difference heatmap
plt.figure(figsize=(17, 13))
sns.heatmap(interaction_diff_df, annot=True, fmt=".3f", cmap="coolwarm", center=0)
plt.title("Difference in SHAP Interaction Values (Pathogenic - Benign)")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

# 6. Extract top 10 most different feature interaction pairs
diff_list = []
features = X_test_sel.columns
import itertools
for i, j in itertools.combinations_with_replacement(range(len(features)), 2):
    f1, f2 = features[i], features[j]
    diff_value = interaction_diff[i, j]
    diff_list.append((f1, f2, diff_value, abs(diff_value)))

# Sort by absolute difference
diff_list_sorted = sorted(diff_list, key=lambda x: x[3], reverse=True)
top_10_diffs = diff_list_sorted[:10]

# 7. Print results
print("\nTop 10 Most Different SHAP Interaction Pairs (Pathogenic - Benign):\n")
for f1, f2, diff, abs_diff in top_10_diffs:
    direction = "↑ pathogenic" if diff > 0 else "↑ benign"
    print(f"{f1} × {f2}: Δ = {diff:.3f} ({direction})")

# %%


# Step 1: Create a list of tuples (f1, f2, diff, abs(diff))
features = list(X_test_sel.columns)
diff_list = []

for i, j in itertools.combinations_with_replacement(range(len(features)), 2):
    f1, f2 = features[i], features[j]
    diff_val = interaction_diff[i, j]
    diff_list.append((f1, f2, diff_val, abs(diff_val)))

# Step 2: Sort by absolute difference
diff_list_sorted = sorted(diff_list, key=lambda x: x[3], reverse=True)

# Step 3: Convert to DataFrame for display/export
top_n = 20  # Change as needed
top_diffs_df = pd.DataFrame(diff_list_sorted[:top_n], columns=["Feature 1", "Feature 2", "Difference", "Absolute Difference"])

# Optional: Add direction for interpretation
top_diffs_df["More Interaction In"] = top_diffs_df["Difference"].apply(lambda x: "Pathogenic" if x > 0 else "Benign")



# %%
import seaborn as sns
plt.figure(figsize=(10, 6))
sns.histplot(actual_predicted_df['probability'], bins=50, kde=True)
plt.title("Distribution of Predicted Probabilities")
plt.xlabel("Predicted Probability")
plt.ylabel("Frequency")
plt.show()

# 2. Scatter Plot with Confidence Intervals (Error Bars)
plt.figure(figsize=(12, 6))
# Create error bars for each sample: lower error and upper error
lower_error = actual_predicted_df['probability'] - actual_predicted_df['lower_ci']
upper_error = actual_predicted_df['upper_ci'] - actual_predicted_df['probability']
# Errorbar expects symmetric errors in a list: [lower_error, upper_error]
plt.errorbar(x=range(len(actual_predicted_df)), y=actual_predicted_df['probability'],
             yerr=[lower_error, upper_error],
             fmt='o', ecolor='gray', alpha=0.6, markersize=10)
plt.title("Predicted Probabilities with 95% Confidence Intervals")
plt.xlabel("Sample Index")
plt.ylabel("Predicted Probability")
plt.show()

# 3. Boxplot of Predicted Probabilities by Actual Class
plt.figure(figsize=(12, 6))
sns.boxplot(x='actual', y='probability', data=actual_predicted_df)
plt.title("Predicted Probability by Actual Class")
plt.xlabel("Actual Class")
plt.ylabel("Predicted Probability")
plt.show()

# %%
# 4. Violin plot for additional detail
import seaborn as sns
plt.figure(figsize=(8, 6))
ax = sns.violinplot(x='actual', y='confidence', data=bad_predictions, inner='quartile', linewidth=2) # or actual predicted
plt.xlabel("Actual Class", fontsize = 14)
plt.ylabel("Confidence", fontsize=14)

# Calculate and display mean confidence
mean_confidence_0 = bad_predictions[bad_predictions['actual'] == 0]['confidence'].mean()
mean_confidence_1 = bad_predictions[bad_predictions['actual'] == 1]['confidence'].mean()

# Annotate the plot with mean confidence values
# plt.text(0, actual_predicted_df['confidence'].max()*1.1, f"Mean: {mean_confidence_0:.2f}", horizontalalignment='center', size='medium', color='black', weight='semibold')
# plt.text(1, actual_predicted_df['confidence'].max()*1.1, f"Mean: {mean_confidence_1:.2f}", horizontalalignment='center', size='medium', color='black', weight='semibold')

plt.show()

# %%

plt.figure(figsize=(8, 6))
ax = sns.violinplot(x='actual', y='confidence', data=actual_predicted_df, inner='quartile', linewidth=2) # or actual predicted
plt.xlabel("Actual Class", fontsize = 14)
plt.ylabel("Confidence", fontsize=14)

# Calculate and display mean confidence
mean_confidence_0 = actual_predicted_df[actual_predicted_df['actual'] == 0]['confidence'].mean()
mean_confidence_1 = actual_predicted_df[actual_predicted_df['actual'] == 1]['confidence'].mean()

# Annotate the plot with mean confidence values
# plt.text(0, actual_predicted_df['confidence'].max()*1.1, f"Mean: {mean_confidence_0:.2f}", horizontalalignment='center', size='medium', color='black', weight='semibold')
# plt.text(1, actual_predicted_df['confidence'].max()*1.1, f"Mean: {mean_confidence_1:.2f}", horizontalalignment='center', size='medium', color='black', weight='semibold')

plt.show()

# %%
