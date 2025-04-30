import pandas as pd
import json
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report, accuracy_score, f1_score, precision_score, recall_score
from imblearn.under_sampling import RandomUnderSampler

# Define paths
processed_data_path = "/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv"
output_folder = "/Users/edatkinson/Repos/Modelling/new_approach/output/"
params_output_path = f"{output_folder}best_xgb_params.json"
metrics_output_path = f"{output_folder}xgb_metrics.txt"

# Load the processed DataFrame
df = pd.read_csv(processed_data_path)
print("Loaded processed DataFrame:")
print(df.head())

# Define features and target
features = df.drop(columns=['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'WTtrinuc', 'mutTrinuc', 'genes'])

target = df['driver_stat']

# Define a parameter grid for XGBoost
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [3, 5, 7],
    'learning_rate': [0.01, 0.1, 0.2],
    'subsample': [0.6, 0.8, 1.0],
    'colsample_bytree': [0.6, 0.8, 1.0],
    # 'gamma': [0, 0.1, 0.2],
    'reg_alpha': [0, 0.1, 1],
    'reg_lambda': [1, 1.5, 2],
}

# Initialise the XGBoost classifier
xgb = XGBClassifier(random_state=42)

# Initialise GridSearchCV
grid_search = GridSearchCV(
    estimator=xgb,
    param_grid=param_grid,
    scoring='f1',  # You can change this to 'accuracy', 'roc_auc', etc.
    cv=3,          # 3-fold cross-validation
    verbose=1,
    n_jobs=-1      # Use all available cores
)

# Perform hyperparameter optimisation
print("Starting hyperparameter optimisation for XGBoost...")
sampler = RandomUnderSampler(random_state=42)
features_undersampled, target_undersampled = sampler.fit_resample(features, target)

print(features_undersampled.shape, target_undersampled.shape)
grid_search.fit(features_undersampled, target_undersampled)

# Get the best parameters and the corresponding score
best_params = grid_search.best_params_
best_score = grid_search.best_score_

print(f"Best Parameters: {best_params}")
print(f"Best Cross-Validation F1 Score: {best_score:.4f}")

# Save best parameters to a JSON file
with open(params_output_path, 'w') as f:
    json.dump(best_params, f, indent=4)
print(f"Best parameters saved to {params_output_path}")

# Evaluate the best model on the entire dataset
best_xgb = grid_search.best_estimator_
predictions = best_xgb.predict(features)

# Compute evaluation metrics
accuracy = accuracy_score(target, predictions)
precision = precision_score(target, predictions)
recall = recall_score(target, predictions)
f1 = f1_score(target, predictions)
classification_rep = classification_report(target, predictions)

# Print evaluation metrics
print("\nEvaluation Metrics with Best Parameters:")
print(f"Accuracy: {accuracy:.4f}")
print(f"Precision: {precision:.4f}")
print(f"Recall: {recall:.4f}")
print(f"F1 Score: {f1:.4f}")
print("\nClassification Report:")
print(classification_rep)

# Save metrics to a text file
with open(metrics_output_path, 'w') as f:
    f.write("Evaluation Metrics with Best Parameters:\n")
    f.write(f"Accuracy: {accuracy:.4f}\n")
    f.write(f"Precision: {precision:.4f}\n")
    f.write(f"Recall: {recall:.4f}\n")
    f.write(f"F1 Score: {f1:.4f}\n\n")
    f.write("Classification Report:\n")
    f.write(classification_rep)
print(f"Metrics saved to {metrics_output_path}")