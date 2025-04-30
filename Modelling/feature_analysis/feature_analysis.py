# %%
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/edatkinson/Repos/Modelling/new_approach')
from prepare_trainingdata import prepare_data
from train_classifier import train_classifier, process_df, alphaFold, dna_shape,gcContent,dna_props,cons_features
import itertools
from xgboost import XGBClassifier
from boruta import BorutaPy
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

# %%
df = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv', sep=',')
df = process_df(df)

# %%

feature_group_names = {
    'alphaFold': alphaFold,
    'dna_shape': dna_shape,
    'gcContent': gcContent,
    'cons_features': cons_features,
    'dna_props': dna_props
}

group_names = list(feature_group_names.keys())

# Initialize list to store results
results = []
feature_importance_total = defaultdict(float)
feature_counts = defaultdict(int)

# Generate all non-empty combinations of feature groups
for k in range(1, len(group_names) + 1):
    for combo_names in itertools.combinations(group_names, k):
        combined_features = [feature for name in combo_names for feature in feature_group_names[name]]
        
        actual_predicted, feature_importance, final_model_tuple, results_tuple, data, chromosome_performance  = train_classifier(
            df, features=combined_features, classifier=XGBClassifier(n_jobs=-1), feature_importance=True
        )

        mean_performance = results_tuple[0]['accuracy'].mean()
        print(combo_names,':',mean_performance)
        results.append((mean_performance, combo_names, feature_importance))

# %%

# %%
for _,_, res in results:
    for _, row in res.iterrows():
        feature = row['feature']
        importance = row['importance']
        try:
            importance = float(importance)
            feature_importance_total[feature] += importance
            feature_counts[feature] += 1
        except Exception as e:
            print(f"Skipping feature '{feature}' due to error: {e}")

       
# %%

average_importance = {
    feature: feature_importance_total[feature] / feature_counts[feature]
    for feature in feature_importance_total
}

for feature, importance in sorted(average_importance.items(), key=lambda x: x[1], reverse=True):
    print(f"{feature}: {importance:.4f}")

# %%

importance_df = pd.DataFrame.from_dict(average_importance, orient='index', columns=['importance'])
importance_df = importance_df.sort_values(by='importance', ascending=False)
importance_df.reset_index(inplace=True)
importance_df.rename(columns={'index': 'feature'}, inplace=True)


plt.figure(figsize=(12, 6))
sns.barplot(data=importance_df[importance_df['importance']>50], x='importance', y='feature', palette='viridis')
# plt.title("Average Feature Importance")
plt.xlabel("Importance", fontsize=15)
plt.ylabel("Feature", fontsize=15)
plt.tight_layout()
plt.show()



# %%

X_train_val, y_train_val, groups_train_val, X_test, y_test, groups_test = prepare_data(df, features=None)

# X_train_val.drop(columns=['genes', 'WTtrinuc', 'mutTrinuc'], inplace=True)

def select_features_with_boruta(X_train, y_train):
    """Performs feature selection using Boruta."""
   
    xgb = XGBClassifier(n_jobs=-1)
    X_train = X_train.replace(np.nan, 0)
    X_train = X_train.drop(columns=['pmed'])
    boruta = BorutaPy(xgb, n_estimators='auto', verbose=2, random_state=42)
    boruta.fit(X_train.values, y_train.values)
    
    selected_features = X_train.columns[boruta.support_]
    tentative_features = X_train.columns[boruta.support_weak_]
    return X_train[selected_features], selected_features.to_list(), tentative_features.to_list()

X_train_val_selected, selected_features, tentative_features = select_features_with_boruta(X_train_val, y_train_val)

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
actual_predicted, feature_importance, final_model_tuple, results_tuple = train_classifier(
            df, features=selected_features, classifier=XGBClassifier(**params), feature_importance=True
        )


# %%

#SHAP
import shap

xgb_model = final_model_tuple[0]

explainer = shap.Explainer(xgb_model)
shap_values = explainer(X_train_val_selected)
shap.summary_plot(shap_values, X_train_val_selected)

# Local explanation for one prediction
shap.plots.waterfall(shap_values[0]) # USE THIS FOR LOOKING AT LOW CONFIDENCE PREDICTIONS


# %%
features = importance_df[importance_df['importance']>50]['feature']

# %%

xgb_model = final_model_tuple[0].fit(X_train_val[features.values], y_train_val)

explainer = shap.Explainer(xgb_model)
shap_values = explainer(X_train_val[features.values])
shap.summary_plot(shap_values, X_train_val[features.values])

# Local explanation for one prediction
shap.plots.waterfall(shap_values[0]) # USE THIS FOR LOOKING AT LOW CONFIDENCE PREDICTIONS



# %%
