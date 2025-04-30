# %%
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/edatkinson/Repos/Modelling/new_approach')
from prepare_trainingdata import prepare_data
from train_classifier import train_classifier, process_df, alphaFold, dna_shape,gcContent,dna_props,cons_features
import itertools
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score

import matplotlib.pyplot as plt

# %%
df = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv', sep=',')


# %%

def analyse_cons(cons_features, df, n_repeats=5):
    results = {
        'phyloP': [],
        'phastCons': []
    }
    phylo_phast = [con for con in cons_features if 'phyloP' in con or 'phastCons' in con]
    phylo_phast = sorted(phylo_phast)
    for con in phylo_phast:
        cv_scores = []
        test_scores = []
        
        for _ in range(n_repeats):
            df_copy = df.copy(deep=True)  # Ensure df is fresh for each iteration
            df_processed = process_df(df_copy) 
            actual_predicted, feature_importance, final_model_tuple, results_tuple, data, chrom_performances = train_classifier(
                df_processed, features=[con], classifier=None, feature_importance=True
            )

            cvperformance = accuracy_score(actual_predicted[0][0][0], actual_predicted[0][0][1])
            tperformance = accuracy_score(actual_predicted[0][1][0], actual_predicted[0][1][1])

            cv_scores.append(cvperformance)
            test_scores.append(tperformance)

        
        avg_cvperformance = np.mean(cv_scores)
        avg_tperformance = np.mean(test_scores)
        if 'phyloP' in con:
            results['phyloP'].append((avg_cvperformance, avg_tperformance, con))
        else:
            results['phastCons'].append((avg_cvperformance, avg_tperformance, con))

        print(f"{con}: Test Accuracy (Avg over {n_repeats} runs) = {avg_tperformance:.4f}")

    return results


res = analyse_cons(cons_features, df, n_repeats=10)  # Increase the number of repeats for more stability


# %%
res_df = pd.DataFrame(res['phyloP'], columns=['cv', 'test', 'num'])
res_df_cons = pd.DataFrame(res['phastCons'], columns=['cv', 'test', 'num'])
# Strip the num column to numbers
res_df['num'] = res_df['num'].str.extract(r'(\d+)')
res_df_cons['num'] = res_df_cons['num'].str.extract(r'(\d+)')


# Convert 'num' column to numeric
res_df['num'] = pd.to_numeric(res_df['num'])
res_df_cons['num'] = pd.to_numeric(res_df_cons['num'])

# Sort the DataFrame by 'num' column
res_df = res_df.sort_values(by='num')
res_df_cons = res_df_cons.sort_values(by='num')
# Convert 'num' column back to string for non-numeric x-axis
res_df['num'] = res_df['num'].astype(str)
res_df_cons['num'] = res_df_cons['num'].astype(str)
# Plot bar chart
plt.figure(figsize=(10, 6))
plt.plot(res_df['num'], res_df['test'], color='blue', alpha=0.7, label='PhyloP', linewidth=2, marker='o')
plt.plot(res_df_cons['num'], res_df_cons['test'], color='red', alpha=0.7, label='PhastCons', linewidth=2, marker='o')
plt.xlabel('N-way', fontsize=15)
plt.ylabel('Accuracy', fontsize=15)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(fontsize = 14)
plt.show()

# plt.figure(figsize=(10, 6))
# plt.plot(res_df['num'], res_df['cv'], color='blue', alpha=0.7, label='phyloP')
# plt.plot(res_df_cons['num'], res_df_cons['cv'], color='red', alpha=0.7, label='phastCons')
# plt.xlabel('Feature Number')
# plt.ylabel('Accuracy')
# plt.legend()
# plt.show()

# %%
