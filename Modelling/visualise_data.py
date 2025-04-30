
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from models.feature_selection import balance_data
from models.cosmic import alphaFold, cons_features, gcContent

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
# %%

# Load the dataset
df = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv', sep=',')
df2 = pd.read_csv('/Users/edatkinson/Repos/Modelling/data/pmed_log_frameshift.csv', sep=',')

#Count how many drivers and passengers there are:
# Balance the data as per your function
# df = balance_data(df)
# df2 = balance_data(df2)

print(df2['phyloP4way'].isnull().sum())




# feat_list = ['phyloP20way', 'k36.Bismap.MultiTrackMappability', 'phyloP17way', 'phastCons470way', 
#              'k24.Umap.MultiTrackMappability', 'k100.Umap.MultiTrackMappability', 'phastCons17way', 
#              'phyloP4way', 'phastCons4way', '14_HelT', '1_Major_Groove_Width', '1_Tilt_rise', 
#              '1_Entropy_(RNA)', '2_Tilt_(RNA)', 'isotropic_temperature']

# feat_list = alphaFold 


feat_list = ['phyloP17way',
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


# ["phyloP4way","phyloP7way", "phyloP17way", "phyloP20way", "phyloP30way",
        # "phyloP100way","phyloP470way","phastCons7way",'phastCons4way', "phastCons17way",
        # "phastCons20way", "phastCons30way", "phastCons100way", "phastCons470way"]
        # "k24.Bismap.MultiTrackMappability", "k36.Umap.MultiTrackMappability",
        # "k36.Bismap.MultiTrackMappability", "k24.Umap.MultiTrackMappability",
        # "k50.Bismap.MultiTrackMappability","k50.Umap.MultiTrackMappability",
        # "k100.Bismap.MultiTrackMappability", "k100.Umap.MultiTrackMappability"]

# feat_list = ['80GCContent', '200GCContent', '200CpGCount', '500GCContent', '1000CpGCount', '40CpG_obs_exp']

# Filter the relevant features and add the 'driver_stat' column
df = df[feat_list + ['driver_stat']]


# Plotting histograms for each feature
plt.figure(figsize=(15, 10))
for i, feature in enumerate(feat_list):
    plt.subplot(5, 5, i + 1)
    sns.histplot(df[df['driver_stat'] == 1][feature], color='blue', kde=True, label='Driver', bins=30)
    sns.histplot(df[df['driver_stat'] == 0][feature], color='red', kde=True, label='Passenger', bins=30)
    plt.title(feature)
    plt.legend()

plt.tight_layout()


# Plotting histograms for each feature
plt.figure(figsize=(15, 10))
for i, feature in enumerate(feat_list):
    plt.subplot(5, 5, i + 1)
    sns.histplot(df2[df2['driver_stat'] == 1][feature], color='blue', kde=True, label='Driver', bins=30)
    sns.histplot(df2[df2['driver_stat'] == 0][feature], color='red', kde=True, label='Passenger', bins=30)
    plt.title(feature)
    plt.legend()

plt.tight_layout()
plt.show()



plt.figure(figsize=(15, 10))
for i, feature in enumerate(feat_list):
    plt.subplot(5, 5, i + 1)
    sns.boxplot(x='driver_stat', y=feature, data=df, hue='driver_stat', palette={0: 'red', 1: 'blue'})
    plt.title(feature)

plt.tight_layout()
plt.show()

# Violin plots for each feature
plt.figure(figsize=(15, 10))
for i, feature in enumerate(feat_list):
    plt.subplot(5, 5, i + 1)
    sns.violinplot(x='driver_stat', y=feature, data=df, hue='driver_stat',palette={0: 'red', 1: 'blue'})
    plt.title(feature)

plt.tight_layout()
plt.show()


# from scipy.stats import ttest_ind

# # Perform t-tests for each feature
# for feature in feat_list:
#     driver_values = df[df['driver_stat'] == 1][feature]
#     passenger_values = df[df['driver_stat'] == 0][feature]
#     stat, p_value = ttest_ind(driver_values, passenger_values, nan_policy='omit')
    
#     print(f"Feature: {feature} - p-value: {p_value}")


# for feature in feat_list:
#     plt.figure(figsize=(6,4))
#     sns.boxplot(x=data['driver_stat'], y=data[feature])
#     plt.title(f'{feature} distribution in Driver vs Passenger')
#     plt.show()

# # Visualise the proportions of the dataset using charts:

# # Total Size, proportions or driver to passenger and ratio of chromosomes in each class.


# for col in data[feat_list]:
#     sns.histplot(data[col], kde=True)
#     plt.title(f"Distribution of {col}")
#     plt.show()


# import seaborn as sns
# import matplotlib.pyplot as plt

# # Compute correlation matrix
# corr_matrix = data[feat_list].corr()

# # Visualize correlation matrix
# plt.figure(figsize=(12, 8))
# sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt='.2f')
# plt.title("Correlation Matrix")
# plt.show()

# # Check for high correlations (e.g., > 0.8 or < -0.8)
# high_corr_features = [(col, row) for col in corr_matrix.columns for row in corr_matrix.index if abs(corr_matrix.loc[row, col]) > 0.8 and col != row]
# print("Highly Correlated Features:", high_corr_features)






# # Drop non-numeric columns and target
# X = df.drop(columns=['chrom', 'pos', 'ref_allele', 'alt_allele', 'driver_stat', 'merged', 'WTtrinuc', 'mutTrinuc'])  # Adjust according to your data
# y = df['driver_stat']

# # Handling missing values (if any)
# X = X.fillna(X.mean())  # Fill with the mean, you can also try other methods

# # Standardize the features (important for PCA)
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(X)

# # Apply PCA
# pca = PCA()
# X_pca = pca.fit_transform(X_scaled)

# # Explained Variance Ratio - This shows how much information each principal component holds
# explained_variance_ratio = pca.explained_variance_ratio_

# # Visualize the explained variance
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, marker='o')
# plt.title("Explained Variance by Each Principal Component")
# plt.xlabel("Principal Components")
# plt.ylabel("Explained Variance Ratio")
# plt.grid(True)
# plt.show()

# # Cumulative explained variance
# cumulative_explained_variance = np.cumsum(explained_variance_ratio)
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(cumulative_explained_variance) + 1), cumulative_explained_variance, marker='o', color='r')
# plt.title("Cumulative Explained Variance by Principal Components")
# plt.xlabel("Principal Components")
# plt.ylabel("Cumulative Explained Variance")
# plt.grid(True)
# plt.show()


# # Project the data onto the first 2 principal components for visualization
# X_pca_2d = X_pca[:, :2]

# # Create a DataFrame for visualization
# pca_df = pd.DataFrame(data=X_pca_2d, columns=['PC1', 'PC2'])
# pca_df['driver_stat'] = y  # Add target variable to color points

# # Visualize the projection onto the first two PCs
# plt.figure(figsize=(10, 6))
# sns.scatterplot(x='PC1', y='PC2', hue='driver_stat', data=pca_df, palette='viridis', alpha=0.7)
# plt.title("PCA: First Two Principal Components")
# plt.xlabel("Principal Component 1")
# plt.ylabel("Principal Component 2")
# plt.grid(True)
# plt.show()

# # Number of components to retain for 90% explained variance
# n_components = np.argmax(cumulative_explained_variance >= 0.90) + 1
# print(f"Number of components explaining at least 90% variance: {n_components}")

# # Apply PCA again, keeping the chosen number of components
# pca = PCA(n_components=n_components)
# X_pca_reduced = pca.fit_transform(X_scaled)
# print(X_pca_reduced.shape)

# Now you can use X_pca_reduced as input features for machine learning models


# feat_list = ['phyloP17way', 'phyloP100way', 'k24.Umap.MultiTrackMappability', 'phyloP7way', 'phyloP470way', 'phastCons17way', '80GCContent', '200GCContent', '200CpGCount', '500GCContent', '1000CpGCount', '9_MGW', '14_HelT', '9_ProT', '1_Shift_shift', '1_Twist_shift', '1_Tilt_(RNA)', '2_Shift_shift', '2_Shift_slide', '2_Tilt_(DNA-protein_complex)', '3_Shift_shift', '4_Roll_roll', '4_Shift_shift', '4_Twist_shift', '4_Stacking_energy_(RNA)', 'atom_site_occupancy', 'phastCons4way', 'phyloP4way']

# %%
