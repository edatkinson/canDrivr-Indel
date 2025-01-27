
import xgboost as xgb
from xgboost import XGBClassifier
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold, LeaveOneGroupOut, cross_val_score
from sklearn.metrics import accuracy_score, classification_report


def prepare_data(data, features=None):
    data = data.drop(columns=['pos','ref_allele','alt_allele', 'WTtrinuc', 'mutTrinuc','mutant_AA','WT_AA'])
    data = data.replace(np.nan, 0)
    groups = data['chrom'] #For LOCO CV
    # groups = np.unique(groups)
    X = data.drop(columns=['driver_stat', 'chrom'])
    y = data['driver_stat']
    X.columns = range(X.shape[1])

    return X, y, groups





# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from xgboost import XGBClassifier
# from sklearn.model_selection import train_test_split, cross_val_score, LeaveOneGroupOut
# from sklearn.feature_selection import SelectFromModel
# from sklearn.preprocessing import StandardScaler, PolynomialFeatures
# from sklearn.metrics import accuracy_score, classification_report
# from scipy.stats import boxcox



# if __name__ == "__main__":

#     data = pd.read_csv('/Users/edatkinson/Repos/canDrivr-Indel/canDrivr/Annotated_data.csv',sep=',') 
#     X, y, groups = prepare_data(data)
    
#     # Split data
#     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)

#     ### --- 2. Data Transformations ---
#     # Standardization (for non-tree models)
#     scaler = StandardScaler()
#     X_train_scaled = scaler.fit_transform(X_train)
#     X_test_scaled = scaler.transform(X_test)

#     ### --- 3. Feature Selection ---
#     # Train initial XGBoost model
#     xgb_model = XGBClassifier(objective='binary:logistic', eval_metric='logloss')
#     xgb_model.fit(X_train_scaled, y_train)

#     # Select important features
#     selector = SelectFromModel(xgb_model, threshold='mean', prefit=True)
#     X_train_selected = selector.transform(X_train_scaled)
#     X_test_selected = selector.transform(X_test_scaled)

#     print(f"Original feature count: {X_train.shape[1]}, Selected feature count: {X_train_selected.shape[1]}")

#     # Retrain with selected features
#     xgb_model_selected = XGBClassifier(objective='binary:logistic', eval_metric='logloss')
#     xgb_model_selected.fit(X_train_selected, y_train)

#     ### --- 4. Evaluation ---
#     y_pred = xgb_model_selected.predict(X_test_selected)
#     y_test_pred_proba = xgb_model_selected.predict_proba(X_test_selected)[:, 1]

#     print("Test Set Accuracy after feature selection and transformations:", accuracy_score(y_test, y_pred))
#     print("Classification Report:\n", classification_report(y_test, y_pred))

#     ## --- 5. Leave-One-Group-Out Cross-Validation (LOCO-CV) ---
#     logo = LeaveOneGroupOut()
#     scores = cross_val_score(xgb_model_selected, (X_train), y_train, cv=logo, groups=groups, scoring='accuracy')

#     print(f"LOGO-CV accuracy scores after feature selection and transformations: {scores}")
#     print(f"Mean LOGO-CV accuracy: {np.mean(scores)}")




