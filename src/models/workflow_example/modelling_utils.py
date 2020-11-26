import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import cohen_kappa_score, roc_auc_score
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE, RandomOverSampler


def getMatchingColNames(operators, features):
    col_names = []
    for col in features.columns:
        if any(operator in col for operator in operators):
            col_names.append(col)
    return col_names

# drop columns with zero variance
def dropZeroVarianceCols(data):
    if not data.empty:
        var_df = data.var()
        keep_col = []
        for col in var_df.index:
            if var_df.loc[col] > 0:
                keep_col.append(col)
        data_drop_cols_var = data.loc[:, keep_col]
    else:
        data_drop_cols_var = data
    return data_drop_cols_var

# normalize based on all participants: return fitted scaler
def getNormAllParticipantsScaler(features, scaler_flag):
    # MinMaxScaler
    if scaler_flag == "minmaxscaler":
        scaler = MinMaxScaler()
    # StandardScaler
    elif scaler_flag == "standardscaler":
        scaler = StandardScaler()
    # RobustScaler
    elif scaler_flag == "robustscaler":
        scaler = RobustScaler()
    else:
        # throw exception
        raise ValueError("The normalization method is not predefined, please check if the PARAMS_FOR_ANALYSIS.NORMALIZED in config.yaml file is correct.")
    scaler.fit(features)
    return scaler

# get metrics: accuracy, precision0, recall0, f10, precision1, recall1, f11, f1_macro, auc, kappa
def getMetrics(pred_y, pred_y_proba, true_y):
    metrics = {}
    count = len(np.unique(true_y))
    label= np.unique(true_y)[0]
    # metrics for all categories
    metrics["accuracy"] = accuracy_score(true_y, pred_y)
    metrics["f1_macro"] = f1_score(true_y, pred_y, average="macro") # unweighted mean
    metrics["auc"] = np.nan if count == 1 else roc_auc_score(true_y, pred_y_proba)
    metrics["kappa"] = cohen_kappa_score(true_y, pred_y)
    # metrics for label 0
    metrics["precision0"] = np.nan if (count == 1 and label == 1) else precision_score(true_y, pred_y, average=None, labels=[0,1], zero_division=0)[0]
    metrics["recall0"] = np.nan if (count == 1 and label == 1) else recall_score(true_y, pred_y, average=None, labels=[0,1])[0]
    metrics["f10"] = np.nan if (count == 1 and label == 1) else f1_score(true_y, pred_y, average=None, labels=[0,1])[0]
    # metrics for label 1
    metrics["precision1"] = np.nan if (count == 1 and label == 0) else precision_score(true_y, pred_y, average=None, labels=[0,1], zero_division=0)[1]
    metrics["recall1"] = np.nan if (count == 1 and label == 0) else recall_score(true_y, pred_y, average=None, labels=[0,1])[1]
    metrics["f11"] = np.nan if (count == 1 and label == 0) else f1_score(true_y, pred_y, average=None, labels=[0,1])[1]

    return metrics

# get feature importances
def getFeatureImportances(model, clf, cols):
    if model == "LogReg":
        # Extract the coefficient of the features in the decision function
        # Calculate the absolute value
        # Normalize it to sum 1
        feature_importances = pd.DataFrame(zip(clf.coef_[0],cols), columns=["Value", "Feature"])
        feature_importances["Value"] = feature_importances["Value"].abs()/feature_importances["Value"].abs().sum()
    elif model == "kNN":
        # Feature importance is not defined for the KNN Classification, return an empty dataframe
        feature_importances = pd.DataFrame(columns=["Value", "Feature"])
    elif model == "SVM":
        # Coefficient of the features are only available for linear kernel
        try:
            # For linear kernel
            # Extract the coefficient of the features in the decision function
            # Calculate the absolute value
            # Normalize it to sum 1
            feature_importances = pd.DataFrame(zip(clf.coef_[0],cols), columns=["Value", "Feature"])
            feature_importances["Value"] = feature_importances["Value"].abs()/feature_importances["Value"].abs().sum()
        except:
            # For nonlinear kernel, return an empty dataframe directly
            feature_importances = pd.DataFrame(columns=["Value", "Feature"])
    elif model == "LightGBM":
        # Extract feature_importances_ and normalize it to sum 1
        feature_importances = pd.DataFrame(zip(clf.feature_importances_,cols), columns=["Value", "Feature"])
        feature_importances["Value"] = feature_importances["Value"]/feature_importances["Value"].sum()
    else:
        # For DT, RF, GB, XGBoost classifier, extract feature_importances_. This field has already been normalized.
        feature_importances = pd.DataFrame(zip(clf.feature_importances_,cols), columns=["Value", "Feature"])

    feature_importances = feature_importances.set_index(["Feature"]).T

    return feature_importances

def createPipeline(model, oversampler_type):

    if oversampler_type == "SMOTE":
        oversampler = SMOTE(sampling_strategy="minority", random_state=0)
    elif oversampler_type == "RandomOverSampler":
        oversampler = RandomOverSampler(sampling_strategy="minority", random_state=0)
    else:
        raise ValueError("RAPIDS pipeline only support 'SMOTE' and 'RandomOverSampler' oversampling methods.")

    if model == "LogReg":
        from sklearn.linear_model import LogisticRegression
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", LogisticRegression(random_state=0))
        ])
    elif model == "kNN":
        from sklearn.neighbors import KNeighborsClassifier
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", KNeighborsClassifier())
        ])
    elif model == "SVM":
        from sklearn.svm import SVC
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", SVC(random_state=0, probability=True))
        ])
    elif model == "DT":
        from sklearn.tree import DecisionTreeClassifier
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", DecisionTreeClassifier(random_state=0))
        ])
    elif model == "RF":
        from sklearn.ensemble import RandomForestClassifier
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", RandomForestClassifier(random_state=0))
        ])
    elif model == "GB":
        from sklearn.ensemble import GradientBoostingClassifier
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", GradientBoostingClassifier(random_state=0))
        ])
    elif model == "XGBoost":
        from xgboost import XGBClassifier
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", XGBClassifier(random_state=0, n_jobs=36))
        ])
    elif model == "LightGBM":
        from lightgbm import LGBMClassifier
        pipeline = Pipeline([
            ("sampling", oversampler),
            ("clf", LGBMClassifier(random_state=0, n_jobs=36))
        ])
    else:
        raise ValueError("RAPIDS pipeline only support LogReg, kNN, SVM, DT, RF, GB, XGBoost, and LightGBM algorithms for classification problems.")

    return pipeline
