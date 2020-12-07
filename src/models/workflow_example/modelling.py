import pandas as pd
import numpy as np
from modelling_utils import getMatchingColNames, getNormAllParticipantsScaler, getMetrics, getFeatureImportances, createPipeline
from sklearn.model_selection import LeaveOneOut, GridSearchCV



def preprocessNumericalFeatures(train_numerical_features, test_numerical_features, scaler, flag):
    # fillna with mean
    if flag == "train":
        numerical_features = train_numerical_features.fillna(train_numerical_features.mean())
    elif flag == "test":
        numerical_features = test_numerical_features.fillna(train_numerical_features.mean())
    else:
        raise ValueError("flag should be 'train' or 'test'")
    # normalize
    if scaler != "notnormalized":
        scaler = getNormAllParticipantsScaler(train_numerical_features, scaler)
        numerical_features = pd.DataFrame(scaler.transform(numerical_features), index=numerical_features.index, columns=numerical_features.columns)

    return numerical_features

def preprocessCategoricalFeatures(categorical_features, mode_categorical_features):
    # fillna with mode
    categorical_features = categorical_features.fillna(mode_categorical_features)
    # one-hot encoding
    categorical_features = categorical_features.apply(lambda col: col.astype("category"))
    if not categorical_features.empty:
        categorical_features = pd.get_dummies(categorical_features)
    return categorical_features

def splitNumericalCategoricalFeatures(features, categorical_feature_colnames):
    numerical_features = features.drop(categorical_feature_colnames, axis=1)
    categorical_features = features[categorical_feature_colnames].copy()
    return numerical_features, categorical_features

def preprocesFeatures(train_numerical_features, test_numerical_features, categorical_features, mode_categorical_features, scaler, flag):
    numerical_features = preprocessNumericalFeatures(train_numerical_features, test_numerical_features, scaler, flag)
    categorical_features = preprocessCategoricalFeatures(categorical_features, mode_categorical_features)
    features = pd.concat([numerical_features, categorical_features], axis=1)
    return features


##############################################################
# Summary of the workflow
# Step 1. Read parameters and data
# Step 2. Nested cross validation
# Step 3. Model evaluation
# Step 4. Save results, parameters, and metrics to CSV files
##############################################################

# For reproducibility
np.random.seed(0)

# Step 1. Read parameters and data
# Read parameters
model = snakemake.params["model"]
scaler = snakemake.params["scaler"]
cv_method = snakemake.params["cv_method"]
categorical_operators = snakemake.params["categorical_operators"]
categorical_colnames_demographic_features = snakemake.params["categorical_demographic_features"]
model_hyperparams = snakemake.params["model_hyperparams"][model]


# Read data and split
data = pd.read_csv(snakemake.input["data"])
index_columns = ["local_segment", "local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"]
if "pid" in data.columns:
    index_columns.append("pid")
data.set_index(index_columns, inplace=True)

data_x, data_y = data.drop("target", axis=1), data[["target"]]

if "pid" in index_columns:
    categorical_feature_colnames = categorical_colnames_demographic_features + getMatchingColNames(categorical_operators, data_x)
else:
    categorical_feature_colnames =  getMatchingColNames(categorical_operators, data_x)



# Step 2. Nested cross validation
cv_class = globals()[cv_method]
inner_cv = cv_class()
outer_cv = cv_class()

fold_id, pid, best_params, true_y, pred_y, pred_y_proba = [], [], [], [], [], []
feature_importances_all_folds = pd.DataFrame()
fold_count = 1

# Outer cross validation
for train_index, test_index in outer_cv.split(data_x):

    # Split train and test, numerical and categorical features
    train_x, test_x = data_x.iloc[train_index], data_x.iloc[test_index]
    train_numerical_features, train_categorical_features = splitNumericalCategoricalFeatures(train_x, categorical_feature_colnames)
    train_y, test_y = data_y.iloc[train_index], data_y.iloc[test_index]
    test_numerical_features, test_categorical_features = splitNumericalCategoricalFeatures(test_x, categorical_feature_colnames)

    # Preprocess: impute and normalize
    mode_categorical_features = train_categorical_features.mode().iloc[0]
    train_x = preprocesFeatures(train_numerical_features, None, train_categorical_features, mode_categorical_features, scaler, "train")
    test_x = preprocesFeatures(train_numerical_features, test_numerical_features, test_categorical_features, mode_categorical_features, scaler, "test")
    train_x, test_x = train_x.align(test_x, join="outer", axis=1, fill_value=0) # in case we get rid off categorical columns

    # Compute number of participants and features
    # values do not change between folds
    if fold_count == 1:
        num_of_rows = train_x.shape[0] + test_x.shape[0]
        num_of_features = train_x.shape[1]

    targets_value_counts = train_y["target"].value_counts()
    if len(targets_value_counts) < 2 or max(targets_value_counts) < 5:
        notes = open(snakemake.log[0], mode="w")
        notes.write(targets_value_counts.to_string())
        notes.close()
        break

    # Inner cross validation
    if min(targets_value_counts) >= 6:
        # SMOTE requires n_neighbors <= n_samples, the default value of n_neighbors is 6
        clf = GridSearchCV(estimator=createPipeline(model, "SMOTE"), param_grid=model_hyperparams, cv=inner_cv, scoring="f1_macro")
    else:
        # RandomOverSampler: over-sample the minority class(es) by picking samples at random with replacement.
        clf = GridSearchCV(estimator=createPipeline(model, "RandomOverSampler"), param_grid=model_hyperparams, cv=inner_cv, scoring="f1_macro")
    clf.fit(train_x, train_y.values.ravel())

    # Collect results and parameters
    best_params = best_params + [clf.best_params_]
    cur_fold_pred = clf.predict(test_x).tolist()
    pred_y = pred_y + cur_fold_pred

    proba_of_two_categories = clf.predict_proba(test_x).tolist()
    pred_y_proba = pred_y_proba + [probabilities[clf.classes_.tolist().index(1)] for probabilities in proba_of_two_categories]

    true_y = true_y + test_y.values.ravel().tolist()
    pid = pid + test_y.index.tolist() # each test partition (fold) in the outer cv is a participant (LeaveOneOut cv)
    feature_importances_current_fold = getFeatureImportances(model, clf.best_estimator_.steps[1][1], train_x.columns)
    feature_importances_all_folds = pd.concat([feature_importances_all_folds, feature_importances_current_fold], sort=False, axis=0)
    fold_id.append(fold_count)
    fold_count = fold_count + 1

# Step 3. Model evaluation
if len(pred_y) > 1:
    metrics = getMetrics(pred_y, pred_y_proba, true_y)
else:
    metrics = {"accuracy": None, "precision0": None, "recall0": None, "f10": None, "precision1": None, "recall1": None, "f11": None, "f1_macro": None, "auc": None, "kappa": None}

# Step 4. Save results, parameters, and metrics to CSV files
fold_predictions = pd.DataFrame({"fold_id": fold_id, "pid": pid, "hyperparameters": best_params, "true_y": true_y, "pred_y": pred_y, "pred_y_proba": pred_y_proba})
fold_metrics = pd.DataFrame({"fold_id":[], "accuracy":[], "precision0": [], "recall0": [], "f10": [], "precision1": [], "recall1": [], "f11": [], "f1_macro": [], "auc": [], "kappa": []})
overall_results = pd.DataFrame({"num_of_rows": [num_of_rows], "num_of_features": [num_of_features], "model": [model], "cv_method": [cv_method], "scaler": [scaler], "accuracy": [metrics["accuracy"]], "precision0": [metrics["precision0"]], "recall0": [metrics["recall0"]], "f10": [metrics["f10"]], "precision1": [metrics["precision1"]], "recall1": [metrics["recall1"]], "f11": [metrics["f11"]], "f1_macro": [metrics["f1_macro"]], "auc": [metrics["auc"]], "kappa": [metrics["kappa"]]})
feature_importances_all_folds.insert(loc=0, column="fold_id", value=fold_id)
feature_importances_all_folds.insert(loc=1, column="pid", value=pid)

fold_predictions.to_csv(snakemake.output["fold_predictions"], index=False)
fold_metrics.to_csv(snakemake.output["fold_metrics"], index=False)
overall_results.to_csv(snakemake.output["overall_results"], index=False)
feature_importances_all_folds.to_csv(snakemake.output["fold_feature_importances"], index=False)
