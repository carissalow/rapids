import pandas as pd
from modeling_utils import dropZeroVarianceCols, getNormAllParticipantsScaler, getMetrics, getFeatureImportances, createPipeline
from sklearn.model_selection import train_test_split, LeaveOneOut, GridSearchCV, cross_val_score, KFold


def getMatchingColNames(operators, features):
    col_names = []
    for col in features.columns:
        if any(operator in col for operator in operators):
            col_names.append(col)
    return col_names

def summarisedNumericalFeatures(col_names, features):
    numerical_features = features.groupby(["pid"])[col_names].var()
    numerical_features.columns = numerical_features.columns.str.replace("daily", "overallvar")
    return numerical_features

def summarisedCategoricalFeatures(col_names, features):
    categorical_features = features.groupby(["pid"])[col_names].agg(lambda x: int(pd.Series.mode(x)[0]))
    categorical_features.columns = categorical_features.columns.str.replace("daily", "overallmode")
    return categorical_features

def summariseFeatures(features, numerical_operators, categorical_operators, cols_var_threshold):
    numerical_col_names = getMatchingColNames(numerical_operators, features)
    categorical_col_names = getMatchingColNames(categorical_operators, features)
    numerical_features = summarisedNumericalFeatures(numerical_col_names, features)
    categorical_features = summarisedCategoricalFeatures(categorical_col_names, features)
    features = pd.concat([numerical_features, categorical_features], axis=1)
    if cols_var_threshold: # double check the categorical features
        features = dropZeroVarianceCols(features)
    return features

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
# Step 1. Read parameters, features and targets
# Step 2. Extract summarised features based on daily features
# Step 3. Create pipeline
# Step 4. Nested cross validation
# Step 5. Model evaluation
# Step 6. Save results, parameters, and metrics to CSV files
##############################################################



# Step 1. Read parameters, features and targets
# Read parameters
model = snakemake.params["model"]
source = snakemake.params["source"]
summarised = snakemake.params["summarised"]
day_segment = snakemake.params["day_segment"]
scaler = snakemake.params["scaler"]
cv_method = snakemake.params["cv_method"]
cols_var_threshold = snakemake.params["cols_var_threshold"]
numerical_operators = snakemake.params["numerical_operators"]
categorical_operators = snakemake.params["categorical_operators"]
categorical_colnames_demographic_features = snakemake.params["categorical_demographic_features"]
model_hyperparams = snakemake.params["model_hyperparams"][model]
rowsnan_colsnan_days_colsvar_threshold = snakemake.params["rowsnan_colsnan_days_colsvar_threshold"] # thresholds for data cleaning
# Read features and targets
demographic_features = pd.read_csv(snakemake.input["demographic_features"], index_col=["pid"])
targets = pd.read_csv(snakemake.input["targets"], index_col=["pid"])
features = pd.read_csv(snakemake.input["cleaned_features"], parse_dates=["local_date"])
# Compute the proportion of missing value cells among all features
nan_ratio = features.isnull().sum().sum() / (features.shape[0] * features.shape[1])


# Step 2. Extract summarised features based on daily features:
# for categorical features: calculate variance across all days
# for numerical features: calculate mode across all days
if summarised == "summarised":
    features = summariseFeatures(features, numerical_operators, categorical_operators, cols_var_threshold)

categorical_feature_colnames = categorical_colnames_demographic_features + getMatchingColNames(categorical_operators, features)

data = pd.concat([features, demographic_features, targets], axis=1, join="inner")
data_x, data_y = data.drop("target", axis=1), data[["target"]]


# Step 3. Create pipeline
pipeline = createPipeline(model)
cv_class = globals()[cv_method]
inner_cv = cv_class()
outer_cv = cv_class()

# Step 4. Nested cross validation
fold_id, pid, best_params, true_y, pred_y, pred_y_prob = [], [], [], [], [], []
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
    train_x, test_x = train_x.align(test_x, join='outer', axis=1, fill_value=0) # in case we get rid off categorical columns

    # Compute number of participants and features
    # values do not change between folds
    if fold_count == 1:
        num_of_participants = train_x.shape[0] + test_x.shape[0]
        num_of_features = train_x.shape[1]

    # Inner cross validation
    clf = GridSearchCV(estimator=pipeline, param_grid=model_hyperparams, cv=inner_cv, scoring="f1_micro")
    clf.fit(train_x, train_y.values.ravel())

    # Collect results and parameters
    best_params = best_params + [clf.best_params_]
    pred_y = pred_y + clf.predict(test_x).tolist()
    pred_y_prob = pred_y_prob + clf.predict_proba(test_x)[:, 1].tolist()
    true_y = true_y + test_y.values.ravel().tolist()
    pid = pid + test_y.index.tolist() # each test partition (fold) in the outer cv is a participant (LeaveOneOut cv)
    feature_importances_current_fold = getFeatureImportances(model, clf.best_estimator_.steps[1][1], train_x.columns)
    feature_importances_all_folds = pd.concat([feature_importances_all_folds, feature_importances_current_fold], sort=False, axis=0)
    fold_id.append(fold_count)
    fold_count = fold_count + 1

# Step 5. Model evaluation
acc, pre1, recall1, f11, auc, kappa = getMetrics(pred_y, pred_y_prob, true_y)

# Step 6. Save results, parameters, and metrics to CSV files
fold_predictions = pd.DataFrame({"fold_id": fold_id, "pid": pid, "hyperparameters": best_params, "true_y": true_y, "pred_y": pred_y, "pred_y_prob": pred_y_prob})
fold_metrics = pd.DataFrame({"fold_id":[], "accuracy":[], "precision1": [], "recall1": [], "f11": [], "auc": [], "kappa": []})
overall_results = pd.DataFrame({"num_of_participants": [num_of_participants], "num_of_features": [num_of_features], "nan_ratio": [nan_ratio], "rowsnan_colsnan_days_colsvar_threshold": [rowsnan_colsnan_days_colsvar_threshold], "model": [model], "cv_method": [cv_method], "source": [source], "scaler": [scaler], "day_segment": [day_segment], "summarised": [summarised], "accuracy": [acc], "precision1": [pre1], "recall1": [recall1], "f11": [f11], "auc": [auc], "kappa": [kappa]})
feature_importances_all_folds.insert(loc=0, column='fold_id', value=fold_id)
feature_importances_all_folds.insert(loc=1, column='pid', value=pid)

fold_predictions.to_csv(snakemake.output["fold_predictions"], index=False)
fold_metrics.to_csv(snakemake.output["fold_metrics"], index=False)
overall_results.to_csv(snakemake.output["overall_results"], index=False)
feature_importances_all_folds.to_csv(snakemake.output["fold_feature_importances"], index=False)
