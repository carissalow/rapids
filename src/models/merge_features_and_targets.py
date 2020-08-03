import pandas as pd
import numpy as np
from modeling_utils import getMatchingColNames, dropZeroVarianceCols


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
    if cols_var_threshold == "True": # double check the categorical features
        features = dropZeroVarianceCols(features)
    elif cols_var_threshold == "Flase":
        pass
    else:
        ValueError("COLS_VAR_THRESHOLD parameter in config.yaml can only be 'True' or 'False'")
    return features


summarised = snakemake.params["summarised"]
cols_var_threshold = snakemake.params["cols_var_threshold"]
numerical_operators = snakemake.params["numerical_operators"]
categorical_operators = snakemake.params["categorical_operators"]
features_exclude_day_idx = snakemake.params["features_exclude_day_idx"]


# Extract summarised features based on daily features:
# for categorical features: calculate variance across all days
# for numerical features: calculate mode across all days
if summarised == "summarised":

    features = pd.read_csv(snakemake.input["cleaned_features"], parse_dates=["local_date"])
    demographic_features = pd.read_csv(snakemake.input["demographic_features"], index_col=["pid"])
    targets = pd.read_csv(snakemake.input["targets"], index_col=["pid"])

    features = summariseFeatures(features, numerical_operators, categorical_operators, cols_var_threshold)
    data = pd.concat([features, demographic_features, targets], axis=1, join="inner")

elif summarised == "notsummarised":

    features = pd.read_csv(snakemake.input["cleaned_features"])
    demographic_features = pd.read_csv(snakemake.input["demographic_features"])
    
    features = features.merge(demographic_features, on="pid", how="left").set_index(["pid", "local_date"])
    targets = pd.read_csv(snakemake.input["targets"], index_col=["pid", "local_date"])
    data = pd.concat([features, targets], axis=1, join="inner")

else:
    raise ValueError("SUMMARISED parameter in config.yaml can only be 'summarised' or 'notsummarised'")

if features_exclude_day_idx and ("day_idx" in data.columns):
    del data["day_idx"]

data.to_csv(snakemake.output[0], index=True)

