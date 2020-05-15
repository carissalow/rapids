import pandas as pd

features = pd.read_csv(snakemake.input["cleaned_features"], parse_dates=["local_date"])

# Compute the proportion of missing value cells among all features
nan_cells_ratio = features.isnull().sum().sum() / (features.shape[0] * features.shape[1])

pd.DataFrame({"nan_cells_ratio": [nan_cells_ratio]}).to_csv(snakemake.output[0], index=False)