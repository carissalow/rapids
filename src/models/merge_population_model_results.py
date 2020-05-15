import pandas as pd

overall_results = pd.read_csv(snakemake.input["overall_results"])
nan_cells_ratio = pd.read_csv(snakemake.input["nan_cells_ratio"])
baseline = pd.read_csv(snakemake.input["baseline"], index_col=["method"])

# add nan cells ratio
overall_results.insert(3, "nan_cells_ratio", nan_cells_ratio["nan_cells_ratio"])

# add baseline
baseline = baseline.stack().to_frame().T
baseline.columns = ['{}_{}'.format(*col) for col in baseline.columns]
baseline = baseline.add_prefix('b_')
results = pd.concat([overall_results, baseline], axis=1)

results.to_csv(snakemake.output[0], index=False)
