import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go


def getCorrMatrixHeatmap(corr_matrix, output_path):
    colnames = corr_matrix.columns
    plot = go.Figure(data=go.Heatmap(z=corr_matrix.values.tolist(),
                                     x=colnames,
                                     y=colnames,
                                     colorscale="Viridis"))
    plot.update_layout(title="Correlation Matrix Heatmap")
    pio.write_html(plot, file=output_path, auto_open=False, include_plotlyjs="cdn")


min_rows_ratio = snakemake.params["min_rows_ratio"]
corr_threshold = snakemake.params["corr_threshold"]

# merge features
features, features_all_sensors = pd.DataFrame(columns=["local_date"]), pd.DataFrame(columns=["local_date"])
pids = set()
last_pid = None
for path in snakemake.input["features"]:
    pid = path.split("/")[2]
    if pid not in pids:
        pids.add(pid)
        features_all_sensors["pid"] = last_pid
        features = pd.concat([features, features_all_sensors], axis=0, ignore_index=True, sort=False)
        features_all_sensors = pd.DataFrame(columns=["local_date"])
    features_per_sensor = pd.read_csv(path)
    features_all_sensors = features_all_sensors.merge(features_per_sensor, on="local_date", how="outer")
    last_pid = pid

features_all_sensors["pid"] = last_pid
features = pd.concat([features, features_all_sensors], axis=0, ignore_index=True, sort=False)
features.set_index(["pid", "local_date"], inplace=True)

# select days based on the input of "phone_valid_sensed_days"
selected_participants_and_days = pd.DataFrame()
for path in snakemake.input["phone_valid_sensed_days"]:
    pid = path.split("/")[2]
    phone_valid_sensed_days = pd.read_csv(path)
    phone_valid_sensed_days = phone_valid_sensed_days[phone_valid_sensed_days["is_valid_sensed_day"] == True]
    phone_valid_sensed_days["pid"] = pid
    selected_participants_and_days = pd.concat([selected_participants_and_days, phone_valid_sensed_days], axis=0)

selected_participants_and_days.set_index(["pid", "local_date"], inplace=True)
features = features.loc[features.index.intersection(selected_participants_and_days.index), :]

# get correlation matrix
features = features.astype(float)
corr_matrix = features.corr(min_periods=min_rows_ratio * features.shape[0])

# replace correlation coefficients less than corr_threshold with NA
corr_matrix[(corr_matrix > -corr_threshold) & (corr_matrix < corr_threshold)] = np.nan

# plot heatmap
getCorrMatrixHeatmap(corr_matrix, snakemake.output[0])
