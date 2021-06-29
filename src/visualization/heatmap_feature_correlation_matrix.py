import numpy as np
import pandas as pd
import plotly.graph_objects as go


def getCorrMatrixHeatmap(corr_matrix, time_segment, html_file):

    feature_names = corr_matrix.columns

    fig = go.Figure(data=go.Heatmap(z=corr_matrix.values.tolist(),
                                     x=feature_names,
                                     y=feature_names,
                                     colorscale="Viridis",
                                     zmin=-1, zmax=1))
    
    fig.update_layout(title="Correlation matrix of features of " + time_segment + " segments.")
    
    html_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))


time_segments_type = snakemake.params["time_segments_type"]
min_rows_ratio = snakemake.params["min_rows_ratio"]
corr_threshold = snakemake.params["corr_threshold"]
corr_method = snakemake.params["corr_method"]

features = pd.read_csv(snakemake.input["all_sensor_features"])


if time_segments_type == "FREQUENCY":
    features["local_segment_label"] = features["local_segment_label"].str[:-4]
if time_segments_type == "EVENT":
    features["local_segment_label"] = "event"

time_segments = set(features["local_segment_label"])

html_file = open(snakemake.output[0], "a", encoding="utf-8")
if features.empty:
    html_file.write("There are no features for any participant.")
else:

    for time_segment in time_segments:
        features_per_segment = features[features["local_segment_label"] == time_segment]
        if features_per_segment.empty:
            html_file.write("There are no features for " + time_segment + " segments.<br>")
        else:
            # drop useless columns
            features_per_segment = features_per_segment.drop(["pid", "local_segment", "local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"], axis=1).astype(float)
            # get correlation matrix
            corr_matrix = features_per_segment.corr(method=corr_method, min_periods=min_rows_ratio * features_per_segment.shape[0])
            # replace correlation coefficients less than corr_threshold with NA
            corr_matrix[(corr_matrix > -corr_threshold) & (corr_matrix < corr_threshold)] = np.nan

            # plot heatmap
            getCorrMatrixHeatmap(corr_matrix, time_segment, html_file)

html_file.close()
