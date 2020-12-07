import pandas as pd
import numpy as np
import plotly.graph_objects as go
from importlib import util
from pathlib import Path
import yaml


def getRowCountHeatmap(data_for_plot, scaled_data_for_plot, pid, time_segment, html_file):

    fig = go.Figure(data=go.Heatmap(z=scaled_data_for_plot.values.tolist(),
                                     x=data_for_plot.columns,
                                     y=data_for_plot.index,
                                     hovertext=data_for_plot.values.tolist(),
                                     hovertemplate="Segment start: %{x}<br>Sensor: %{y}<br>Row count: %{hovertext}<extra></extra>",
                                     zmin=0, zmax=1,
                                     colorscale='Viridis'))

    fig.update_layout(title="Heatmap of sensor row count for " + time_segment + " segments. Pid: " + pid +". Label: " + label + "<br>y-axis shows the included sensors.<br>x-axis shows the start (date and time) of a time segment.<br>z-axis (color) shows row count per sensor per segment instance.")
    fig["layout"].update(margin=dict(t=160))
    
    html_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))




# import filter_data_by_segment from src/features/utils/utils.py
spec = util.spec_from_file_location("util", str(Path(snakemake.scriptdir).parent / "features" / "utils" / "utils.py"))
mod = util.module_from_spec(spec)
spec.loader.exec_module(mod)
filter_data_by_segment = getattr(mod,  "filter_data_by_segment")





phone_data_yield = pd.read_csv(snakemake.input["phone_data_yield"], index_col=["local_segment_start_datetime"], parse_dates=["local_segment_start_datetime"])
# make sure the phone_data_yield file contains "phone_data_yield_rapids_ratiovalidyieldedminutes" and "phone_data_yield_rapids_ratiovalidyieldedhours" columns
if ("phone_data_yield_rapids_ratiovalidyieldedminutes" not in phone_data_yield.columns) or ("phone_data_yield_rapids_ratiovalidyieldedhours" not in phone_data_yield.columns):
    raise ValueError("Please make sure [PHONE_DATA_YIELD][RAPIDS][COMPUTE] is True AND [PHONE_DATA_YIELD][RAPIDS][FEATURES] contains [ratiovalidyieldedminutes, ratiovalidyieldedhours].")
phone_data_yield = phone_data_yield[["local_segment_label", "phone_data_yield_rapids_ratiovalidyieldedminutes", "phone_data_yield_rapids_ratiovalidyieldedhours"]]

time_segments = pd.read_csv(snakemake.input["time_segments_labels"], header=0)["label"]
pid = snakemake.params["pid"]

with open(snakemake.input["participant_file"], "r", encoding="utf-8") as f:
    participant_file = yaml.safe_load(f)
label = participant_file["PHONE"]["LABEL"]

sensor_names = []
sensors_row_count = dict(zip(time_segments, [pd.DataFrame()] * len(time_segments)))

for sensor_path in snakemake.input["all_sensors"]:
    sensor_data = pd.read_csv(sensor_path, usecols=["assigned_segments"])
    sensor_name = sensor_path.split("/")[-1].replace("_with_datetime.csv", "")
    sensor_names.append(sensor_name)
    
    if not sensor_data.empty:
        for time_segment in time_segments:
            sensor_data_per_segment = filter_data_by_segment(sensor_data, time_segment)

            if not sensor_data_per_segment.empty:
                # extract local start datetime of the segment from "local_segment" column
                sensor_data_per_segment["local_segment_start_datetime"] = pd.to_datetime(sensor_data_per_segment["local_segment"].apply(lambda x: x.split("#")[1].split(",")[0]))
                sensor_row_count = sensor_data_per_segment.groupby("local_segment_start_datetime")[["local_segment"]].count().rename(columns={"local_segment": sensor_name})
                sensors_row_count[time_segment] = pd.concat([sensors_row_count[time_segment], sensor_row_count], axis=1, sort=False)

# add phone data yield features and plot heatmap
html_file = open(snakemake.output[0], "a", encoding="utf-8")
sensor_names.extend(["ratiovalidyieldedminutes", "ratiovalidyieldedhours"])
for time_segment in time_segments:
    if not phone_data_yield.empty:
        phone_data_yield_per_segment = phone_data_yield[phone_data_yield["local_segment_label"] == time_segment].rename(columns={"phone_data_yield_rapids_ratiovalidyieldedminutes": "ratiovalidyieldedminutes","phone_data_yield_rapids_ratiovalidyieldedhours": "ratiovalidyieldedhours"}).round(3)
        if not phone_data_yield_per_segment.empty:
            sensors_row_count[time_segment] = pd.concat([sensors_row_count[time_segment], phone_data_yield_per_segment], axis=1, sort=True)
    
    # consider all the sensors
    data_for_plot = sensors_row_count[time_segment].transpose().reindex(pd.Index(sensor_names))

    if data_for_plot.empty:
        html_file.write("There are no records of selected sensors in database for " + time_segment + " segments. Pid: " + pid + ". Label: " + label + ".<br>")
    else:
        # except for phone data yield sensor, scale each sensor (row) to the range of [0, 1]
        scaled_data_for_plot = data_for_plot.copy()
        scaled_data_for_plot.loc[sensor_names[:-2]] = scaled_data_for_plot.fillna(np.nan).loc[sensor_names[:-2]].apply(lambda x: (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)) if np.nanmax(x) != np.nanmin(x) else (x / np.nanmin(x)), axis=1)

        getRowCountHeatmap(data_for_plot, scaled_data_for_plot, pid, time_segment, html_file)

html_file.close()
