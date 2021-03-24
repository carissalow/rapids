import pandas as pd
import numpy as np
import plotly.graph_objects as go
import yaml





def getPhoneDataYieldHeatmap(data_for_plot, y_axis_labels, time_segment, type, time, html_file):

    fig = go.Figure(data=go.Heatmap(z=data_for_plot.values.tolist(),
                                       x=data_for_plot.columns.tolist(),
                                       y=y_axis_labels,
                                       hovertext=data_for_plot.values.tolist(),
                                       hovertemplate="Time since first segment: %{x}<br>Participant: %{y}<br>Ratiovalidyielded" + type + ": %{z}<extra></extra>" if time == "RELATIVE_TIME" else "Time: %{x}<br>Participant: %{y}<br>Ratiovalidyielded" + type + ": %{z}<extra></extra>",
                                       zmin=0, zmax=1,
                                       colorscale="Viridis"))

    if time == "RELATIVE_TIME":
        fig.update_layout(title="Heatmap of valid yielded " + type + " ratio for " + time_segment + " segments.<br>y-axis shows participant information (format: pid.label).<br>x-axis shows the time since their first segment.<br>z-axis (color) shows valid yielded " + type + " ratio during a segment instance.")
    else:
        fig.update_layout(title="Heatmap of valid yielded " + type + " ratio for " + time_segment + " segments.<br>y-axis shows participant information (format: pid.label).<br>x-axis shows the time.<br>z-axis (color) shows valid yielded " + type + " ratio during a segment instance.")

    fig["layout"]["xaxis"].update(side="bottom")
    fig["layout"].update(xaxis_title="Time Since First Segment" if time == "RELATIVE_TIME" else "Time")
    fig["layout"].update(margin=dict(t=160))
    
    html_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))






time = snakemake.params["time"]
y_axis_labels, phone_data_yield_minutes, phone_data_yield_hours = [], {}, {}
for phone_data_yield_path, participant_file_path, time_segments_path in zip(snakemake.input["phone_data_yield"], snakemake.input["participant_file"], snakemake.input["time_segments_labels"]):
    
    # set pid.label as y_axis_label
    pid = phone_data_yield_path.split("/")[3]
    time_segments = pd.read_csv(time_segments_path, header=0)["label"]

    with open(participant_file_path, "r", encoding="utf-8") as f:
        participant_file = yaml.safe_load(f)
    label = participant_file["PHONE"]["LABEL"]

    y_axis_label = pid + "." + label
    y_axis_labels.append(y_axis_label)

    
    phone_data_yield = pd.read_csv(phone_data_yield_path, index_col=["local_segment_start_datetime"], parse_dates=["local_segment_start_datetime"])
    # make sure the phone_data_yield file contains "phone_data_yield_rapids_ratiovalidyieldedminutes" and "phone_data_yield_rapids_ratiovalidyieldedhours" columns
    if ("phone_data_yield_rapids_ratiovalidyieldedminutes" not in phone_data_yield.columns) or ("phone_data_yield_rapids_ratiovalidyieldedhours" not in phone_data_yield.columns):
        raise ValueError("Please make sure [PHONE_DATA_YIELD][RAPIDS][COMPUTE] is True AND [PHONE_DATA_YIELD][RAPIDS][FEATURES] contains [ratiovalidyieldedminutes, ratiovalidyieldedhours].")

    if not phone_data_yield.empty:

        for time_segment in time_segments:
            phone_data_yield_per_segment = phone_data_yield[phone_data_yield["local_segment_label"] == time_segment]

            if not phone_data_yield_per_segment.empty:

                if time == "RELATIVE_TIME":
                    # set number of minutes after the first start date time of local segments as x_axis_label
                    phone_data_yield_per_segment.index = phone_data_yield_per_segment.index - phone_data_yield_per_segment.index.min()
                elif time == "ABSOLUTE_TIME":
                    pass
                else:
                    raise ValueError("[HEATMAP_PHONE_DATA_YIELD_PER_PARTICIPANT_PER_TIME_SEGMENT][TIME] can only be RELATIVE_TIME or ABSOLUTE_TIME")

                phone_data_yield_minutes_per_segment = phone_data_yield_per_segment[["phone_data_yield_rapids_ratiovalidyieldedminutes"]].rename(columns={"phone_data_yield_rapids_ratiovalidyieldedminutes": y_axis_label})
                phone_data_yield_hours_per_segment = phone_data_yield_per_segment[["phone_data_yield_rapids_ratiovalidyieldedhours"]].rename(columns={"phone_data_yield_rapids_ratiovalidyieldedhours": y_axis_label})

                if time_segment not in phone_data_yield_minutes.keys():
                    phone_data_yield_minutes[time_segment] = phone_data_yield_minutes_per_segment
                    phone_data_yield_hours[time_segment] = phone_data_yield_hours_per_segment
                else:
                    phone_data_yield_minutes[time_segment] = pd.concat([phone_data_yield_minutes[time_segment], phone_data_yield_minutes_per_segment], axis=1, sort=True)
                    phone_data_yield_hours[time_segment] = pd.concat([phone_data_yield_hours[time_segment], phone_data_yield_hours_per_segment], axis=1, sort=True)


html_file = open(snakemake.output[0], "a", encoding="utf-8")
if len(phone_data_yield_minutes.keys()) == 0:
    html_file.write("There is no sensor data for the sensors in [PHONE_DATA_YIELD][SENSORS].")
for time_segment in phone_data_yield_minutes.keys():
    minutes_data_for_plot = phone_data_yield_minutes[time_segment].transpose().reindex(pd.Index(y_axis_labels)).round(3)
    hours_data_for_plot = phone_data_yield_hours[time_segment].transpose().reindex(pd.Index(y_axis_labels)).round(3)

    getPhoneDataYieldHeatmap(minutes_data_for_plot, y_axis_labels, time_segment, "minutes", time, html_file)
    getPhoneDataYieldHeatmap(hours_data_for_plot, y_axis_labels, time_segment, "hours", time, html_file)

html_file.close()
