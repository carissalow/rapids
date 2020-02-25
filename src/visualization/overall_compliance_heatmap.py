import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.figure_factory as ff
from dateutil import tz
import datetime

def getOneRow(data_per_participant, last_seven_dates, col_name, row):
    data = pd.read_csv(data_per_participant, index_col=["local_date"])
    if col_name == "num_sensors":
        data["num_sensors"] = data.max(axis=1)
    for date in last_seven_dates:
        if date in data.index:
            row.append(data.loc[date][col_name])
        else:
            row.append(0)
    return row

def getOverallComplianceHeatmap(sensors_with_data, valid_sensed_hours, last_seven_dates, bin_size, min_bins_per_hour, output_path):
    plot = ff.create_annotated_heatmap(z=sensors_with_data[last_seven_dates].values,
                                       x=[date.replace("-", "/") for date in last_seven_dates],
                                       y=[pid + "<br>" + label for pid, label in zip(sensors_with_data["pid"].to_list(), sensors_with_data["label"].to_list())],
                                       annotation_text=valid_sensed_hours[last_seven_dates].values,
                                       hovertemplate='Date: %{x}<br>Participant: %{y}<br>Number of sensors with data: %{z}<extra></extra>',
                                       colorscale="Viridis",
                                       colorbar={"tick0": 0,"dtick": 1},
                                       showscale=True)
    plot.update_layout(title="Overall compliance heatmap for last seven days.<br>Bin's color shows how many sensors logged at least one row of data for that day.<br>Bin's text shows the valid hours of that day.(A valid hour has at least one row of any sensor in "+ str(min_bins_per_hour) +" out of " + str(int(60 / bin_size)) + " bins of " + str(bin_size) + " minutes)")
    plot["layout"]["xaxis"].update(side="bottom")
    pio.write_html(plot, file=output_path, auto_open=False)


phone_sensed_bins = snakemake.input["phone_sensed_bins"]
phone_valid_sensed_days = snakemake.input["phone_valid_sensed_days"]
pid_files = snakemake.input["pid_files"]
local_timezone = snakemake.params["local_timezone"]
bin_size = snakemake.params["bin_size"]
min_bins_per_hour = snakemake.params["min_bins_per_hour"]


cur_date = datetime.datetime.now().astimezone(tz.gettz(local_timezone)).date()
last_seven_dates = []
for date_offset in range(6,-1,-1):
    last_seven_dates.append((cur_date - datetime.timedelta(days=date_offset)).strftime("%Y-%m-%d"))


sensors_with_data_records, valid_sensed_hours_records = [], []
for sensors_with_data_individual, valid_sensed_hours_individual, pid_file in zip(phone_sensed_bins, phone_valid_sensed_days, pid_files):
    
    with open(pid_file) as external_file:
        external_file_content = external_file.readlines()
    device_id = external_file_content[0].split(",")[-1].strip()
    label = external_file_content[2].strip()
    pid = pid_file.split("/")[-1]

    sensors_with_data_records.append(getOneRow(sensors_with_data_individual, last_seven_dates, "num_sensors", [pid, label, device_id]))
    valid_sensed_hours_records.append(getOneRow(valid_sensed_hours_individual, last_seven_dates, "valid_hours", [pid, label, device_id]))

sensors_with_data = pd.DataFrame(data=sensors_with_data_records, columns=["pid", "label", "device_id"] + last_seven_dates)
valid_sensed_hours = pd.DataFrame(data=valid_sensed_hours_records, columns=["pid", "label", "device_id"] + last_seven_dates)

if sensors_with_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no sensor data for all participants")
    empty_html.close()
else:
    getOverallComplianceHeatmap(sensors_with_data, valid_sensed_hours, last_seven_dates, bin_size, min_bins_per_hour, snakemake.output[0])