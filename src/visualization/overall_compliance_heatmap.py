import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.graph_objects as go
from dateutil import tz
import datetime

def getOneRow(data_per_participant, last_certain_dates, col_name, row, expected_num_of_days, only_show_valid_days):

    data = pd.read_csv(data_per_participant, index_col=["local_date"])

    if col_name == "num_sensors":
        data["num_sensors"] = data.max(axis=1)
    
    if only_show_valid_days and col_name == "valid_sensed_hours":
        # replace invalid days' valid sensed hours with np.nan to let our heatmap only shows valid days
        data.loc[data[data["is_valid_sensed_day"] == False].index, "valid_sensed_hours"] = np.nan

    if expected_num_of_days == -1:
        # show all days
        data.index = pd.to_datetime(data.index)
        start_date = data.index.min()
        # upsample data into one day bins
        data = data.resample("1D").sum()
        data["date_idx"] = (data.index - start_date).days
        data.set_index("date_idx", inplace=True, drop=True)
        row = row + data[col_name].tolist()
    else:
        # only show last certain days
        for date in last_certain_dates:
            if date in data.index:
                row.append(data.loc[date][col_name])
            else:
                row.append(0)

    return row

def getOverallComplianceHeatmap(sensors_with_data, valid_sensed_hours, last_certain_dates, bin_size, min_bins_per_hour, expected_num_of_days, output_path):
    plot = go.Figure(data=go.Heatmap(z=valid_sensed_hours[last_certain_dates].values,
                                       x=[date.replace("-", "/") for date in last_certain_dates] if expected_num_of_days != -1 else last_certain_dates,
                                       y=[pid + "." + label for pid, label in zip(sensors_with_data["pid"].to_list(), sensors_with_data["label"].to_list())],
                                       text=sensors_with_data[last_certain_dates].values,
                                       hovertemplate="Date: %{x}<br>Participant: %{y}<br>Valid sensed hours: %{z}<br>Number of sensors with data: %{text}<extra></extra>" if expected_num_of_days != -1 else "Date_idx: %{x}<br>Participant: %{y}<br>Valid sensed hours: %{z}<br>Number of sensors with data: %{text}<extra></extra>",
                                       colorscale="Viridis",
                                       colorbar={"tick0": 0,"dtick": 1},
                                       showscale=True))
    if expected_num_of_days != -1:
        plot.update_layout(title="Overall compliance heatmap for last " + str(expected_num_of_days) + " days.<br>Bin's color shows valid sensed hours for that day.<br>A valid hour has at least one row of any sensor in "+ str(min_bins_per_hour) +" out of " + str(int(60 / bin_size)) + " bins of " + str(bin_size) + " minutes")
    else:
        plot.update_layout(title="Overall compliance heatmap for all days.<br>Bin's color shows valid sensed hours for that day.<br>A valid hour has at least one row of any sensor in "+ str(min_bins_per_hour) +" out of " + str(int(60 / bin_size)) + " bins of " + str(bin_size) + " minutes")
    
    plot["layout"]["xaxis"].update(side="bottom")
    pio.write_html(plot, file=output_path, auto_open=False, include_plotlyjs="cdn")


phone_sensed_bins = snakemake.input["phone_sensed_bins"]
phone_valid_sensed_days = snakemake.input["phone_valid_sensed_days"]
pid_files = snakemake.input["pid_files"]
only_show_valid_days = snakemake.params["only_show_valid_days"]
local_timezone = snakemake.params["local_timezone"]
bin_size = snakemake.params["bin_size"]
min_bins_per_hour = snakemake.params["min_bins_per_hour"]
expected_num_of_days = int(snakemake.params["expected_num_of_days"])

if expected_num_of_days < -1:
    raise ValueError("EXPECTED_NUM_OF_DAYS of OVERALL_COMPLIANCE_HEATMAP section in config.yaml must be larger or equal to -1.")

last_certain_dates = []
if expected_num_of_days != -1:
    # get the list of dates to show
    cur_date = datetime.datetime.now().astimezone(tz.gettz(local_timezone)).date()
    for date_offset in range(expected_num_of_days-1, -1, -1):
        last_certain_dates.append((cur_date - datetime.timedelta(days=date_offset)).strftime("%Y-%m-%d"))

sensors_with_data_records, valid_sensed_hours_records = [], []
for sensors_with_data_individual, valid_sensed_hours_individual, pid_file in zip(phone_sensed_bins, phone_valid_sensed_days, pid_files):
    
    with open(pid_file, encoding="ISO-8859-1") as external_file:
        external_file_content = external_file.readlines()
    device_id = external_file_content[0].split(",")[-1].strip()
    label = external_file_content[2].strip()
    pid = pid_file.split("/")[-1]

    sensors_with_data_records.append(getOneRow(sensors_with_data_individual, last_certain_dates, "num_sensors", [pid, label, device_id], expected_num_of_days, only_show_valid_days))
    valid_sensed_hours_records.append(getOneRow(valid_sensed_hours_individual, last_certain_dates, "valid_sensed_hours", [pid, label, device_id], expected_num_of_days, only_show_valid_days))

if expected_num_of_days == -1:
    # get the date_idx of all days
    total_num_of_days = max([len(x) for x in sensors_with_data_records]) - 3
    last_certain_dates = [date_idx for date_idx in range(total_num_of_days)]

sensors_with_data = pd.DataFrame(data=sensors_with_data_records, columns=["pid", "label", "device_id"] + last_certain_dates).replace(0, np.nan)
valid_sensed_hours = pd.DataFrame(data=valid_sensed_hours_records, columns=["pid", "label", "device_id"] + last_certain_dates).replace(0, np.nan)

if sensors_with_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no sensor data for all participants")
    empty_html.close()
else:
    getOverallComplianceHeatmap(sensors_with_data, valid_sensed_hours, last_certain_dates, bin_size, min_bins_per_hour, expected_num_of_days, snakemake.output[0])
