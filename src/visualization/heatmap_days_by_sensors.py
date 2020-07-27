import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go
from datetime import datetime, timedelta

def getRowCountHeatmap(row_count_sensors_normalized, row_count_sensors, pid, output_path):
    plot = go.Figure(data=go.Heatmap(z=row_count_sensors_normalized.T.values.tolist(),
                                     x=[datetime.strftime(idx[0], "%Y/%m/%d")+"("+str(idx[1])+")" for idx in row_count_sensors.index],
                                     y=row_count_sensors.columns.tolist(),
                                     hovertext=row_count_sensors.T.values.tolist(),
                                     hovertemplate="Date: %{x}<br>Sensor: %{y}<br>Row count: %{hovertext}<extra></extra>",
                                     colorscale="Viridis"))
    plot.update_layout(title="Row count heatmap for " + pid)
    pio.write_html(plot, file=output_path, auto_open=False, include_plotlyjs="cdn")



phone_valid_sensed_days = pd.read_csv(snakemake.input["phone_valid_sensed_days"], parse_dates=["local_date"], index_col=["local_date"])
phone_valid_sensed_days = phone_valid_sensed_days[phone_valid_sensed_days["is_valid_sensed_day"] == True]

row_count_sensors = pd.DataFrame()
for sensor_path in snakemake.input["sensors"]:
    # plugin_studentlife_audio_android => conversion; plugin_google_activity_recognition => AR; applications_foreground => apps
    sensor_name = sensor_path.split("/")[-1].replace("_with_datetime.csv", "").replace("plugin_studentlife_audio_android", "conversion").replace("plugin_google_activity_recognition", "AR").replace("applications_foreground", "apps")
    sensor_data = pd.read_csv(sensor_path, encoding="ISO-8859-1", parse_dates=["local_date"], dtype={"label": str})
    if sensor_data.empty:
        row_count_sensor = pd.DataFrame(columns=[sensor_name])
    else:
        row_count_sensor = sensor_data[["timestamp", "local_date"]].groupby(["local_date"]).count().rename(columns={"timestamp": sensor_name})
    row_count_sensors = row_count_sensors.join(row_count_sensor, how="outer")

row_count_sensors.index = pd.to_datetime(row_count_sensors.index)
row_count_sensors = row_count_sensors.join(phone_valid_sensed_days[["valid_sensed_hours"]], how="outer")

# set date_idx based on the first date
reference_date = row_count_sensors.index.min()
last_date = row_count_sensors.index.max()
row_count_sensors["date_idx"] = (row_count_sensors.index - reference_date).days
row_count_sensors["local_date"] = row_count_sensors.index
row_count_sensors.set_index(["local_date", "date_idx"], inplace=True)


expected_num_of_days = int(snakemake.params["expected_num_of_days"])
if expected_num_of_days < -1:
    raise ValueError("EXPECTED_NUM_OF_DAYS of HEATMAP_DAYS_BY_SENSORS section in config.yaml must be larger or equal to -1.")
# if expected_num_of_days = -1, return all dates
expected_num_of_days = (last_date - reference_date).days if expected_num_of_days == -1 else expected_num_of_days

# add empty rows to make sure different participants have the same date_idx range
date_idx_range = [idx for idx in range(expected_num_of_days)]
date_range = [reference_date + timedelta(days=idx) for idx in date_idx_range]
all_dates = pd.DataFrame({"local_date": date_range, "date_idx": date_idx_range})
all_dates.set_index(["local_date", "date_idx"], inplace=True)

row_count_sensors = row_count_sensors.merge(all_dates, left_index=True, right_index=True, how="right")

# normalize each sensor (column)
if row_count_sensors.count().max() > 1:
    row_count_sensors_normalized = (row_count_sensors-row_count_sensors.min())/(row_count_sensors.max()-row_count_sensors.min())
else:
    row_count_sensors_normalized = row_count_sensors

pid = sensor_path.split("/")[2]
getRowCountHeatmap(row_count_sensors_normalized, row_count_sensors, pid, snakemake.output[0])
