import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go
import datetime

def getComplianceMatrix(dates, compliance_bins):
    compliance_matrix = []
    for date in dates:
        date_bins = compliance_bins[compliance_bins["local_date"] == date]["count"].tolist()
        compliance_matrix.append(date_bins)
    return compliance_matrix


def getHourlyRowCountHeatmap(dates, hourly_row_count, sensor_name, pid, output_path):
    plot = go.Figure(data=go.Heatmap(z=hourly_row_count,
                                     x=[x for x in range(24)],
                                     y=[datetime.datetime.strftime(date, '%Y/%m/%d') for date in dates],
                                     colorscale='Viridis'))
    plot.update_layout(title="Row count heatmap for " + sensor_name + " of " + pid)
    pio.write_html(plot, file=output_path, auto_open=False)



sensor_data = pd.read_csv(snakemake.input[0])
sensor_name = snakemake.params["table"]
pid = snakemake.params["pid"]

# check if we have sensor data
if sensor_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no "+ sensor_name + " data for "+pid)
    empty_html.close()
else:
    start_date = sensor_data["local_date"][0]
    end_date = sensor_data.at[sensor_data.index[-1],"local_date"]

    # Make local hour double digit
    sensor_data["local_hour"] = sensor_data["local_hour"].map("{0:0=2d}".format)

    # Group and count by local_date and local_hour
    sensor_data_hourly_bins = sensor_data.groupby(["local_date","local_hour"]).agg(count=("timestamp","count")).reset_index()

    # Add first and last day boundaries for resampling
    sensor_data_hourly_bins = sensor_data_hourly_bins.append([pd.Series([start_date, "00", 0], sensor_data_hourly_bins.columns),
                                                            pd.Series([end_date, "23", 0], sensor_data_hourly_bins.columns)])

    # Rebuild local date hour for resampling
    sensor_data_hourly_bins["local_date_hour"] = pd.to_datetime(sensor_data_hourly_bins["local_date"] + \
                                                " " + sensor_data_hourly_bins["local_hour"] + ":00:00")

    resampled_hourly_bins = pd.DataFrame(sensor_data_hourly_bins.resample("1H", on="local_date_hour")["count"].sum())

    # Extract list of dates for creating the heatmap
    resampled_hourly_bins.reset_index(inplace=True)
    resampled_hourly_bins["local_date"] = resampled_hourly_bins["local_date_hour"].dt.date
    dates = resampled_hourly_bins["local_date"].drop_duplicates().tolist()

    # Create heatmap
    hourly_row_count = getComplianceMatrix(dates, resampled_hourly_bins)
    getHourlyRowCountHeatmap(dates, hourly_row_count, sensor_name, pid, snakemake.output[0])
