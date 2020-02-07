import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.graph_objects as go
import datetime

def getComplianceMatrix(dates, compliance_bins):
    compliance_matrix = []
    for date in dates:
        date_bins = compliance_bins[compliance_bins["local_date"] == date]["count"].tolist()
        compliance_matrix.append(date_bins)
    return compliance_matrix


def getRowCountHeatmap(dates, row_count_per_bin, sensor_name, pid, output_path, bin_size):
    bins_per_hour = int(60 / bin_size)
    x_axis_labels = ["{0:0=2d}".format(x // bins_per_hour) + ":" + \
                    "{0:0=2d}".format(x % bins_per_hour * bin_size) for x in range(24 * bins_per_hour)]
    plot = go.Figure(data=go.Heatmap(z=row_count_per_bin,
                                     x=x_axis_labels,
                                     y=[datetime.datetime.strftime(date, '%Y/%m/%d') for date in dates],
                                     colorscale="Viridis"))
    plot.update_layout(title="Row count heatmap for " + sensor_name + " of " + pid)
    pio.write_html(plot, file=output_path, auto_open=False)



sensor_data = pd.read_csv(snakemake.input[0], encoding="ISO-8859-1")
sensor_name = snakemake.params["table"]
pid = snakemake.params["pid"]
bin_size = snakemake.params["bin_size"]

# check if we have sensor data
if sensor_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no "+ sensor_name + " data for "+pid)
    empty_html.close()
else:
    start_date = sensor_data["local_date"][0]
    end_date = sensor_data.at[sensor_data.index[-1],"local_date"]

    sensor_data["local_date_time"] = pd.to_datetime(sensor_data["local_date_time"])    
    sensor_data = sensor_data[["local_date_time"]]
    sensor_data["count"] = 1

    # Add first and last day boundaries for resampling
    sensor_data = sensor_data.append([pd.Series([datetime.datetime.strptime(start_date + " 00:00:00", "%Y-%m-%d %H:%M:%S"), 0], sensor_data.columns),
                                                            pd.Series([datetime.datetime.strptime(end_date + " 23:59:59", "%Y-%m-%d %H:%M:%S"),  0], sensor_data.columns)])

    # Resample into bins with the size of bin_size
    resampled_bins = pd.DataFrame(sensor_data.resample(str(bin_size) + "T", on="local_date_time")["count"].sum())
    
    # Extract list of dates for creating the heatmap
    resampled_bins.reset_index(inplace=True)
    resampled_bins["local_date"] = resampled_bins["local_date_time"].dt.date
    dates = resampled_bins["local_date"].drop_duplicates().tolist()

    # Create heatmap
    row_count_per_bin = getComplianceMatrix(dates, resampled_bins)
    row_count_per_bin = np.asarray(row_count_per_bin)
    row_count_per_bin = np.where(row_count_per_bin == 0, np.nan, row_count_per_bin)
    getRowCountHeatmap(dates, row_count_per_bin, sensor_name, pid, snakemake.output[0], bin_size)
