import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.graph_objects as go
import datetime

def getDatesComplianceMatrix(phone_sensed_bins):
    dates = phone_sensed_bins.index
    compliance_matrix = []
    for date in dates:
        compliance_matrix.append(phone_sensed_bins.loc[date, :].tolist())
    return dates, compliance_matrix

def getComplianceHeatmap(dates, compliance_matrix, pid, output_path, bin_size):
    bins_per_hour = int(60 / bin_size)
    x_axis_labels = ["{0:0=2d}".format(x // bins_per_hour) + ":" + \
                    "{0:0=2d}".format(x % bins_per_hour * bin_size) for x in range(24 * bins_per_hour)]
    plot = go.Figure(data=go.Heatmap(z=compliance_matrix,
                                     x=x_axis_labels,
                                     y=[datetime.datetime.strftime(date, '%Y/%m/%d') for date in dates],
                                     colorscale='Viridis',
                                     colorbar={'tick0': 0,'dtick': 1}))
    plot.update_layout(title="Compliance heatmap.<br>Five-minute bins showing how many sensors logged at least one row of data in that period for " + pid + "<br>Label: " + label + ", device_id: " + device_id)
    pio.write_html(plot, file=output_path, auto_open=False, include_plotlyjs="cdn")

# get current patient id
pid = snakemake.params["pid"]
bin_size = snakemake.params["bin_size"]

with open(snakemake.input["pid_file"], encoding="ISO-8859-1") as external_file:
    external_file_content = external_file.readlines()
device_id = external_file_content[0].split(",")[-1]
label = external_file_content[2]

phone_sensed_bins = pd.read_csv(snakemake.input["sensor"], parse_dates=["local_date"], index_col="local_date")

if phone_sensed_bins.empty:
    empty_html = open(snakemake.output[0], "w", encoding="ISO-8859-1")
    empty_html.write("There is no sensor data for " + pid + "<br>Label: " + label + ", device_id: " + device_id)
    empty_html.close()
else:
    # resample to impute missing dates
    phone_sensed_bins = phone_sensed_bins.resample("1D").asfreq().fillna(0)
    # get dates and compliance_matrix
    dates, compliance_matrix = getDatesComplianceMatrix(phone_sensed_bins)
    # convert compliance_matrix from list to np.array and replace 0 with np.nan
    compliance_matrix = np.asarray(compliance_matrix)
    compliance_matrix = np.where(compliance_matrix == 0, np.nan, compliance_matrix)
    # get heatmap
    getComplianceHeatmap(dates, compliance_matrix, pid, snakemake.output[0], bin_size)