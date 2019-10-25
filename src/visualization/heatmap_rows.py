import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.graph_objects as go

def getHourlyRowCount(dates, sensor_data):
    hourly_row_count = []
    for date in dates:
        num_rows = []
        daily_rows = sensor_data[sensor_data["local_date"] == date]
        for hour in range(24):
            hourly_rows = daily_rows[daily_rows["local_hour"] == hour]
            num_rows.append(hourly_rows.shape[0])
        hourly_row_count.append(num_rows)
    return hourly_row_count

def getHourlyRowCountHeatmap(dates, hourly_row_count, sensor_name, pid, output_path):
    plot = go.Figure(data=go.Heatmap(z=hourly_row_count,x=[x for x in range(24)],y=dates,colorscale='Viridis'))
    plot.update_layout(title="Hourly row count heatmap for " + pid + " for sensor " + sensor_name)
    pio.write_html(plot, file=output_path, auto_open=False)



sensor_data = pd.read_csv(snakemake.input[0])
# get current sensor name
sensor_name = snakemake.params["table"]
# get current patient id
pid = snakemake.params["pid"]
# get sorted date list
dates = list(set(sensor_data["local_date"]))
dates.sort()
# get num of rows per hour per day
hourly_row_count = getHourlyRowCount(dates, sensor_data)
# get heatmap
getHourlyRowCountHeatmap(dates, hourly_row_count, sensor_name, pid, snakemake.output[0])
