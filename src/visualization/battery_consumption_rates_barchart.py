import pandas as pd
import datetime
import plotly.io as pio
import plotly.graph_objects as go

def getBatteryConsumptionRatesBarChart(battery_data, pid):
    plot = go.Figure(go.Bar(
                    x=battery_data["battery_daily_avgconsumptionrate"],
                    y=battery_data["local_date"].apply(lambda x: x.strftime("%Y/%m/%d")).tolist(),
                    orientation='h'))
    plot.update_layout(title="Daily battery consumption rates bar chart for " + pid + "<br>Label: " + label + ", device_id: " + device_id,
                    xaxis_title="battery drains % per hour",
                    )
    return plot



battery_data = pd.read_csv(snakemake.input["sensor"], parse_dates=["local_date"])
pid = snakemake.params["pid"]

with open(snakemake.input["pid_file"], encoding="ISO-8859-1") as external_file:
    external_file_content = external_file.readlines()
device_id = external_file_content[0].split(",")[-1]
label = external_file_content[2]

if battery_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no battery data for " + pid + "<br>Label: " + label + ", device_id: " + device_id)
    empty_html.close()
else:
    battery_data.set_index(["local_date"], inplace=True)
    battery_data = battery_data.resample("1D").asfreq().fillna(0).reset_index()
    plot = getBatteryConsumptionRatesBarChart(battery_data, pid)
    pio.write_html(plot, file=snakemake.output[0], auto_open=False)