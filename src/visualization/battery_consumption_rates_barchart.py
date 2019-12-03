import pandas as pd
import datetime
import plotly.io as pio
import plotly.graph_objects as go

def getBatteryConsumptionRatesBarChart(battery_data, pid):
    plot = go.Figure(go.Bar(
                    x=battery_data["battery_daily_avgconsumptionrate"],
                    y=battery_data["local_date"].apply(lambda x: x.strftime("%Y/%m/%d")).tolist(),
                    orientation='h'))
    plot.update_layout(title="Daily battery consumption rates bar chart for " + pid,
                    xaxis_title="battery drains % per hour",
                    )
    return plot



battery_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date"])
pid = snakemake.params["pid"]
if battery_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no battery data for " + pid)
    empty_html.close()
else:
    battery_data.set_index(["local_date"], inplace=True)
    battery_data = battery_data.resample("1D").asfreq().fillna(0).reset_index()
    plot = getBatteryConsumptionRatesBarChart(battery_data, pid)
    pio.write_html(plot, file=snakemake.output[0], auto_open=False)