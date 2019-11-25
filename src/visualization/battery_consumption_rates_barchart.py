import pandas as pd
import datetime
import plotly.io as pio
import plotly.graph_objects as go

def getBatteryConsumptionRatesBarChart(battery_data, pid):
    plot = go.Figure(go.Bar(
                    x=battery_data["battery_daily_avgconsumptionrate"],
                    y=battery_data["local_date"].apply(lambda x: x.replace("-","/")).tolist(),
                    orientation='h'))
    plot.update_layout(title="Daily battery consumption rates bar chart for " + pid,
                    xaxis_title="battery drains % per hour",
                    )
    return plot
    


battery_data = pd.read_csv(snakemake.input[0])
pid = snakemake.params["pid"]
if battery_data.empty:
    empty_html = open(snakemake.output[0], "w")
    empty_html.write("There is no battery data for " + pid)
    empty_html.close()
else:
    plot = getBatteryConsumptionRatesBarChart(battery_data, pid)
    pio.write_html(plot, file=snakemake.output[0], auto_open=False)