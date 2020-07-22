import pandas as pd
import plotly.express as px
import plotly.io as pio


# merge "phone_valid_sensed_days" for all participants
selected_participants_and_days = pd.DataFrame()
for path in snakemake.input["phone_valid_sensed_days"]:
    phone_valid_sensed_days = pd.read_csv(path)
    phone_valid_sensed_days = phone_valid_sensed_days[phone_valid_sensed_days["is_valid_sensed_day"] == True]
    selected_participants_and_days = pd.concat([selected_participants_and_days, phone_valid_sensed_days], axis=0)

# plot histogram
fig = px.histogram(selected_participants_and_days, x="valid_sensed_hours")
fig.update_layout(title="Phone Valid Hours Histogram")
pio.write_html(fig, file=snakemake.output[0], auto_open=False, include_plotlyjs="cdn")