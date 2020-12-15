import pandas as pd
from pandas.core import indexing
import yaml


def extract_empatica_data(sensor_data_file, output_file, start_date, end_date, timezone, sensor):
    print(sensor_data_file)
    print(output_file)
    print(start_date)
    print(end_date)
    print(timezone)
    print(sensor)
    data = pd.read_csv(sensor_data_file)
    print(data)

    # extract
    print(output_file)
    data.to_csv(output_file, index = False)


sensor_data_file = snakemake.input[0]
output_file = snakemake.output[0]
with open(snakemake.input[1], "r", encoding="utf-8") as f:
    participant_file = yaml.safe_load(f)

start_date = participant_file["EMPATICA"]["START_DATE"]
end_date = participant_file["EMPATICA"]["END_DATE"]
timezone = snakemake.params["data_configuration"]["TIMEZONE"]["VALUE"]
sensor = snakemake.params["sensor"]

extract_empatica_data(sensor_data_file, output_file, start_date, end_date, timezone, sensor)

