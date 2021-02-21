import pandas as pd
from pandas.core import indexing
import yaml

import csv
from collections import OrderedDict


def processAcceleration(x, y, z):
    x = float(x)
    y = float(y)
    z = float(z)
    return {'x': x, 'y': y, 'z': z}


def readFile(file, dtype):
    dict = OrderedDict()

    with open(file, 'rt') as csvfile:
        if dtype in ('electrodermal_activity', 'temperature', 'heartrate', 'blood_volume_pulse'):
            reader = csv.reader(csvfile, delimiter='\n')
        elif dtype == 'accelerometer':
            reader = csv.reader(csvfile, delimiter=',')
        i = 0
        for row in reader:
            if i == 0:
                timestamp = float(row[0])
            elif i == 1:
                hertz = float(row[0])
            else:
                if i == 2:
                    pass
                else:
                    timestamp = timestamp + 1.0 / hertz
                if dtype in ('electrodermal_activity', 'temperature', 'heartrate', 'blood_volume_pulse'):
                    dict[timestamp] = row[0]
                elif dtype == 'accelerometer':
                    dict[timestamp] = processAcceleration(row[0], row[1], row[2])
            i += 1
    return dict


def extract_empatica_data(sensor_data_file, output_file, start_date, end_date, timezone, sensor):
    # read sensor data
    if sensor in ('electrodermal_activity', 'temperature', 'heartrate', 'blood_volume_pulse'):
        ddict = readFile(sensor_data_file, sensor)
        df = pd.DataFrame.from_dict(ddict, orient='index', columns=[sensor])
        df[sensor] = df[sensor].astype(float)
        df.index.name = 'timestamp'

    elif sensor == 'accelerometer':
        ddict = readFile(sensor_data_file, sensor)
        df = pd.DataFrame.from_dict(ddict, orient='index', columns=['x', 'y', 'z'])
        df['x'] = df['x'].astype(float)
        df['y'] = df['y'].astype(float)
        df['z'] = df['z'].astype(float)
        df.index.name = 'timestamp'

    elif sensor == 'inter_beat_interval':
        df = pd.read_csv(sensor_data_file, names=['timestamp', sensor], header=None)
        timestampstart = float(df['timestamp'][0])
        df['timestamp'] = (df['timestamp'][1:len(df)]).astype(float) + timestampstart
        df = df.drop([0])
        df[sensor] = df[sensor].astype(float)
        df = df.set_index('timestamp')

    else:
        raise ValueError(
            "sensor can only be one of ['electrodermal_activity','temperature','heartrate','blood_volume_pulse','accelerometer','inter_beat_interval'].")

    # filter based on given start and end date
    start_date_utc = pd.Timestamp(start_date, tz=timezone).timestamp()
    end_date_utc = pd.Timestamp(end_date, tz=timezone).timestamp()
    df = df[start_date_utc:end_date_utc]

    # format timestamps
    df.index *= 1000
    df.index = df.index.astype(int)

    # output csv file
    df.to_csv(output_file)


sensor_data_file = snakemake.input[0]
output_file = snakemake.output[0]
with open(snakemake.input[1], "r", encoding="utf-8") as f:
    participant_file = yaml.safe_load(f)

start_date = participant_file["EMPATICA"]["START_DATE"]
end_date = participant_file["EMPATICA"]["END_DATE"]
timezone = snakemake.params["data_configuration"]["TIMEZONE"]["VALUE"]
sensor = snakemake.params["sensor"]

extract_empatica_data(sensor_data_file, output_file, start_date, end_date, timezone, sensor)
