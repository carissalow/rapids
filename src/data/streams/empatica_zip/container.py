from zipfile import ZipFile
import warnings
from pathlib import Path
import pandas as pd
from pandas.core import indexing
import yaml
import csv
from collections import OrderedDict
from io import BytesIO, StringIO

def processAcceleration(x, y, z):
    x = float(x)
    y = float(y)
    z = float(z)
    return {'x': x, 'y': y, 'z': z}


def readFile(file, dtype):
    dict = OrderedDict()
    # file is an in-memory buffer
    with file as csvfile:
        if dtype in ('EMPATICA_ELECTRODERMAL_ACTIVITY', 'EMPATICA_TEMPERATURE', 'EMPATICA_HEARTRATE', 'EMPATICA_BLOOD_VOLUME_PULSE'):
            reader = csv.reader(csvfile, delimiter='\n')
        elif dtype == 'EMPATICA_ACCELEROMETER':
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
                if dtype in ('EMPATICA_ELECTRODERMAL_ACTIVITY', 'EMPATICA_TEMPERATURE', 'EMPATICA_HEARTRATE', 'EMPATICA_BLOOD_VOLUME_PULSE'):
                    dict[timestamp] = row[0]
                elif dtype == 'EMPATICA_ACCELEROMETER':
                    dict[timestamp] = processAcceleration(row[0], row[1], row[2])
            i += 1
    return dict


def extract_empatica_data(data,  sensor):
    sensor_data_file = BytesIO(data).getvalue().decode('utf-8')
    sensor_data_file = StringIO(sensor_data_file)
    column = sensor.replace("EMPATICA_", "").lower()
    # read sensor data
    if sensor in ('EMPATICA_ELECTRODERMAL_ACTIVITY', 'EMPATICA_TEMPERATURE', 'EMPATICA_HEARTRATE', 'EMPATICA_BLOOD_VOLUME_PULSE'):
        ddict = readFile(sensor_data_file, sensor)
        df = pd.DataFrame.from_dict(ddict, orient='index', columns=[column])
        df[column] = df[column].astype(float)
        df.index.name = 'timestamp'

    elif sensor == 'EMPATICA_ACCELEROMETER':
        ddict = readFile(sensor_data_file, sensor)
        df = pd.DataFrame.from_dict(ddict, orient='index', columns=['x', 'y', 'z'])
        df['x'] = df['x'].astype(float)
        df['y'] = df['y'].astype(float)
        df['z'] = df['z'].astype(float)
        df.index.name = 'timestamp'

    elif sensor == 'EMPATICA_INTER_BEAT_INTERVAL':
        df = pd.read_csv(sensor_data_file, names=['timestamp', column], header=None)
        timestampstart = float(df['timestamp'][0])
        df['timestamp'] = (df['timestamp'][1:len(df)]).astype(float) + timestampstart
        df = df.drop([0])
        df[column] = df[column].astype(float)
        df = df.set_index('timestamp')

    else:
        raise ValueError(
            "sensor has an invalid name: {}".format(sensor))

    # format timestamps
    df.index *= 1000
    df.index = df.index.astype(int)
    return(df)

def pull_data(data_configuration, device, sensor, container, columns_to_download):
    sensor_csv = container + '.csv'
    warning = True
    participant_data = pd.DataFrame(columns=columns_to_download.values())
    participant_data.set_index('timestamp', inplace=True)

    available_zipfiles = list((Path(data_configuration["FOLDER"]) / Path(device)).rglob("*.zip"))
    if len(available_zipfiles) == 0:
        warnings.warn("There were no zip files in: {}. If you were expecting data for this participant the [EMPATICA][DEVICE_IDS] key in their participant file is missing the pid".format((Path(data_configuration["FOLDER"]) / Path(device))))

    for zipfile in available_zipfiles:
        print("Extracting {} data from {} for {}".format(sensor, zipfile, device))
        with ZipFile(zipfile, 'r') as zipFile:
            listOfFileNames = zipFile.namelist()
            for fileName in listOfFileNames:
                if fileName == sensor_csv:
                    participant_data = pd.concat([participant_data, extract_empatica_data(zipFile.read(fileName),  sensor)], axis=0)
                    warning = False
            if warning:
                warnings.warn("We could not find a zipped file for {} in {} (we tried to find {})".format(sensor, zipFile, sensor_csv))
    
    participant_data.sort_index(inplace=True, ascending=True)
    participant_data.reset_index(inplace=True)
    participant_data.drop_duplicates(subset='timestamp', keep='first',inplace=True)
    participant_data["device_id"] = device
    return(participant_data)

# print(pull_data({'FOLDER': 'data/external/empatica'}, "e01", "EMPATICA_accelerometer", {'TIMESTAMP': 'timestamp', 'DEVICE_ID': 'device_id', 'DOUBLE_VALUES_0': 'x', 'DOUBLE_VALUES_1': 'y', 'DOUBLE_VALUES_2': 'z'}))