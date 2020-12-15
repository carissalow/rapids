from zipfile import ZipFile
import warnings
sensor_short_name = {"accelerometer":"acc",
                "temperature":"temp",
                "tags":"tags",
                "heartrate":"hr",
                "inter_beat_interval":"ibi",
                "blood_volume_pulse":"bvp",
                "electrodermal_activity":"eda"}

sensor_csv = sensor_short_name[snakemake.params["sensor"]] + '.csv'
warning = True
with ZipFile(snakemake.input[0], 'r') as zipFile:
    listOfFileNames = zipFile.namelist()
    for fileName in listOfFileNames:
        if fileName.endswith(sensor_csv):
            with open(snakemake.output[0], 'wb') as outputFile:
                outputFile.write(zipFile.read(fileName))
                warning = False
if(warning):
    warnings.warn("We could not find a zipped file for {} in {} (we tried to find {})".format(snakemake.params["sensor"], snakemake.input[0], sensor_csv))
