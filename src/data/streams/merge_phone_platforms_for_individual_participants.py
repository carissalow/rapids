import pandas as pd
import numpy as np

pid = snakemake.params["pid"]

merged_platforms = pd.DataFrame()
for platforms_file_path in snakemake.input["platforms_files"]:
    platforms_per_sensor = pd.read_csv(platforms_file_path)
    if merged_platforms.shape[0] < platforms_per_sensor.shape[0]:
        merged_platforms = platforms_per_sensor
    elif merged_platforms.shape[0] == platforms_per_sensor.shape[0]:
        if (merged_platforms["os"] == platforms_per_sensor["os"]).all():
            merged_platforms["timestamp"] = np.minimum(merged_platforms["timestamp"], platforms_per_sensor["timestamp"])
        else:
            raise ValueError("Phone platforms of " + pid + " inferred from different sensors are not consistent. Please check it manually!")

merged_platforms.to_csv(snakemake.output[0], index=False)
