import pandas as pd
import numpy as np

pid = snakemake.params["pid"]

merged_platforms = pd.DataFrame(columns=["timestamp", "device_id", "os"])
for platforms_file_path in snakemake.input["platforms_files"]:
    platforms_per_sensor = pd.read_csv(platforms_file_path)
    # Unable to merge automatically if a device_id was used multiple times (e.g. device_id#1, device_id#2, device_id#1)
    if platforms_per_sensor["device_id"].nunique() != platforms_per_sensor.shape[0]:
        raise ValueError("Due to the complexity of " + pid + "'s devices, please merge platforms across all sensors manually!")
    # Merge platforms across all sensors based on device_id
    merged_platforms = merged_platforms.merge(platforms_per_sensor, on=["device_id", "os"], how="outer", sort=False)
    # Keep the smaller timestamp per device_id
    merged_platforms.insert(0, "timestamp", merged_platforms[["timestamp_x", "timestamp_y"]].min(axis=1).astype(int))
    merged_platforms.drop(["timestamp_x", "timestamp_y"], axis=1, inplace=True)

# Sort by timestamp, and drop consecutive duplicates of os column
merged_platforms.sort_values(by=["timestamp"], inplace=True)
merged_platforms = merged_platforms[merged_platforms["os"] != merged_platforms["os"].shift()]

merged_platforms.to_csv(snakemake.output[0], index=False)
