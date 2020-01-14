import pandas as pd
import numpy as np

light_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
metrics = snakemake.params["metrics"]

light_features = pd.DataFrame(columns=["local_date"] + ["light_" + day_segment + "_" + x for x in metrics])
if not light_data.empty:
    if day_segment != "daily":
        light_data =light_data[light_data["local_day_segment"] == day_segment]
    
    if not light_data.empty:
        light_features = pd.DataFrame()
        if "count" in metrics:
            light_features["light_" + day_segment + "_count"] = light_data.groupby(["local_date"]).count()["timestamp"]
        
        # get light ambient luminance related features
        if "maxlux" in metrics:
            light_features["light_" + day_segment + "_maxlux"] = light_data.groupby(["local_date"])["double_light_lux"].max()
        if "minlux" in metrics:
            light_features["light_" + day_segment + "_minlux"] = light_data.groupby(["local_date"])["double_light_lux"].min()
        if "avglux" in metrics:
            light_features["light_" + day_segment + "_avglux"] = light_data.groupby(["local_date"])["double_light_lux"].mean()
        if "medianlux" in metrics:
            light_features["light_" + day_segment + "_medianlux"] = light_data.groupby(["local_date"])["double_light_lux"].median()
        if "stdlux" in metrics:
            light_features["light_" + day_segment + "_stdlux"] = light_data.groupby(["local_date"])["double_light_lux"].std()
        
        light_features = light_features.fillna(0).reset_index()

light_features.to_csv(snakemake.output[0], index=False)