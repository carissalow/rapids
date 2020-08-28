import pandas as pd
from importlib import import_module, util
from pathlib import Path

# import filter_data_by_segment from src/features/utils/utils.py
spec = util.spec_from_file_location("util", str(Path(snakemake.scriptdir).parent / "utils" / "utils.py"))
mod = util.module_from_spec(spec)
spec.loader.exec_module(mod)
filter_data_by_segment = getattr(mod,  "filter_data_by_segment")
rapids_log_tag = getattr(mod,  "rapids_log_tag")

location_data = pd.read_csv(snakemake.input["location_data"][0])
day_segments_labels = pd.read_csv(snakemake.input["day_segments_labels"], header=0)
mypath = snakemake.params["mypath"]
provider = snakemake.params["provider"]
provider_key = snakemake.params["provider_key"]
location_features = pd.DataFrame(columns=["local_segment"])

if "FEATURES" not in provider:
        raise ValueError("Provider config[LOCATION][PROVIDERS][{}] is missing a FEATURES attribute in config.yaml".format(provider_key))

if provider["COMPUTE"] == True:
        code_path = provider["SRC_FOLDER"] + ".main"
        feature_module = import_module(code_path)
        feature_function = getattr(feature_module,  provider["SRC_FOLDER"] + "_location_features")
        
        for day_segment in day_segments_labels["label"]:
                print("{} Processing {} {}".format(rapids_log_tag, provider_key, day_segment))
                features = feature_function(location_data, day_segment, provider, filter_data_by_segment=filter_data_by_segment)
                location_features = location_features.merge(features, how="outer")
else:
        for feature in provider["FEATURES"]:
                location_features[feature] = None

segment_colums = pd.DataFrame()
segment_colums[["local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"]] = location_features["local_segment"].str.split(pat="#", expand=True)
for i in range(segment_colums.shape[1]):
        location_features.insert(1 + i, segment_colums.columns[i], segment_colums[segment_colums.columns[i]])
location_features.to_csv(snakemake.output[0], index=False)