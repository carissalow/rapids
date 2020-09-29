import pandas as pd
from importlib import import_module, util
from pathlib import Path

# import fetch_provider_features from src/features/utils/utils.py
spec = util.spec_from_file_location("util", str(Path(snakemake.scriptdir).parent / "utils" / "utils.py"))
mod = util.module_from_spec(spec)
spec.loader.exec_module(mod)
fetch_provider_features = getattr(mod,  "fetch_provider_features")

battery_episodes_file = snakemake.input["battery_episodes"]
day_segments_file = snakemake.input["day_segments_labels"]
provider = snakemake.params["provider"]
provider_key = snakemake.params["provider_key"]

sensor_features = fetch_provider_features(provider, provider_key, "battery", battery_episodes_file, day_segments_file)

sensor_features.to_csv(snakemake.output[0], index=False)