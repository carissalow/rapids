import pandas as pd
from applications_foreground.applications_foreground_base import base_applications_foreground_features

apps_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"], encoding="ISO-8859-1")
day_segment = snakemake.params["day_segment"]
single_categories = snakemake.params["single_categories"]
multiple_categories_with_genres = snakemake.params["multiple_categories"]
single_apps = snakemake.params["single_apps"]
excluded_categories = snakemake.params["excluded_categories"]
excluded_apps = snakemake.params["excluded_apps"]
requested_features = snakemake.params["features"]
apps_features = pd.DataFrame(columns=["local_date"])

single_categories = list(set(single_categories) - set(excluded_categories))
multiple_categories = list(multiple_categories_with_genres.keys() - set(excluded_categories))
apps = list(set(single_apps) - set(excluded_apps))
type_count = len(single_categories) + len(multiple_categories) + len(apps)

params = {}
params["multiple_categories_with_genres"] = multiple_categories_with_genres
params["single_categories"] = single_categories
params["multiple_categories"] = multiple_categories
params["apps"] = apps

# exclude categories in the excluded_categories list
if "system_apps" in excluded_categories:
    apps_data = apps_data[apps_data["is_system_app"] == 0]
apps_data = apps_data[~apps_data["genre"].isin(excluded_categories)]
# exclude apps in the excluded_apps list
apps_data = apps_data[~apps_data["application_name"].isin(excluded_apps)]

apps_features = apps_features.merge(base_applications_foreground_features(apps_data, day_segment, requested_features, params), on="local_date", how="outer")

assert len(requested_features) * type_count + 1 == apps_features.shape[1], "The number of features in the output dataframe (=" + str(apps_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your application foreground feature extraction functions"

apps_features.to_csv(snakemake.output[0], index=False)
