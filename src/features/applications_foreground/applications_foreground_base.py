import pandas as pd
import itertools
from scipy.stats import entropy


def compute_features(filtered_data, apps_type, requested_features, apps_features, day_segment):        
    if "timeoffirstuse" in requested_features:
        time_first_event = filtered_data.sort_values(by="timestamp", ascending=True).drop_duplicates(subset="local_date", keep="first").set_index("local_date")
        apps_features["apps_" + day_segment + "_timeoffirstuse" + apps_type] = time_first_event["local_hour"] * 60 + time_first_event["local_minute"]
    if "timeoflastuse" in requested_features:
        time_last_event = filtered_data.sort_values(by="timestamp", ascending=False).drop_duplicates(subset="local_date", keep="first").set_index("local_date")
        apps_features["apps_" + day_segment + "_timeoflastuse" + apps_type] = time_last_event["local_hour"] * 60 + time_last_event["local_minute"]
    if "frequencyentropy" in requested_features:
        apps_with_count = filtered_data.groupby(["local_date","application_name"]).count().sort_values(by="timestamp", ascending=False).reset_index()
        apps_features["apps_" + day_segment + "_frequencyentropy" + apps_type] = apps_with_count.groupby("local_date")["timestamp"].agg(entropy)
    if "count" in requested_features:
        apps_features["apps_" + day_segment + "_count" + apps_type] = filtered_data.groupby(["local_date"]).count()["timestamp"]
        apps_features.fillna(value={"apps_" + day_segment + "_count" + apps_type: 0}, inplace=True)
    return apps_features


def base_applications_foreground_features(apps_data, day_segment, requested_features, params):
    multiple_categories_with_genres = params["multiple_categories_with_genres"]
    single_categories = params["single_categories"]
    multiple_categories = params["multiple_categories"]
    apps = params["apps"]

    # deep copy the apps_data for the top1global computation
    apps_data_global = apps_data.copy()

    apps_features = pd.DataFrame(columns=["local_date"] + ["apps_" + day_segment + "_" + x for x in ["".join(feature) for feature in itertools.product(requested_features, single_categories + multiple_categories + apps)]])
    if not apps_data.empty:
        if day_segment != "daily":
            apps_data =apps_data[apps_data["local_day_segment"] == day_segment]
        
        if not apps_data.empty:
            apps_features = pd.DataFrame()
            # single category
            for sc in single_categories:
                if sc == "all":
                    apps_features = compute_features(apps_data, "all", requested_features, apps_features, day_segment)
                else:
                    filtered_data = apps_data[apps_data["genre"].isin([sc])]
                    apps_features = compute_features(filtered_data, sc, requested_features, apps_features, day_segment)
            # multiple category
            for mc in multiple_categories:
                filtered_data = apps_data[apps_data["genre"].isin(multiple_categories_with_genres[mc])]
                apps_features = compute_features(filtered_data, mc, requested_features, apps_features, day_segment)
            # single apps
            for app in apps:
                col_name = app
                if app == "top1global":
                    # get the most used app
                    apps_with_count = apps_data_global.groupby(["local_date","package_name"]).count().sort_values(by="timestamp", ascending=False).reset_index()
                    app = apps_with_count.iloc[0]["package_name"]
                    col_name = "top1global"
                
                filtered_data = apps_data[apps_data["package_name"].isin([app])]
                apps_features = compute_features(filtered_data, col_name, requested_features, apps_features, day_segment)

            apps_features = apps_features.reset_index()
    
    return apps_features
