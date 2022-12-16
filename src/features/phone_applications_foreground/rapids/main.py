import pandas as pd
import numpy as np
import itertools
from scipy.stats import entropy


def compute_features(filtered_data, apps_type, requested_features, apps_features, time_segment):        
    if "timestamp" in filtered_data.columns:
        timestamp_column = "timestamp"
    else:
        timestamp_column = "start_timestamp"

    # There is the rare occasion that filtered_data is empty (found in testing)
    if "timeoffirstuse" in requested_features:
        time_first_event = filtered_data.sort_values(by=timestamp_column, ascending=True).drop_duplicates(subset="local_segment", keep="first").set_index("local_segment")
        if time_first_event.empty:
            apps_features["timeoffirstuse" + apps_type] = np.nan
        elif timestamp_column == "timestamp":
            apps_features["timeoffirstuse" + apps_type] = time_first_event["local_hour"] * 60 + time_first_event["local_minute"]
        else:
            apps_features["timeoffirstuse" + apps_type] = time_first_event["local_start_date_time"].dt.hour * 60 + time_first_event["local_start_date_time"].dt.minute

    if "timeoflastuse" in requested_features:
        time_last_event = filtered_data.sort_values(by=timestamp_column, ascending=False).drop_duplicates(subset="local_segment", keep="first").set_index("local_segment")
        if time_last_event.empty:
            apps_features["timeoflastuse" + apps_type] = np.nan
        elif timestamp_column == "timestamp":
            apps_features["timeoflastuse" + apps_type] = time_last_event["local_hour"] * 60 + time_last_event["local_minute"]
        else:
            apps_features["timeoflastuse" + apps_type] = time_last_event["local_start_date_time"].dt.hour * 60 + time_last_event["local_start_date_time"].dt.minute

    if "frequencyentropy" in requested_features:
        apps_with_count = filtered_data.groupby(["local_segment","application_name"]).count().sort_values(by=timestamp_column, ascending=False).reset_index()
        if (len(apps_with_count.index) < 2 ):
            apps_features["frequencyentropy" + apps_type] = np.nan
        else:    
            apps_features["frequencyentropy" + apps_type] = apps_with_count.groupby("local_segment")[timestamp_column].agg(entropy)

    if "countevent" in requested_features:
        apps_features["countevent" + apps_type] = filtered_data.groupby(["local_segment"]).count()[timestamp_column]

    if "countepisode" in requested_features:
        apps_features["countepisode" + apps_type] = filtered_data.groupby(["local_segment"]).count()["start_timestamp"]

    if "minduration" in requested_features:
        apps_features["minduration" + apps_type] = filtered_data.groupby(by = ["local_segment"])["duration"].min()
            
    if "maxduration" in requested_features:
        apps_features["maxduration" + apps_type] = filtered_data.groupby(by = ["local_segment"])["duration"].max()
            
    if "meanduration" in requested_features:
        apps_features["meanduration" + apps_type] = filtered_data.groupby(by = ["local_segment"])["duration"].mean()
            
    if "sumduration" in requested_features:
        apps_features["sumduration" + apps_type] = filtered_data.groupby(by = ["local_segment"])["duration"].sum()
    
    apps_features.index.names = ["local_segment"]
    return apps_features

def process_app_features(data, requested_features, time_segment, provider, filter_data_by_segment):
    
    excluded_categories = provider["EXCLUDED_CATEGORIES"]
    excluded_apps = provider["EXCLUDED_APPS"]
    single_categories = provider["SINGLE_CATEGORIES"]
    multiple_categories = {}
    if isinstance(provider["MULTIPLE_CATEGORIES"], dict):
        for mcategory_name, mcategory_content in provider["MULTIPLE_CATEGORIES"].items():
            if len(mcategory_content) > 0 and mcategory_name not in excluded_categories:
                multiple_categories[mcategory_name] = mcategory_content
    custom_categories = {}
    if isinstance(provider["CUSTOM_CATEGORIES"], dict):
        for owncategory_name, owncategory_content in provider["CUSTOM_CATEGORIES"].items():
            if len(owncategory_content) > 0 and owncategory_name not in excluded_categories:
                custom_categories[owncategory_name] = owncategory_content
    single_apps = provider["SINGLE_APPS"]
    single_categories = list(set(single_categories) - set(excluded_categories))
    single_apps = list(set(single_apps) - set(excluded_apps))

    # exclude categories in the excluded_categories list
    if "system_apps" in excluded_categories:
        data = data[data["is_system_app"] == 0]
    data = data[~data["genre"].isin(excluded_categories)]
    # exclude apps in the excluded_apps list
    data = data[~data["package_name"].isin(excluded_apps)]
            
    features = pd.DataFrame(columns=["local_segment"] + ["".join(feature) for feature in itertools.product(requested_features, single_categories + list(custom_categories.keys()) + list(multiple_categories.keys()) + single_apps)])
    if not data.empty:
        # deep copy the data for the top1global computation
        data_global = data.copy()
        
        data = filter_data_by_segment(data, time_segment)

        if not data.empty:
            features = pd.DataFrame()
            # single category
            single_categories.sort()
            for sc in single_categories:
                if sc == "all":
                    features = compute_features(data, "all", requested_features, features, time_segment)
                else:
                    filtered_data = data[data["genre"].isin([sc])]
                    features = compute_features(filtered_data, sc, requested_features, features, time_segment)
            # own categories
            for owncategory_name, owncategory_content in custom_categories.items():
                filtered_data = data[data["package_name"].isin(owncategory_content)]
                features = compute_features(filtered_data, owncategory_name, requested_features, features, time_segment)
            # multiple categories
            for mcategory_name, mcategory_content in multiple_categories.items():
                filtered_data = data[data["genre"].isin(mcategory_content)]
                features = compute_features(filtered_data, mcategory_name, requested_features, features, time_segment)
            # single apps
            for app in single_apps:
                col_name = app
                if app == "top1global":
                    # get the most used app
                    apps_with_count = data_global.groupby(["package_name"]).count().sort_values(by="timestamp", ascending=False).reset_index()
                    app = apps_with_count.iloc[0]["package_name"]
                    col_name = "top1global"
                filtered_data = data[data["package_name"].isin([app])]
                features = compute_features(filtered_data, col_name, requested_features, features, time_segment)
 
            features = features.reset_index()

    return features

def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    
    app_foreground_data = pd.read_csv(sensor_data_files["sensor_data"])
    app_episodes_requirement = provider["INCLUDE_EPISODE_FEATURES"]
    
    # if INCLUDE_EPISODE_FEATURES = True, we compute all requested episodes and events features using episode data; otherwise, we compute only events features using event data
    if app_episodes_requirement:
        app_foreground_data = app_foreground_data.drop(app_foreground_data[ (app_foreground_data['duration'] < provider["IGNORE_EPISODES_SHORTER_THAN"]) | (app_foreground_data['duration'] > provider["IGNORE_EPISODES_LONGER_THAN"])].index)
        requested_features = provider["FEATURES"]["APP_EPISODES"] + provider["FEATURES"]["APP_EVENTS"]
    else:
        requested_features = provider["FEATURES"]["APP_EVENTS"]
        
    features = process_app_features(app_foreground_data, requested_features, time_segment, provider, filter_data_by_segment)
    features.fillna(value={feature_name: 0 for feature_name in features.columns if feature_name.startswith(("countevent", "countepisode", "minduration", "maxduration", "meanduration", "sumduration"))}, inplace=True)
    
    return features