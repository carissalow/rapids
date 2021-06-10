import pandas as pd
import numpy as np
import itertools
from scipy.stats import entropy


def compute_features(filtered_data, apps_type, requested_features, apps_features, time_segment):        
    # There is the rare occasion that filtered_data is empty (found in testing)
    if "timeoffirstuse" in requested_features:
        time_first_event = filtered_data.sort_values(by="timestamp", ascending=True).drop_duplicates(subset="local_segment", keep="first").set_index("local_segment")
        if time_first_event.empty:
            apps_features["timeoffirstuse" + apps_type] = np.nan
        else:
            apps_features["timeoffirstuse" + apps_type] = time_first_event["local_hour"] * 60 + time_first_event["local_minute"]
    if "timeoflastuse" in requested_features:
        time_last_event = filtered_data.sort_values(by="timestamp", ascending=False).drop_duplicates(subset="local_segment", keep="first").set_index("local_segment")
        if time_last_event.empty:
            apps_features["timeoflastuse" + apps_type] = np.nan
        else:
            apps_features["timeoflastuse" + apps_type] = time_last_event["local_hour"] * 60 + time_last_event["local_minute"]
    if "frequencyentropy" in requested_features:
        apps_with_count = filtered_data.groupby(["local_segment","application_name"]).count().sort_values(by="timestamp", ascending=False).reset_index()
        if (len(apps_with_count.index) < 2 ):
            apps_features["frequencyentropy" + apps_type] = np.nan
        else:    
            apps_features["frequencyentropy" + apps_type] = apps_with_count.groupby("local_segment")["timestamp"].agg(entropy)
    if "count" in requested_features:
        apps_features["count" + apps_type] = filtered_data.groupby(["local_segment"]).count()["timestamp"]
        apps_features.fillna(value={"count" + apps_type: 0}, inplace=True)

    if "minduration" in requested_features:
        grouped_data = filtered_data.groupby(by = ['local_segment'])['duration'].min()
        if grouped_data.empty:
            apps_features["minduration" + apps_type] = np.nan
        else:
            apps_features["minduration" + apps_type] = grouped_data
            
    if "maxduration" in requested_features:
        grouped_data = filtered_data.groupby(by = ['local_segment'])['duration'].max()
        if grouped_data.empty:
            apps_features["maxduration" + apps_type] = np.nan
        else:
            apps_features["maxduration" + apps_type] = grouped_data
            
    if "meanduration" in requested_features:
        grouped_data = filtered_data.groupby(by = ['local_segment'])['duration'].mean()
        if grouped_data.empty:
            apps_features["meanduration" + apps_type] = np.nan
        else:
            apps_features["meanduration" + apps_type] = grouped_data
            
    if "sumduration" in requested_features:
        grouped_data = filtered_data.groupby(by = ['local_segment'])['duration'].sum()
        if grouped_data.empty:
            apps_features["sumduration" + apps_type] = np.nan
        else:
            apps_features["sumduration" + apps_type] = grouped_data

    return apps_features

def process_app_features(data, requested_features, time_segment, provider, filter_data_by_segment):
    
    excluded_categories = provider["EXCLUDED_CATEGORIES"]
    excluded_apps = provider["EXCLUDED_APPS"]
    multiple_categories_with_genres = provider["MULTIPLE_CATEGORIES"]
    single_categories = provider["SINGLE_CATEGORIES"]
    multiple_categories = provider["MULTIPLE_CATEGORIES"]
    single_apps = provider["SINGLE_APPS"]

    single_categories = list(set(single_categories) - set(excluded_categories))
    multiple_categories = list(multiple_categories_with_genres.keys() - set(excluded_categories))
    single_apps = list(set(single_apps) - set(excluded_apps))

    # exclude categories in the excluded_categories list
    if "system_apps" in excluded_categories:
        data = data[data["is_system_app"] == 0]
    data = data[~data["genre"].isin(excluded_categories)]
    # exclude apps in the excluded_apps list
    data = data[~data["package_name"].isin(excluded_apps)]
            
    features = pd.DataFrame(columns=["local_segment"] + ["".join(feature) for feature in itertools.product(requested_features, single_categories + multiple_categories + single_apps)])
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
            # multiple category
            for mc in multiple_categories:
                filtered_data = data[data["genre"].isin(multiple_categories_with_genres[mc])]
                features = compute_features(filtered_data, mc, requested_features, features, time_segment)
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
    
    apps_events_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_events_features = provider["FEATURES"]["APP_EVENTS"]
    
    app_episodes_requirement = provider["INCLUDE_EPISODE_FEATURES"]
    
    features = process_app_features(apps_events_data, requested_events_features, time_segment, provider, filter_data_by_segment)
    
    if app_episodes_requirement:
        episode_data = pd.read_csv(sensor_data_files["episode_data"])
        requested_episodes_features = provider["FEATURES"]["APP_EPISODES"]

        episode_data = episode_data.drop(episode_data[ (episode_data['duration'] < provider["IGNORE_EPISODES_SHORTER_THAN"]) & (episode_data['duration'] > provider["IGNORE_EPISODES_LONGER_THAN"])].index)
        
        episodes_features = process_app_features(episode_data, requested_episodes_features, time_segment, provider, filter_data_by_segment)
        
        features = pd.merge(episodes_features, features, how='outer', on='local_segment')
        
    return features