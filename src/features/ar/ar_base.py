import pandas as pd
import numpy as np
import scipy.stats as stats
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

def base_ar_features(ar_data, ar_deltas, day_segment, requested_features):
    # name of the features this function can compute
    base_features_names = ["count","mostcommonactivity","countuniqueactivities","activitychangecount","sumstationary","summobile","sumvehicle"]
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    ar_features = pd.DataFrame(columns = ["local_date"] + ["ar_" + day_segment + "_" + x for x in features_to_compute])
    if not ar_data.empty:
        ar_deltas = splitOvernightEpisodes(ar_deltas, [],["activity"])

        if day_segment != "daily":
            ar_deltas = splitMultiSegmentEpisodes(ar_deltas, day_segment, [])

        ar_data.local_date_time = pd.to_datetime(ar_data.local_date_time)
        resampledData = ar_data.set_index(ar_data.local_date_time)
        resampledData.drop(columns=["local_date_time"], inplace=True)

        if day_segment != "daily":
            resampledData = resampledData.loc[resampledData["local_day_segment"] == day_segment]

        if not resampledData.empty:
            ar_features = pd.DataFrame()

            # finding the count of samples of the day
            if "count" in features_to_compute:
                ar_features["ar_" + day_segment + "_count"] = resampledData["activity_type"].resample("D").count()

            # finding most common activity of the day
            if "mostcommonactivity" in features_to_compute:
                ar_features["ar_" + day_segment + "_mostcommonactivity"] = resampledData["activity_type"].resample("D").apply(lambda x: stats.mode(x)[0] if len(stats.mode(x)[0]) != 0 else None)

            # finding different number of activities during a day
            if "countuniqueactivities" in features_to_compute:
                ar_features["ar_" + day_segment + "_countuniqueactivities"] = resampledData["activity_type"].resample("D").nunique()
            
            # finding Number of times activity changed
            if "activitychangecount" in features_to_compute:
                resampledData["activity_type_shift"] = resampledData["activity_type"].shift().fillna(resampledData["activity_type"].head(1))
                resampledData["different_activity"] = np.where(resampledData["activity_type"]!=resampledData["activity_type_shift"],1,0)
                ar_features["ar_" + day_segment + "_activitychangecount"] = resampledData["different_activity"].resample("D").sum()


            deltas_features = {"sumstationary":["still","tilting"], 
                            "summobile":["on_foot","walking","running","on_bicycle"],
                            "sumvehicle":["in_vehicle"]}
            
            for column, activity_labels in deltas_features.items():
                if column in features_to_compute:
                    filtered_data = ar_deltas[ar_deltas["activity"].isin(pd.Series(activity_labels))]
                    if not filtered_data.empty:
                        ar_features["ar_" + day_segment + "_" + column] = ar_deltas[ar_deltas["activity"].isin(pd.Series(activity_labels))].groupby(["local_start_date"])["time_diff"].sum().fillna(0)
                    else:
                        ar_features["ar_" + day_segment + "_" + column] = 0
        
            ar_features.index.names = ["local_date"]
            ar_features = ar_features.reset_index()

    return ar_features
