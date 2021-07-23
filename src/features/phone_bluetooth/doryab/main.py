import pandas as pd
import numpy as np
from sklearn.cluster import KMeans

def deviceFeatures(devices, ownership, common_devices, features_to_compute, features):
    if devices.shape[0] == 0:
        device_value_counts = pd.DataFrame(columns=["local_segment", "bt_address", "scans"], dtype=int)
    else:
        device_value_counts = devices.groupby(["local_segment"])["bt_address"].value_counts().to_frame("scans").reset_index()

    if "countscans" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].sum().to_frame("countscans" + ownership), how="outer")
    if "uniquedevices" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["bt_address"].nunique().to_frame("uniquedevices" + ownership), how="outer")
    if "meanscans" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].mean().to_frame("meanscans" + ownership), how="outer")
    if "stdscans" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].std().to_frame("stdscans" + ownership), how="outer")
    # Most frequent device within segments, across segments, and across dataset
    if "countscansmostfrequentdevicewithinsegments" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].max().to_frame("countscansmostfrequentdevicewithinsegments" + ownership), how="outer")
    if "countscansmostfrequentdeviceacrosssegments" in features_to_compute:
        common_device = common_devices['most_segments']
        features = features.join(device_value_counts.query("bt_address in @common_device").groupby("local_segment")["scans"].max().to_frame("countscansmostfrequentdeviceacrosssegments" + ownership), how="outer")
    if "countscansmostfrequentdeviceacrossdataset" in features_to_compute:
        common_device = common_devices['most_dataset']
        features = features.join(device_value_counts.query("bt_address in @common_device").groupby("local_segment")["scans"].max().to_frame("countscansmostfrequentdeviceacrossdataset" + ownership), how="outer")
    # Least frequent device within segments, across segments, and across dataset
    if "countscansleastfrequentdevicewithinsegments" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].min().to_frame("countscansleastfrequentdevicewithinsegments" + ownership), how="outer")
    if "countscansleastfrequentdeviceacrosssegments" in features_to_compute:
        common_device = common_devices['least_segments']
        features = features.join(device_value_counts.query("bt_address in @common_device").groupby("local_segment")["scans"].min().to_frame("countscansleastfrequentdeviceacrosssegments" + ownership), how="outer")
    if "countscansleastfrequentdeviceacrossdataset" in features_to_compute:
        common_device = common_devices['least_dataset']
        features = features.join(device_value_counts.query("bt_address in @common_device").groupby("local_segment")["scans"].min().to_frame("countscansleastfrequentdeviceacrossdataset" + ownership), how="outer")

    return(features)

def deviceFrequency(bt_data):
    bt_data = bt_data[["local_date", "bt_address"]].dropna(subset=["bt_address"])
    bt_data = bt_data.groupby("bt_address").agg({"local_date": pd.Series.nunique, "bt_address" : 'count'})
    bt_data = bt_data.rename(columns={"local_date" : "days_scanned", "bt_address" : "scans"})
    bt_data["scans_per_day"] = bt_data["scans"] / bt_data["days_scanned"]
    return bt_data

def ownership_based_on_clustering(bt_frequency):
    bt_frequency = bt_frequency.reset_index()
    for col in ["scans_per_day", "days_scanned", "scans"]:
        col_zscore = col + '_z'
        bt_frequency[col_zscore] = (bt_frequency[col] - bt_frequency[col].mean()) / bt_frequency[col].std(ddof=0)

    bt_frequency = bt_frequency.dropna(how='any')
    if len(bt_frequency) == 0:
        bt_frequency["own_device"] = None
        return bt_frequency[["bt_address", "own_device"]]

    avgfreq_z = bt_frequency["scans_per_day_z"]
    numdays_z = bt_frequency["days_scanned_z"]
    score = avgfreq_z + numdays_z
    maxscore = np.max(score)
    minscore = np.min(score)
    midscore = (maxscore + minscore) / 2
    initial_k2 = np.array([[maxscore], [minscore]], np.int32)
    initial_k3 = np.array([[maxscore], [midscore], [minscore]], np.int32)
    X_array = score.values
    X = np.reshape(X_array, (len(score), 1))

    # K = 2, devices I own VS devices other people own
    kmeans_k2 = KMeans(n_clusters=2, init = initial_k2, n_init = 1).fit(X)
    labels_k2 = kmeans_k2.labels_
    centers_k2 = [c[0] for c in kmeans_k2.cluster_centers_]
    diff_k2 = [(X_array[xi] - centers_k2[labels_k2[xi]])**2 for xi in range(0, len(X_array))]
    sum_dist_k2 = sum(diff_k2)

    # By default, model with K = 2 is chosen
    labels = labels_k2
    centers = centers_k2
    numclust = 2
    if len(X_array) > 2:
        # K = 3, devices I own VS devices my partner/roommate owns (can also be other devices I own though) VS devices other people own
        kmeans_k3 = KMeans(n_clusters=3, init=initial_k3,  n_init = 1).fit(X)
        labels_k3 = kmeans_k3.labels_
        centers_k3 = [c[0] for c in kmeans_k3.cluster_centers_]
        diff_k3 = [(X_array[xi] - centers_k3[labels_k3[xi]])**2 for xi in range(0, len(X_array))]
        sum_dist_k3 = sum(diff_k3)
        # Model with K = 3 is chosen if sum of squared distances between clustered points and cluster centers is smaller or equal to what we get with K = 2
        if sum_dist_k3 <= sum_dist_k2:
            labels = labels_k3
            centers = centers_k3
            numclust = 3
    
    maxcluster = np.where(labels == np.argmax(centers), 1, 0)
    bt_frequency["own_device"] = maxcluster
    return bt_frequency[["bt_address", "own_device"]]

def mostLeastScannedDevices(devices):
    device_counts = devices["bt_address"].value_counts().sort_index(ascending=False).sort_values(ascending=False)
    return ("","") if (len(device_counts) == 0) else (device_counts.idxmax(), device_counts.idxmin())

def validate_requested_features(provider):
    base_features = {"DEVICES": set(["countscans", "uniquedevices", "meanscans", "stdscans"]),
                        "SCANS_MOST_FREQUENT_DEVICE": set(["withinsegments", "acrosssegments", "acrossdataset"]),
                        "SCANS_LEAST_FREQUENT_DEVICE": set(["withinsegments", "acrosssegments", "acrossdataset"])}

    # Check we have three arrays of features
    ownership_keys = [x.lower() for x in provider["FEATURES"].keys()]
    if set(ownership_keys) != set(["own", "others", "all"]):
        raise ValueError("[PHONE_BLUETOOTH][DORYAB][FEATURES] config key must have three types called ALL, OWN and OTHERS, instead you provided {}".format(ownership_keys))
    
    # Check each array contains valid features
    for ownership_key in provider["FEATURES"].keys():
        for type_key in provider["FEATURES"][ownership_key]:
            if len(provider["FEATURES"][ownership_key][type_key]) > 0 and not set(provider["FEATURES"][ownership_key][type_key]) <= base_features[type_key]:
                raise ValueError("[PHONE_BLUETOOTH][DORYAB][FEATURES][{}][{}] config key only supports features called [{}], instead you provided [{}]".format(ownership_key, type_key, ",".join(base_features[type_key]), ",".join(set(provider["FEATURES"][ownership_key][type_key]) - base_features[type_key])))

def doryab_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    bt_data = pd.read_csv(sensor_data_files["sensor_data"])
    feature_prefix = {"DEVICES":"", "SCANS_MOST_FREQUENT_DEVICE":"countscansmostfrequentdevice", "SCANS_LEAST_FREQUENT_DEVICE":"countscansleastfrequentdevice"}
    validate_requested_features(provider)

    device_ownership = ownership_based_on_clustering(deviceFrequency(bt_data))
    bt_data = bt_data.merge(device_ownership, how="left", on="bt_address")
    bt_data["own_device"].fillna(0, inplace=True)
    dataset_most_common_device, dataset_least_common_device = mostLeastScannedDevices(bt_data)
    segment_bt_data = filter_data_by_segment(bt_data, time_segment)  
    features = pd.DataFrame(columns=['local_segment']).set_index("local_segment")
    for ownership in provider["FEATURES"].keys():

        features_to_compute = []
        for type_key in provider["FEATURES"][ownership]:
            features_to_compute = features_to_compute + [feature_prefix[type_key] + feature for feature in provider["FEATURES"][ownership][type_key]]

        if ownership == "OWN":
            owner_segment_bt_data = segment_bt_data.query("own_device == 1")
        elif ownership == "OTHERS":
            owner_segment_bt_data = segment_bt_data.query("own_device == 0")
        else: #ALL
            owner_segment_bt_data = segment_bt_data

        segment_most_common_device, segment_least_common_device = mostLeastScannedDevices(owner_segment_bt_data)
        common_devices = {"most_dataset": dataset_most_common_device, "least_dataset": dataset_least_common_device, 
                        "most_segments": segment_most_common_device, "least_segments": segment_least_common_device}

        features = deviceFeatures(owner_segment_bt_data, ownership.lower(), common_devices, features_to_compute, features)
    features = features.reset_index()

    # Impute all NaN except for std dev
    for column in features:
        if column not in ["stdscansall", "stdscansown", "stdscansothers"]:
            features[column].fillna(0.0, inplace=True)
    return features
