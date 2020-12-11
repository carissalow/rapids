import pandas as pd
import numpy as np
from sklearn.cluster import KMeans

def deviceFeatures(devices, ownership, features_to_compute, features):
    if devices.shape[0] == 0:
        device_value_counts = pd.DataFrame(columns=["local_segment", "bt_address", "scans"], dtype=int)
    else:
        device_value_counts = devices.groupby(["local_segment"])["bt_address"].value_counts().to_frame("scans").reset_index()

    if "countscans" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].sum().to_frame("countscans" + ownership), how="outer")
    if "uniquedevices" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["bt_address"].nunique().to_frame("uniquedevices" + ownership), how="outer")
    if "countscansmostuniquedevice" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].max().to_frame("countscansmostuniquedevice" + ownership), how="outer")
    if "countscansleastuniquedevice" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].min().to_frame("countscansleastuniquedevice" + ownership), how="outer")
    if "meanscans" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].mean().to_frame("meanscans" + ownership), how="outer")
    if "stdscans" in features_to_compute:
        features = features.join(device_value_counts.groupby("local_segment")["scans"].std().to_frame("stdscans" + ownership), how="outer")
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

    # K = 3, devices I own VS devices my partner/roommate owns (can also be other devices I own though) VS devices other people own
    kmeans_k3 = KMeans(n_clusters=3, init=initial_k3,  n_init = 1).fit(X)
    labels_k3 = kmeans_k3.labels_
    centers_k3 = [c[0] for c in kmeans_k3.cluster_centers_]
    diff_k3 = [(X_array[xi] - centers_k3[labels_k3[xi]])**2 for xi in range(0, len(X_array))]
    sum_dist_k3 = sum(diff_k3)

    if sum_dist_k2 < sum_dist_k3: # K = 2 is better
        labels = labels_k2
        centers = centers_k2
        numclust = 2
    else:
        labels = labels_k3
        centers = centers_k3
        numclust = 3
    
    maxcluster = np.where(labels == np.argmax(centers), 1, 0)
    bt_frequency["own_device"] = maxcluster
    return bt_frequency[["bt_address", "own_device"]]
    

def doryab_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    bt_data = pd.read_csv(sensor_data_files["sensor_data"])
    base_features = set(["countscans", "uniquedevices", "countscansmostuniquedevice", "countscansleastuniquedevice", "meanscans", "stdscans"])
    ownership_keys = [x.lower() for x in provider["FEATURES"].keys()]
    if set(ownership_keys) != set(["own", "others", "all"]):
        raise ValueError("[PHONE_BLUETOOTH][DORYAB][FEATURES] config key can only have three lists called ALL, OWN and OTHERS, instead you provided {}".format(ownership_keys))
    
    device_ownership = ownership_based_on_clustering(deviceFrequency(bt_data)).set_index("bt_address")
    bt_data = bt_data.set_index("bt_address").join(device_ownership, how="left").reset_index()
    bt_data["own_device"].fillna(0, inplace=True)
    segment_bt_data = filter_data_by_segment(bt_data, time_segment)
    features = pd.DataFrame(columns=['local_segment']).set_index("local_segment")
    for ownership in provider["FEATURES"].keys():
        features_to_compute = list(set(provider["FEATURES"][ownership]) & base_features)
        if ownership == "OWN":
            owner_segment_bt_data = segment_bt_data.query("own_device == 1")
        elif ownership == "OTHERS":
            owner_segment_bt_data = segment_bt_data.query("own_device == 0")
        else: #ALL
            owner_segment_bt_data = segment_bt_data
        features = deviceFeatures(owner_segment_bt_data, ownership.lower(), features_to_compute, features)
        
    features = features.reset_index()
    return features
