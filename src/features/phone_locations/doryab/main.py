import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle
from sklearn.cluster import DBSCAN
from math import radians, cos, sin, asin, sqrt

def doryab_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    location_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    dbscan_eps = provider["DBSCAN_EPS"]
    dbscan_minsamples = provider["DBSCAN_MINSAMPLES"]
    threshold_static = provider["THRESHOLD_STATIC"]
    maximum_gap_allowed = provider["MAXIMUM_GAP_ALLOWED"]
    sampling_frequency = provider["SAMPLING_FREQUENCY"]
    
    minutes_data_used = provider["MINUTES_DATA_USED"]
    if(minutes_data_used):
            requested_features.append("minutesdataused")

    # name of the features this function can compute
    base_features_names = ["locationvariance","loglocationvariance","totaldistance","averagespeed","varspeed","circadianmovement","numberofsignificantplaces","numberlocationtransitions","radiusgyration","timeattop1location","timeattop2location","timeattop3location","movingtostaticratio","outlierstimepercent","maxlengthstayatclusters","minlengthstayatclusters","meanlengthstayatclusters","stdlengthstayatclusters","locationentropy","normalizedlocationentropy","minutesdataused"]    
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    

    if location_data.empty:
        location_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    else:
        location_data = filter_data_by_segment(location_data, time_segment)

        if location_data.empty:
            location_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
        else:
            location_features = pd.DataFrame()

            if sampling_frequency == 0:
                sampling_frequency = getSamplingFrequency(location_data)

            if "minutesdataused" in features_to_compute:
                for localDate in location_data["local_segment"].unique():
                    location_features.loc[localDate,"minutesdataused"] = getMinutesData(location_data[location_data["local_segment"]==localDate])

            location_features.index.name = 'local_segment'
            
            location_data = location_data[(location_data['double_latitude']!=0.0) & (location_data['double_longitude']!=0.0)]

            if location_data.empty:
                location_features = pd.DataFrame(columns=["local_date"] + ["location_" + time_segment + "_" + x for x in features_to_compute])
                location_features = location_features.reset_index(drop=True)
                return location_features

            if "locationvariance" in features_to_compute:
                location_features["locationvariance"] = location_data.groupby(['local_segment'])['double_latitude'].var() + location_data.groupby(['local_segment'])['double_longitude'].var()
            
            if "loglocationvariance" in features_to_compute:
                location_features["loglocationvariance"] = (location_data.groupby(['local_segment'])['double_latitude'].var() + location_data.groupby(['local_segment'])['double_longitude'].var()).apply(lambda x: np.log10(x) if x > 0 else None)

            
            preComputedDistanceandSpeed = pd.DataFrame()
            for localDate in location_data['local_segment'].unique():
                distance, speeddf = get_all_travel_distances_meters_speed(location_data[location_data['local_segment']==localDate],threshold_static,maximum_gap_allowed)
                preComputedDistanceandSpeed.loc[localDate,"distance"] = distance.sum()
                preComputedDistanceandSpeed.loc[localDate,"avgspeed"] = speeddf[speeddf['speedTag'] == 'Moving']['speed'].mean()
                preComputedDistanceandSpeed.loc[localDate,"varspeed"] = speeddf[speeddf['speedTag'] == 'Moving']['speed'].var()

            if "totaldistance" in features_to_compute:
                for localDate in location_data['local_segment'].unique():
                    location_features.loc[localDate,"totaldistance"] = preComputedDistanceandSpeed.loc[localDate,"distance"]

            if "averagespeed" in features_to_compute:
                for localDate in location_data['local_segment'].unique():
                    location_features.loc[localDate,"averagespeed"] = preComputedDistanceandSpeed.loc[localDate,"avgspeed"]

            if "varspeed" in features_to_compute:
                for localDate in location_data['local_segment'].unique():
                    location_features.loc[localDate,"varspeed"] = preComputedDistanceandSpeed.loc[localDate,"varspeed"]

            if "circadianmovement" in features_to_compute:
                for localDate in location_data['local_segment'].unique():
                    location_features.loc[localDate,"circadianmovement"] = circadian_movement(location_data[location_data['local_segment']==localDate])

            newLocationData = cluster_and_label(location_data, eps= distance_to_degrees(dbscan_eps), min_samples=dbscan_minsamples)

            if "numberofsignificantplaces" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"numberofsignificantplaces"] = number_of_significant_places(newLocationData[newLocationData['local_segment']==localDate])

            if "numberlocationtransitions" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"numberlocationtransitions"] = number_location_transitions(newLocationData[newLocationData['local_segment']==localDate])
            
            if "radiusgyration" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"radiusgyration"] = radius_of_gyration(newLocationData[newLocationData['local_segment']==localDate],sampling_frequency)

            if "timeattop1location" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"timeattop1"] = time_at_topn_clusters_in_group(newLocationData[newLocationData['local_segment']==localDate],1,sampling_frequency)

            if "timeattop2location" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"timeattop2"] = time_at_topn_clusters_in_group(newLocationData[newLocationData['local_segment']==localDate],2,sampling_frequency)
            
            if "timeattop3location" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"timeattop3"] = time_at_topn_clusters_in_group(newLocationData[newLocationData['local_segment']==localDate],3,sampling_frequency)

            if "movingtostaticratio" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"movingtostaticratio"] =  (newLocationData[newLocationData['local_segment']==localDate].shape[0]*sampling_frequency) / (location_data[location_data['local_segment']==localDate].shape[0] * sampling_frequency)

            if "outlierstimepercent" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"outlierstimepercent"] = outliers_time_percent(newLocationData[newLocationData['local_segment']==localDate],sampling_frequency)

            preComputedmaxminCluster = pd.DataFrame()
            for localDate in newLocationData['local_segment'].unique():
                    smax, smin, sstd,smean = len_stay_at_clusters_in_minutes(newLocationData[newLocationData['local_segment']==localDate],sampling_frequency)
                    preComputedmaxminCluster.loc[localDate,"maxlengthstayatclusters"] = smax 
                    preComputedmaxminCluster.loc[localDate,"minlengthstayatclusters"] = smin 
                    preComputedmaxminCluster.loc[localDate,"stdlengthstayatclusters"] = sstd
                    preComputedmaxminCluster.loc[localDate,"meanlengthstayatclusters"] = smean
            
            if "maxlengthstayatclusters" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"maxlengthstayatclusters"] = preComputedmaxminCluster.loc[localDate,"maxlengthstayatclusters"]
            
            if "minlengthstayatclusters" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"minlengthstayatclusters"] = preComputedmaxminCluster.loc[localDate,"minlengthstayatclusters"]

            if "stdlengthstayatclusters" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"stdlengthstayatclusters"] = preComputedmaxminCluster.loc[localDate,"stdlengthstayatclusters"]

            if "meanlengthstayatclusters" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"meanlengthstayatclusters"] = preComputedmaxminCluster.loc[localDate,"meanlengthstayatclusters"]

            if "locationentropy" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"locationentropy"] = location_entropy(newLocationData[newLocationData['local_segment']==localDate])

            if "normalizedlocationentropy" in features_to_compute:
                for localDate in newLocationData['local_segment'].unique():
                    location_features.loc[localDate,"normalizedlocationentropy"] = location_entropy_normalized(newLocationData[newLocationData['local_segment']==localDate])
            
            location_features = location_features.reset_index()

    return location_features


def getMinutesData(locationData):

    return locationData[['local_hour','local_minute']].drop_duplicates(inplace = False).shape[0]

def distance_to_degrees(d):
    #Just an approximation, but speeds up clustering by a huge amount and doesnt introduce much error
    #over small distances
    d = d / 1852
    d = d / 60
    return d

def get_all_travel_distances_meters_speed(locationData,threshold,maximum_gap_allowed):
    
    lat_lon_temp = pd.DataFrame()

    lat_lon_temp['_lat_before'] = locationData.double_latitude
    lat_lon_temp['_lat_after'] = locationData.double_latitude.shift(-1)
    lat_lon_temp['_lon_before'] = locationData.double_longitude
    lat_lon_temp['_lon_after'] = locationData.double_longitude.shift(-1)
    lat_lon_temp['time_before'] = pd.to_datetime(locationData['local_time'], format="%H:%M:%S")
    lat_lon_temp['time_after'] = lat_lon_temp['time_before'].shift(-1)
    lat_lon_temp['time_diff'] = lat_lon_temp['time_after'] - lat_lon_temp['time_before']
    lat_lon_temp['timeInSeconds'] = lat_lon_temp['time_diff'].apply(lambda x: x.total_seconds())

    lat_lon_temp = lat_lon_temp[lat_lon_temp['timeInSeconds'] <= maximum_gap_allowed]

    if lat_lon_temp.empty:
        return pd.Series(), pd.DataFrame({"speed": [], "speedTag": []})
    
    lat_lon_temp['distances'] = lat_lon_temp.apply(haversine, axis=1)  # meters
    lat_lon_temp['speed']  = (lat_lon_temp['distances'] / lat_lon_temp['timeInSeconds'] )
    lat_lon_temp['speed'] = lat_lon_temp['speed'].replace(np.inf, np.nan) * 3.6
    distances = lat_lon_temp['distances']
    lat_lon_temp = lat_lon_temp.dropna()
    lat_lon_temp['speedTag'] = np.where(lat_lon_temp['speed'] >= threshold,"Moving","Static")

    return distances,lat_lon_temp[['speed','speedTag']]


def vincenty_row(x):
    """
    :param x: A row from a dataframe
    :return: The distance in meters between
    """

    try:
       return vincenty((x['_lat_before'], x['_lon_before']),(x['_lat_after'], x['_lon_after'])).meters

    except:
       return 0

def haversine(x):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = x['_lon_before'], x['_lat_before'],x['_lon_after'], x['_lat_after']
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r* 1000


def circadian_movement_energies(locationData):
    time = (locationData["timestamp"].values / 1000.0)  # seconds
    ylat = locationData["double_latitude"].values
    ylong = locationData["double_longitude"].values
    hours_intervals = np.arange(23.5, 24.51, 0.01)  # hours
    seconds_intervals = hours_intervals * 60 * 60  # seconds
    frequency = 1 / seconds_intervals

    power_latitude = LombScargle(time, ylat).power(frequency=frequency, normalization='psd')
    power_longitude = LombScargle(time, ylong).power(frequency=frequency, normalization='psd')

    energy_latitude = np.sum(power_latitude)
    energy_longitude = np.sum(power_longitude)
    return (energy_latitude, energy_longitude)

def circadian_movement(locationData):
    
    energy_latitude, energy_longitude = circadian_movement_energies(locationData)
    return np.log10(energy_latitude + energy_longitude)

def cluster_and_label(df,**kwargs):
    """

    :param df:   a df with columns "latitude", "longitude", and "datetime"
                                     or
               a df with comlumns "latitude","longitude" and a datetime index
    :param kwargs: arguments for sklearn's DBSCAN
    :return: a new df of labeled locations with moving points removed, where the cluster
             labeled as "1" is the largest, "2" the second largest, and so on
    """
    location_data = df
    if not isinstance(df.index, pd.DatetimeIndex):
        location_data = df.set_index("local_date_time")

    stationary = remove_moving(location_data,1)

    #return degrees(arcminutes=nautical(meters= d))
    #nautical miles = m ÷ 1,852
    clusterer = DBSCAN(**kwargs)

    counts_df = stationary[["double_latitude" ,"double_longitude"]].groupby(["double_latitude" ,"double_longitude"]).size().reset_index()
    counts = counts_df[0]
    lat_lon = counts_df[["double_latitude","double_longitude"]].values
    cluster_results = clusterer.fit_predict(lat_lon, sample_weight= counts)

    #Need to extend labels back to original df without weights
    counts_df["location_label"] = cluster_results
    # remove the old count column
    del counts_df[0]

    merged = pd.merge(stationary,counts_df, on = ["double_latitude" ,"double_longitude"])

    #Now compute the label mapping:
    cluster_results = merged["location_label"].values
    valid_clusters = cluster_results[np.where(cluster_results != -1)]
    label_map = rank_count_map(valid_clusters)

    #And remap the labels:
    merged.index = stationary.index
    stationary = stationary.assign(location_label = merged["location_label"].map(label_map).values)
    stationary.loc[:, "location_label"] = merged["location_label"].map(label_map)
    return stationary

def rank_count_map(clusters):
    """ Returns a function which will map each element of a list 'l' to its rank,
    such that the most common element maps to 1

    Is used in this context to sort the cluster labels so that cluster with rank 1 is the most
    visited.

    If return_dict, return a mapping dict rather than a function

    If a function, if the value can't be found label as -1

    """
    labels, counts = tuple(np.unique(clusters, return_counts = True))
    sorted_by_count = [x for (y,x) in sorted(zip(counts, labels), reverse = True)]
    label_to_rank = {label : rank + 1 for (label, rank) in [(sorted_by_count[i],i) for i in range(len(sorted_by_count))]}

    return lambda x: label_to_rank.get(x, -1)


def remove_moving(df, v):

    if not df.index.is_monotonic:
        df = df.sort_index()

    lat_lon_temp = pd.DataFrame()

    lat_lon_temp['_lat_before'] = df.double_latitude.shift()
    lat_lon_temp['_lat_after'] =  df.double_latitude.shift(-1)
    lat_lon_temp['_lon_before'] = df.double_longitude.shift()
    lat_lon_temp['_lon_after'] =  df.double_longitude.shift(-1)

    #
    distance = lat_lon_temp.apply( haversine, axis = 1) / 1000
    time = ((pd.to_datetime(df.reset_index().local_date_time.shift(-1),format="%Y-%m-%d %H:%M:%S") - pd.to_datetime(df.reset_index().local_date_time.shift(),format="%Y-%m-%d %H:%M:%S")) / np.timedelta64(1,'s')).fillna(-1) / (60.*60)
    time.index = distance.index.copy()
    
    return df[(distance / time) < v]

def number_of_significant_places(locationData):
    
    uniquelst = locationData[locationData["location_label"] >= 1]["location_label"].unique()
    return len(uniquelst)

def number_location_transitions(locationData):
    
    # ignores transitions from moving to static and vice-versa, but counts transitions from outliers to major location clusters
    df = pd.DataFrame()

    df['boolCol'] = (locationData.location_label == locationData.location_label.shift())

    return df[df['boolCol'] == False].shape[0] - 1

def radius_of_gyration(locationData,sampling_frequency):
    if locationData is None or len(locationData) == 0:
        return None
    # Center is the centroid, not the home location
    valid_clusters = locationData[locationData["location_label"] != -1]
    centroid_all_clusters = (valid_clusters.groupby('location_label')[['double_latitude','double_longitude']].mean()).mean()
    clusters_centroid = valid_clusters.groupby('location_label')[['double_latitude','double_longitude']].mean()
    
    rog = 0
    for labels in clusters_centroid.index:
        lat_lon_dict = dict()
        lat_lon_dict['_lon_before'] = clusters_centroid.loc[labels].double_longitude
        lat_lon_dict['_lat_before'] = clusters_centroid.loc[labels].double_latitude
        lat_lon_dict['_lon_after'] = centroid_all_clusters.double_longitude
        lat_lon_dict['_lat_after'] = centroid_all_clusters.double_latitude
        
        distance = haversine(lat_lon_dict) ** 2
        
        time_in_cluster = locationData[locationData["location_label"]==labels].shape[0]* sampling_frequency
        rog = rog + (time_in_cluster * distance)
    
    time_all_clusters = valid_clusters.shape[0] * sampling_frequency
    if time_all_clusters == 0:
        return 0
    final_rog = (1/time_all_clusters) * rog

    return np.sqrt(final_rog)

def time_at_topn_clusters_in_group(locationData,n,sampling_frequency):  # relevant only for global location features since, top3_clusters = top3_clusters_in_group for local
    
    if locationData is None or len(locationData) == 0:
        return None

    locationData = locationData[locationData["location_label"] >= 1]  # remove outliers/ cluster noise
    valcounts = locationData["location_label"].value_counts().to_dict()
    sorted_valcounts = sorted(valcounts.items(), key=lambda kv: (-kv[1], kv[0]))

    if len(sorted_valcounts) >= n:
        topn = sorted_valcounts[n-1]
        topn_time = topn[1] * sampling_frequency
    else:
        topn_time = None

    return topn_time

def outliers_time_percent(locationData,sampling_frequency):
    if locationData is None or len(locationData) == 0:
        return None
    clusters = locationData["location_label"]
    numoutliers = clusters[(clusters == -1)].sum() * sampling_frequency
    numtotal = len(clusters) * sampling_frequency
    return (float(-1*numoutliers) / numtotal)


def moving_time_percent(locationData):
    if locationData is None or len(locationData) == 0:
        return None
    lbls = locationData["location_label"]
    nummoving = lbls.isnull().sum()
    numtotal = len(lbls)
    # print (nummoving)
    # print(numtotal)
    return (float(nummoving) / numtotal)

def len_stay_at_clusters_in_minutes(locationData,sampling_frequency):
    if locationData is None or len(locationData) == 0:
        return  (None, None, None,None)

    lenstays = []
    count = 0
    prev_loc_label = None
    for row in locationData.iterrows():
        cur_loc_label = row[1]["location_label"]
        if np.isnan(cur_loc_label):
            continue
        elif prev_loc_label == None:
            prev_loc_label = int(cur_loc_label)
            count += 1
        else:
            if prev_loc_label == int(cur_loc_label):
                count += 1
            else:
                lenstays.append(count)
                prev_loc_label = int(cur_loc_label)
                count = 0 + 1
    if count > 0:  # in case of no transition
        lenstays.append(count)
    lenstays = np.array(lenstays) * sampling_frequency

    if len(lenstays) > 0:
        smax = np.max(lenstays)
        smin = np.min(lenstays)
        sstd = np.std(lenstays)
        smean = np.mean(lenstays)
    else:
        smax = None
        smin = None
        sstd = None
        smean = None
    return (smax, smin, sstd, smean)


def location_entropy(locationData):
    if locationData is None or len(locationData) == 0:
        return None

    clusters = locationData[locationData["location_label"] >= 1]  # remove outliers/ cluster noise
    if len(clusters) > 0:
        # Get percentages for each location
        percents = clusters["location_label"].value_counts(normalize=True)
        entropy = -1 * percents.map(lambda x: x * np.log(x)).sum()
        return entropy
    else:
        return None

def location_entropy_normalized(locationData):
    if locationData is None or len(locationData) == 0:
        return None

    locationData = locationData[locationData["location_label"] >= 1]  # remove outliers/ cluster noise
    entropy = location_entropy(locationData)
    unique_clusters = locationData["location_label"].unique()
    num_clusters = len(unique_clusters)
    if num_clusters == 0 or len(locationData) == 0 or entropy is None:
        return None
    else:
        return entropy / num_clusters


def getSamplingFrequency(locationData):

    return ((pd.to_datetime(locationData['local_time'], format="%H:%M:%S") - pd.to_datetime(locationData['local_time'].shift(periods=1), format="%H:%M:%S")).apply(lambda x: x.total_seconds())/60).median()