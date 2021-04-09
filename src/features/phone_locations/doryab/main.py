import pandas as pd
import numpy as np
import warnings
from astropy.timeseries import LombScargle
from sklearn.cluster import DBSCAN,OPTICS
from math import radians, cos, sin, asin, sqrt

def doryab_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    location_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    accuracy_limit = provider["ACCURACY_LIMIT"]
    dbscan_eps = provider["DBSCAN_EPS"]
    dbscan_minsamples = provider["DBSCAN_MINSAMPLES"]
    threshold_static = provider["THRESHOLD_STATIC"]
    maximum_gap_allowed = provider["MAXIMUM_ROW_GAP"]
    maximum_row_duration = provider["MAXIMUM_ROW_DURATION"]
    cluster_on = provider["CLUSTER_ON"]
    clustering_algorithm = provider["CLUSTERING_ALGORITHM"]
    radius_from_home = provider["RADIUS_FOR_HOME"]
    
    minutes_data_used = provider["MINUTES_DATA_USED"]
    if(minutes_data_used):
            requested_features.append("minutesdataused")

    # name of the features this function can compute
    base_features_names = ["locationvariance","loglocationvariance","totaldistance","averagespeed","varspeed","circadianmovement","numberofsignificantplaces","numberlocationtransitions","radiusgyration","timeattop1location","timeattop2location","timeattop3location","movingtostaticratio","outlierstimepercent","maxlengthstayatclusters","minlengthstayatclusters","meanlengthstayatclusters","stdlengthstayatclusters","locationentropy","normalizedlocationentropy","minutesdataused","timeathome"]    
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    if clustering_algorithm == "DBSCAN":
        hyperparameters = {'eps' : distance_to_degrees(dbscan_eps), 'min_samples': dbscan_minsamples}
    elif clustering_algorithm == "OPTICS":
        hyperparameters = {'max_eps': distance_to_degrees(dbscan_eps), 'min_samples': 2, 'metric':'euclidean', 'cluster_method' : 'dbscan'} 
    else:
        raise ValueError("config[PHONE_LOCATIONS][DORYAB][CLUSTERING ALGORITHM] only accepts DBSCAN or OPTICS but you provided ",clustering_algorithm)
    
    rows_before_accuracy_filter = len(location_data)
    location_data.query("accuracy < @accuracy_limit", inplace=True)
    if rows_before_accuracy_filter > 0 and len(location_data) == 0:
        warnings.warn("Cannot compute Doryab location features because there are no rows with an accuracy value lower than ACCURACY_LIMIT: {}".format(accuracy_limit))

    if location_data.empty:
        location_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    else:
        if cluster_on == "PARTICIPANT_DATASET":
            location_data = cluster_and_label(location_data,clustering_algorithm,threshold_static,**hyperparameters)
            location_data = filter_data_by_segment(location_data, time_segment)
        elif cluster_on == "TIME_SEGMENT":
            location_data = filter_data_by_segment(location_data, time_segment)
            location_data = cluster_and_label(location_data,clustering_algorithm,threshold_static,**hyperparameters)
        else:
            raise ValueError("config[PHONE_LOCATIONS][DORYAB][CLUSTER_ON] only accepts PARTICIPANT_DATASET or TIME_SEGMENT but you provided ",cluster_on)

        if location_data.empty:
            location_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
        else:
            location_features = pd.DataFrame()

            if "minutesdataused" in features_to_compute:
                for localDate in location_data["local_segment"].unique():
                    location_features.loc[localDate,"minutesdataused"] = getMinutesData(location_data[location_data["local_segment"]==localDate])

            location_features.index.name = 'local_segment'
            
            location_data = location_data[(location_data['double_latitude']!=0.0) & (location_data['double_longitude']!=0.0)]

            if location_data.empty:
                location_features = pd.DataFrame(columns=["local_segment"] + ["location_" + time_segment + "_" + x for x in features_to_compute])
                location_features = location_features.reset_index(drop=True)
                return location_features

            location_data['timeInSeconds'] = (location_data.timestamp.diff(-1)* -1)/1000    
            if "locationvariance" in features_to_compute:
                location_features["locationvariance"] = location_data.groupby(['local_segment'])['double_latitude'].var() + location_data.groupby(['local_segment'])['double_longitude'].var()
            
            if "loglocationvariance" in features_to_compute:
                location_features["loglocationvariance"] = (location_data.groupby(['local_segment'])['double_latitude'].var() + location_data.groupby(['local_segment'])['double_longitude'].var()).apply(lambda x: np.log10(x) if x > 0 else None)

            
            preComputedDistanceandSpeed = pd.DataFrame()
            for localDate in location_data['local_segment'].unique():
                speeddf = get_all_travel_distances_meters_speed(location_data[location_data['local_segment']==localDate],threshold_static,maximum_gap_allowed)
                preComputedDistanceandSpeed.loc[localDate,"distance"] = speeddf['distances'].sum() 
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

            
            stationaryLocations = location_data[location_data['stationary_or_not'] == 1]

            if "numberofsignificantplaces" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"numberofsignificantplaces"] = number_of_significant_places(stationaryLocations[stationaryLocations['local_segment']==localDate])

            if "numberlocationtransitions" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"numberlocationtransitions"] = number_location_transitions(stationaryLocations[stationaryLocations['local_segment']==localDate])
            
            if "radiusgyration" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"radiusgyration"] = radius_of_gyration(stationaryLocations[stationaryLocations['local_segment']==localDate])

            preComputedTimeArray = pd.DataFrame()
            for localDate in stationaryLocations["local_segment"].unique():
                top1,top2,top3,smax,smin,sstd,smean = len_stay_timeattopn(stationaryLocations[stationaryLocations["local_segment"]==localDate],maximum_gap_allowed,maximum_row_duration)
                preComputedTimeArray.loc[localDate,"timeattop1"] = top1
                preComputedTimeArray.loc[localDate,"timeattop2"] = top2
                preComputedTimeArray.loc[localDate,"timeattop3"] = top3
                preComputedTimeArray.loc[localDate,"maxlengthstayatclusters"] = smax
                preComputedTimeArray.loc[localDate,"minlengthstayatclusters"] = smin
                preComputedTimeArray.loc[localDate,"stdlengthstayatclusters"] = sstd
                preComputedTimeArray.loc[localDate,"meanlengthstayatclusters"] = smean

            if "timeattop1location" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"timeattop1"] = preComputedTimeArray.loc[localDate,"timeattop1"]

            if "timeattop2location" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"timeattop2"] = preComputedTimeArray.loc[localDate,"timeattop2"]
            
            if "timeattop3location" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"timeattop3"] = preComputedTimeArray.loc[localDate,"timeattop3"]

            if "movingtostaticratio" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"movingtostaticratio"] =  (stationaryLocations[stationaryLocations['local_segment']==localDate]['timeInSeconds'].sum()) / (location_data[location_data['local_segment']==localDate]['timeInSeconds'].sum())

            if "outlierstimepercent" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"outlierstimepercent"] = outlier_time_percent_new(stationaryLocations[stationaryLocations['local_segment']==localDate])
            
            if "maxlengthstayatclusters" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"maxlengthstayatclusters"] = preComputedTimeArray.loc[localDate,"maxlengthstayatclusters"]
            
            if "minlengthstayatclusters" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"minlengthstayatclusters"] = preComputedTimeArray.loc[localDate,"minlengthstayatclusters"]

            if "stdlengthstayatclusters" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"stdlengthstayatclusters"] = preComputedTimeArray.loc[localDate,"stdlengthstayatclusters"]

            if "meanlengthstayatclusters" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"meanlengthstayatclusters"] = preComputedTimeArray.loc[localDate,"meanlengthstayatclusters"]

            if "locationentropy" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"locationentropy"] = location_entropy(stationaryLocations[stationaryLocations['local_segment']==localDate])

            if "normalizedlocationentropy" in features_to_compute:
                for localDate in stationaryLocations['local_segment'].unique():
                    location_features.loc[localDate,"normalizedlocationentropy"] = location_entropy_normalized(stationaryLocations[stationaryLocations['local_segment']==localDate])
            
            if "timeathome" in features_to_compute:
                calculationDf = stationaryLocations[['local_segment','distancefromhome','timeInSeconds']].copy()
                calculationDf.loc[calculationDf.timeInSeconds >= maximum_gap_allowed,'timeInSeconds'] = maximum_row_duration
                location_features["timeathome"] = calculationDf[calculationDf["distancefromhome"] <= radius_from_home].groupby("local_segment")["timeInSeconds"].sum()/60

            location_features = location_features.reset_index()

    return location_features

def len_stay_timeattopn(locationData,maximum_gap_allowed,maximum_row_duration):
    if locationData is None or len(locationData) == 0:
        return  (None, None, None,None, None, None, None)

    calculationDf = locationData[locationData["location_label"] >= 1][['location_label','timeInSeconds']].copy()
    calculationDf.loc[calculationDf.timeInSeconds >= maximum_gap_allowed,'timeInSeconds'] = maximum_row_duration
    timeArray = calculationDf.groupby('location_label')['timeInSeconds'].sum().reset_index()['timeInSeconds'].sort_values(ascending=False)/60
    
    if len(timeArray) == 3:
        return (timeArray[0],timeArray[1],timeArray[2],timeArray.max(),timeArray.min(),timeArray.std(),timeArray.mean())
    elif len(timeArray)==2:
        return (timeArray[0],timeArray[1],None,timeArray.max(),timeArray.min(),timeArray.std(),timeArray.mean())
    elif len(timeArray)==1:
        return (timeArray[0],None,None,timeArray.max(),timeArray.min(),timeArray.std(),timeArray.mean())
    else:
        return (None,None,None,timeArray.max(),timeArray.min(),timeArray.std(),timeArray.mean())
    

def getMinutesData(locationData):

    return locationData[['local_hour','local_minute']].drop_duplicates(inplace = False).shape[0]

def distance_to_degrees(d):
    #Just an approximation, but speeds up clustering by a huge amount and doesnt introduce much error
    #over small distances
    d = d / 1852
    d = d / 60
    return d

def get_all_travel_distances_meters_speed(locationData,threshold,maximum_gap_allowed):
    
    lat_lon_temp = locationData[locationData['timeInSeconds'] <= maximum_gap_allowed][['double_latitude','double_longitude','timeInSeconds']]
    
    if lat_lon_temp.empty:
        return pd.DataFrame({"speed": [], "speedTag": [],"distances": []})
    
    lat_lon_temp['distances'] = haversine(lat_lon_temp['double_longitude'],lat_lon_temp['double_latitude'],lat_lon_temp['double_longitude'].shift(-1),lat_lon_temp['double_latitude'].shift(-1)) 
    lat_lon_temp['speed']  = (lat_lon_temp['distances'] / lat_lon_temp['timeInSeconds'] )  # meter/second
    lat_lon_temp['speed'] = lat_lon_temp['speed'].replace(np.inf, np.nan) * 3.6  
    
    lat_lon_temp['speedTag'] = np.where(lat_lon_temp['speed'] >= threshold,"Moving","Static")
    
    return lat_lon_temp[['speed','speedTag','distances']]


def vincenty_row(x):
    """
    :param x: A row from a dataframe
    :return: The distance in meters between
    """

    try:
       return vincenty((x['_lat_before'], x['_lon_before']),(x['_lat_after'], x['_lon_after'])).meters

    except:
       return 0

def haversine(lon1,lat1,lon2,lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = np.radians([lon1, lat1, lon2, lat2])

    # haversine formula 
    a = np.sin((lat2-lat1)/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2-lon1)/2.0)**2

    r = 6371 # Radius of earth in kilometers. Use 3956 for miles

    return (r * 2 * np.arcsin(np.sqrt(a)) * 1000)


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

def cluster_and_label(location_data, clustering_algorithm, threshold_static, **kwargs):
    """

    :param location_data:   
        a data frame with "latitude", "longitude", and "datetime" columns
                                     or
        a data frame with "latitude", "longitude", and a datetime index
    :param kwargs: arguments for sklearn's DBSCAN or OPTICS
    :return: 
        a new data frame of labeled locations with moving points removed, where the cluster 
        labeled as "0" is the largest, "1" is the second largest, and so on
    """
    if not location_data.empty:
        if not isinstance(location_data.index, pd.DatetimeIndex):
            location_data = location_data.set_index("local_date_time")

        stationary = mark_moving(location_data, threshold_static)

        counts_df = stationary[["double_latitude", "double_longitude"]].groupby(["double_latitude", "double_longitude"]).size().reset_index()
        counts = counts_df[0]
        lat_lon = counts_df[["double_latitude", "double_longitude"]].values

        if counts_df.shape[0] == 1:
            cluster_results = np.array([-1])
        elif clustering_algorithm == "DBSCAN":
            clusterer = DBSCAN(**kwargs)
            cluster_results = clusterer.fit_predict(lat_lon, sample_weight=counts)
        else:
            clusterer = OPTICS(**kwargs)
            cluster_results = clusterer.fit_predict(lat_lon)

        #Need to extend labels back to original df without weights
        counts_df["location_label"] = cluster_results
        # remove the old count column
        del counts_df[0]

        merged = pd.merge(stationary, counts_df, on=["double_latitude", "double_longitude"])

        #Now compute the label mapping:
        cluster_results = merged["location_label"].values
        valid_clusters = cluster_results[np.where(cluster_results != -1)]
        label_map = rank_count_map(valid_clusters)

        #And remap the labels:
        merged.index = stationary.index
        stationary = stationary.assign(location_label=merged["location_label"].map(label_map).values)
        stationary.loc[:, "location_label"] = merged["location_label"].map(label_map)
        return stationary
    else:
        return location_data

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


def mark_moving(df, threshold_static):

    if not df.index.is_monotonic:
        df = df.sort_index()

    distance = haversine(df.double_longitude,df.double_latitude,df.double_longitude.shift(-1),df.double_latitude.shift(-1))/ 1000
    time = (df.timestamp.diff(-1) * -1) / (1000*60*60)
    
    df['stationary_or_not'] = np.where((distance / time) < threshold_static,1,0)   # 1 being stationary,0 for moving 

    return df

def number_of_significant_places(locationData):
    
    uniquelst = locationData[locationData["location_label"] >= 1]["location_label"].unique()
    return len(uniquelst)

def number_location_transitions(locationData):
    
    # ignores transitions from moving to static and vice-versa, but counts transitions from outliers to major location clusters
    df = pd.DataFrame()

    df['boolCol'] = (locationData.location_label == locationData.location_label.shift())

    return df[df['boolCol'] == False].shape[0] - 1

def radius_of_gyration(locationData):
    if locationData is None or len(locationData) == 0:
        return None
    # Center is the centroid, not the home location
    valid_clusters = locationData[locationData["location_label"] != -1]
    centroid_all_clusters = (valid_clusters.groupby('location_label')[['double_latitude','double_longitude']].mean()).mean()
    clusters_centroid = valid_clusters.groupby('location_label')[['double_latitude','double_longitude']].mean()
    
    rog = 0
    for labels in clusters_centroid.index:    
        distance = haversine(clusters_centroid.loc[labels].double_longitude,clusters_centroid.loc[labels].double_latitude,
                    centroid_all_clusters.double_longitude,centroid_all_clusters.double_latitude) ** 2
        
        time_in_cluster = locationData[locationData["location_label"]==labels]['timeInSeconds'].sum()
        rog = rog + (time_in_cluster * distance)
    
    time_all_clusters = valid_clusters['timeInSeconds'].sum()
    if time_all_clusters == 0:
        return 0
    final_rog = (1/time_all_clusters) * rog

    return np.sqrt(final_rog)

def outlier_time_percent_new(locationData):
    if locationData is None or len(locationData)==0:
        return None

    clustersDf = locationData[["location_label","timeInSeconds"]]
    numoutliers = clustersDf[clustersDf["location_label"]== -1]["timeInSeconds"].sum()
    numtotal = clustersDf.timeInSeconds.sum()

    return numoutliers/numtotal

def location_entropy(locationData):
    if locationData is None or len(locationData) == 0:
        return None

    clusters = locationData[locationData["location_label"] >= 1]  # remove outliers/ cluster noise
    if len(clusters) > 0:
        # Get percentages for each location
        percents = clusters.groupby(['location_label'])['timeInSeconds'].sum() / clusters['timeInSeconds'].sum()
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
    elif np.log(num_clusters)==0:
        return None
    else:
        return entropy / np.log(num_clusters)
