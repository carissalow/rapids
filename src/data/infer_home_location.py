import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN,OPTICS
from math import radians, cos, sin, asin, sqrt

def filterDatafromDf(origDf):
    
    return origDf[origDf['local_hour']<=6]

def distance_to_degrees(d):
    #Just an approximation, but speeds up clustering by a huge amount and doesnt introduce much error
    #over small distances
    d = d / 1852
    d = d / 60
    return d

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

        merged = pd.merge(stationary,counts_df, on=["double_latitude" ,"double_longitude"])

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

# Infer a participants home location

origDf = pd.read_csv(snakemake.input[0])
filteredDf = filterDatafromDf(origDf)
if filteredDf.empty:
    filteredDf.to_csv(snakemake.output[0])
else:
    dbscan_eps = snakemake.params["dbscan_eps"]
    dbscan_minsamples = snakemake.params["dbscan_minsamples"]
    threshold_static = snakemake.params["threshold_static"]
    clustering_algorithm = snakemake.params["clustering_algorithm"]

    if clustering_algorithm == "DBSCAN":
        hyperparameters = {'eps' : distance_to_degrees(dbscan_eps), 'min_samples': dbscan_minsamples}
    elif clustering_algorithm == "OPTICS":
        hyperparameters = {'max_eps': distance_to_degrees(dbscan_eps), 'min_samples': 2, 'metric':'euclidean', 'cluster_method' : 'dbscan'} 
    else:
        raise ValueError("config[PHONE_LOCATIONS][HOME_INFERENCE][CLUSTERING ALGORITHM] only accepts DBSCAN or OPTICS but you provided ",clustering_algorithm)

    filteredDf = cluster_and_label(filteredDf,clustering_algorithm,threshold_static,**hyperparameters)

    origDf['home_latitude'] = filteredDf[filteredDf['location_label']==1][['double_latitude','double_longitude']].mean()['double_latitude']
    origDf['home_longitude'] = filteredDf[filteredDf['location_label']==1][['double_latitude','double_longitude']].mean()['double_longitude']

    distanceFromHome = haversine(origDf.double_longitude,origDf.double_latitude,origDf.home_longitude,origDf.home_latitude)

    finalDf = origDf.drop(['home_latitude','home_longitude'], axis=1)
    finalDf.insert(len(finalDf.columns)-1,'distancefromhome',distanceFromHome)
    finalDf.to_csv(snakemake.output[0], index=False)


