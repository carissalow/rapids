import pandas as pd
import numpy as np
import scipy.stats as stats
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

day_segment = snakemake.params["segment"]
features = snakemake.params["features"]

#Read csv into a pandas dataframe
data = pd.read_csv(snakemake.input['gar_events'],parse_dates=['local_date_time'])
ar_deltas = pd.read_csv(snakemake.input['gar_deltas'],parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
columns = list("ar_" + str(day_segment) + "_" + column for column in features)

if data.empty:
    finalDataset = pd.DataFrame(columns = columns)
else:
    finalDataset = pd.DataFrame()
    ar_deltas = splitOvernightEpisodes(ar_deltas, [],['activity'])

    if day_segment != "daily":
        ar_deltas = splitMultiSegmentEpisodes(ar_deltas, day_segment, [])

    data.local_date_time = pd.to_datetime(data.local_date_time)
    resampledData = data.set_index(data.local_date_time)
    resampledData.drop(columns=['local_date_time'],inplace=True)

    if(day_segment!='daily'):
        resampledData = resampledData.loc[resampledData['local_day_segment'] == str(day_segment)]

    if resampledData.empty:
        finalDataset = pd.DataFrame(columns = columns)
    else:
        #Finding the count of samples of the day
        if("count" in features):
            finalDataset["ar_" + str(day_segment) + "_count"] = resampledData['activity_type'].resample('D').count()

        #Finding most common activity of the day
        if("mostcommonactivity" in features):
            finalDataset["ar_" + str(day_segment) + "_mostcommonactivity"] = resampledData['activity_type'].resample('D').apply(lambda x: stats.mode(x)[0] if len(stats.mode(x)[0]) != 0 else None)

        #finding different number of activities during a day
        if("countuniqueactivities" in features):
            finalDataset["ar_" + str(day_segment) + "_countuniqueactivities"] = resampledData['activity_type'].resample('D').nunique()
        
        #finding Number of times activity changed
        if("activitychangecount" in features):
            resampledData['activity_type_shift'] = resampledData['activity_type'].shift().fillna(resampledData['activity_type'].head(1))
            resampledData['different_activity'] = np.where(resampledData['activity_type']!=resampledData['activity_type_shift'],1,0)
            finalDataset["ar_" + str(day_segment) + "_activitychangecount"] = resampledData['different_activity'].resample('D').sum()


        deltas_features = {'sumstationary':['still','tilting'], 
                        'summobile':['on_foot','running','on_bicycle'],
                        'sumvehicle':['in_vehicle']}
        
        for column, activity_labels in deltas_features.items():
            if column in features:
                finalDataset["ar_" + str(day_segment) + "_"+str(column)] = (ar_deltas[ar_deltas['activity'].isin(pd.Series(activity_labels))]  
                        .groupby(['local_start_date'])['time_diff']  
                        .agg({"ar_" + str(day_segment) + "_" + str(column) :'sum'}))  
    
finalDataset.index.names = ['local_date']
finalDataset.to_csv(snakemake.output[0])