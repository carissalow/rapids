import pandas as pd
import numpy as np
import scipy.stats as stats
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

day_segment = snakemake.params["segment"]

#Read csv into a pandas dataframe
data = pd.read_csv(snakemake.input['gar_events'],parse_dates=['local_date_time'])
ar_deltas = pd.read_csv(snakemake.input['gar_deltas'],parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
columns = ['count','most_common_activity','count_unique_activities','activity_change_count','sumstationary','summobile','sumvehicle']
columns = list("ar_" + str(day_segment) + "_" + column for column in columns)


if data.empty:
    finalDataset = pd.DataFrame(columns = columns)
else:

    ar_deltas = splitOvernightEpisodes(ar_deltas, [],['activity'])

    if day_segment != "daily":
        ar_deltas = splitMultiSegmentEpisodes(ar_deltas, day_segment, [])

    data.local_date_time = pd.to_datetime(data.local_date_time)
    resampledData = data.set_index(data.local_date_time)
    resampledData.drop(columns=['local_date_time'],inplace=True)

    if(day_segment!='daily'):
        resampledData = resampledData.loc[resampledData['local_day_segment'] == str(day_segment)]
    
    count = resampledData['activity_type'].resample('D').count()

    #Finding most common activity of the day
    mostCommonActivity = resampledData['activity_type'].resample('D').apply(lambda x:stats.mode(x)[0])

    #finding different number of activities during a day
    uniqueActivities = resampledData['activity_type'].resample('D').nunique()
    
    #finding Number of times activity changed
    resampledData['activity_type_shift'] = resampledData['activity_type'].shift().fillna(resampledData['activity_type'].head(1),inplace=True)
    resampledData['different_activity'] = np.where(resampledData['activity_type']!=resampledData['activity_type_shift'],1,0)
    countChanges = resampledData['different_activity'].resample('D').sum()
    finalDataset = pd.concat([count, mostCommonActivity, uniqueActivities, countChanges],axis=1)

    deltas_metrics = {'sumstationary':['still','tilting'], 
                    'summobile':['on_foot','running','on_bicycle'],
                    'sumvehicle':['on_vehicle']}
  
    for column, activity_labels in deltas_metrics.items():
        metric = (ar_deltas[ar_deltas['activity'].isin(pd.Series(activity_labels))]  
                .groupby(['local_start_date'])['time_diff']  
                .agg({"ar_" + str(day_segment) + "_" + str(column) :'sum'}))  
        finalDataset = finalDataset.merge(metric,how='outer',left_index=True,right_index=True)
    
finalDataset.fillna(0,inplace=True)
finalDataset.index.names = ['local_date']
finalDataset.columns=columns
finalDataset.to_csv(snakemake.output[0])