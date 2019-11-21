import pandas as pd
import numpy as np
import scipy.stats as stats

day_segment = snakemake.params["segment"]

#Read csv into a pandas dataframe
data = pd.read_csv(snakemake.input[0])
columns = ['count','most_common_activity','count_unique_activities','activity_change_count']
columns = list("ar_" + str(day_segment) + "_" + column for column in columns)

if data.empty:
    finalDataset = pd.DataFrame(columns = columns)
else:
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

finalDataset.index.names = ['local_date']
finalDataset.columns=columns
finalDataset.to_csv(snakemake.output[0])