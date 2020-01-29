import pandas as pd
import numpy as np
import datetime as dt
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

day_segment = snakemake.params["day_segment"]
all_steps = snakemake.params["metrics_all_steps"]
sedentary_bout = snakemake.params["metrics_sedentary_bout"]
active_bout = snakemake.params["metrics_active_bout"]
threshold_active_bout = snakemake.params['threshold_active_bout']

#Read csv into a pandas dataframe
data = pd.read_csv(snakemake.input['steps_data'],parse_dates=['local_date_time'])
columns = list("step_" + str(day_segment) + "_" + column for column in (all_steps + sedentary_bout + active_bout))

if (day_segment != 'daily'):
    data = data.loc[data['local_day_segment'] == str(day_segment)]
    
if data.empty:
    finalDataset = pd.DataFrame(columns = columns)
else:
    finalDataset = pd.DataFrame()

    #Preprocessing:
    data.local_date_time = pd.to_datetime(data.local_date_time)
    resampledData = data.set_index(data.local_date_time)
    resampledData.index.names = ['datetime']

    resampledData['time_diff_minutes'] = resampledData['local_date_time'].diff().fillna(pd.Timedelta(seconds=0)).dt.total_seconds().div(60).astype(int)

    #Sedentary Bout when you have less than 10 steps in a minute
    #Active Bout when you have greater or equal to 10 steps in a minute
    resampledData['active_sedentary'] = np.where(resampledData['steps']<int(threshold_active_bout),'sedentary','active')
    
    activeData = resampledData[resampledData['active_sedentary']=="active"]
    sedentaryData = resampledData[resampledData['active_sedentary']=="sedentary"]

    #Time Calculations of sedentary/active bouts:
    resampledData['active_sedentary_groups'] = (resampledData.active_sedentary != resampledData.active_sedentary.shift()).cumsum().values

    #Get the total minutes for each episode
    minutesGroupedBy = resampledData.groupby(['local_date','active_sedentary','active_sedentary_groups'])['time_diff_minutes'].sum()
    
    #Get Stats for all episodes in terms of minutes
    statsMinutes = minutesGroupedBy.groupby(['local_date','active_sedentary']).agg([max,min,np.mean,np.std])
    mux = pd.MultiIndex.from_product([statsMinutes.index.levels[0], statsMinutes.index.levels[1]],names=['local_date','active_sedentary'])
    statsMinutes = statsMinutes.reindex(mux, fill_value=None).reset_index()
    statsMinutes.set_index('local_date',inplace = True)
    
    #Descriptive Statistics Features:
    if("sumallsteps" in all_steps):
        finalDataset["step_" + str(day_segment) + "_sumallsteps"] = resampledData['steps'].resample('D').sum()
    
    if("maxallsteps" in all_steps):
        finalDataset["step_" + str(day_segment) + "_maxallsteps"] = resampledData['steps'].resample('D').max()

    if("minallsteps" in all_steps):
        finalDataset["step_" + str(day_segment) + "_minallsteps"] = resampledData['steps'].resample('D').min()
    
    if("avgallsteps" in all_steps):
        finalDataset["step_" + str(day_segment) + "_avgallsteps"] = resampledData['steps'].resample('D').mean()
    
    if("stdallsteps" in all_steps):
        finalDataset["step_" + str(day_segment) + "_stdallsteps"] = resampledData['steps'].resample('D').std()
    
    if("countsedentarybout" in sedentary_bout):
        finalDataset["step_" + str(day_segment) + "_countsedentarybout"] = sedentaryData['active_sedentary'].resample('D').count()

    if("countactivebout" in active_bout):
        finalDataset["step_" + str(day_segment) + "_countactivebout"] = activeData['active_sedentary'].resample('D').count()

    if("maxdurationsedentarybout" in sedentary_bout):
        finalDataset["step_" + str(day_segment) + "_maxdurationsedentarybout"] = statsMinutes[statsMinutes['active_sedentary']=='sedentary']['max']
    
    if("mindurationsedentarybout" in sedentary_bout):
        finalDataset["step_" + str(day_segment) + "_mindurationsedentarybout"] = statsMinutes[statsMinutes['active_sedentary']=='sedentary']['min']

    if("avgdurationsedentarybout" in sedentary_bout):
        finalDataset["step_" + str(day_segment) + "_avgdurationsedentarybout"] = statsMinutes[statsMinutes['active_sedentary']=='sedentary']['mean']

    if("stddurationsedentarybout" in sedentary_bout):
        finalDataset["step_" + str(day_segment) + "_stddurationsedentarybout"] = statsMinutes[statsMinutes['active_sedentary']=='sedentary']['std']
    
    if("maxdurationactivebout" in active_bout):
        finalDataset["step_" + str(day_segment) + "_maxdurationactivebout"] = statsMinutes[statsMinutes['active_sedentary']== 'active']['max']
    
    if("mindurationactivebout" in active_bout):
        finalDataset["step_" + str(day_segment) + "_mindurationactivebout"] = statsMinutes[statsMinutes['active_sedentary']== 'active']['min']
    
    if("avgdurationactivebout" in active_bout):
        finalDataset["step_" + str(day_segment) + "_avgdurationactivebout"] = statsMinutes[statsMinutes['active_sedentary']== 'active']['mean']
    
    if("stddurationactivebout" in active_bout):
        finalDataset["step_" + str(day_segment) + "_stddurationactivebout"] = statsMinutes[statsMinutes['active_sedentary']== 'active']['std']

finalDataset.index.names = ['local_date']
finalDataset.to_csv(snakemake.output[0])    
