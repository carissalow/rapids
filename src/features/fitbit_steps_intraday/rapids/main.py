import pandas as pd
import numpy as np

# calculate Gini coefficient, a measure of normalized variability. Adapted from: https://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
def calculateGiniCoefficient(x):
    x = np.array(x) 
    x = np.sort(x) 
    n = len(x) 
    i = np.arange(1, n+1) 
    gini = (np.sum(((2*i)-n-1)*x))/(n*(np.sum(x)))
    return gini

# extract lower and upper bounds of cadence band from feature name
def getCadenceBandBounds(feature):
    band = feature.split("totalminutes")[1].split("cadence")[0]
    if band == "0":
        bounds = [0, 1]
    elif band == "120plus":
        bounds = [120, np.inf]
    else:
        bounds = [int(x) for x in band.split("to")]
        bounds[1]+=1
    return bounds

# get device wear time based on definition from Tudor-Locke et al., 2017, which is used to compute uncensored average cadence feature
def getDeviceWearTimes(steps_data, threshold_device_nonwear_time):
    steps_data = steps_data.assign(iszerostepcount = np.where(steps_data["steps"] == 0, 1, 0), isnotzerostepcount = np.where(steps_data["steps"] != 0, 1, 0))
 
    # put consecutive rows into the same group if they have the same values of "iszerostepcount", "local_timezone", and "local_segment"
    steps_data["group_idx"] = (steps_data[["iszerostepcount", "local_timezone", "local_segment"]].shift() != steps_data[["iszerostepcount", "local_timezone", "local_segment"]]).any(axis = 1).cumsum()
    weartimes = steps_data.groupby(["local_segment", "local_timezone", "group_idx"]).agg(sumsteps = ("steps", "sum"), zerostepduration = ("iszerostepcount", "sum"), nonzerostepduration = ("isnotzerostepcount", "sum"))

    # create an indicator for groups of consecutive zero step counts with duration greater than or equal to the threshold (default: 60) if threshold was not set to 0
    if threshold_device_nonwear_time == 0:
        weartimes = weartimes.assign(isnonweartime = 0)
    else:
        weartimes = weartimes.assign(isnonweartime = np.where(weartimes["zerostepduration"] >= threshold_device_nonwear_time, 1, 0))
    
    # sum duration for rows where the isnonweartime indicator does not equal 1 to obtain device wear times
    weartimes = weartimes.assign(weartimeduration = np.where(weartimes["isnonweartime"] != 1, weartimes["zerostepduration"] + weartimes["nonzerostepduration"], 0))

    return weartimes

def statsFeatures(steps_data, features_to_compute, features_type, steps_features, *args, **kwargs):
    if features_type == "steps" or features_type == "sumsteps":
        col_name = "steps"
        reference_hour = kwargs["reference_hour"]
    elif features_type == "durationsedentarybout" or features_type == "durationactivebout":
        col_name = "duration"
    else:
        raise ValueError("features_type can only be one of ['steps', 'sumsteps', 'durationsedentarybout', 'durationactivebout'].")

    if "count" + features_type.replace("duration", "episode") in features_to_compute:
        steps_features["count" + features_type.replace("duration", "episode")] = steps_data.groupby(["local_segment"])[col_name].count()
    if "sum" + features_type in features_to_compute:
        steps_features["sum" + features_type] = steps_data.groupby(["local_segment"])[col_name].sum()
    if "max" + features_type in features_to_compute:
        steps_features["max" + features_type] = steps_data.groupby(["local_segment"])[col_name].max()
    if "min" + features_type in features_to_compute:
        steps_features["min" + features_type] = steps_data.groupby(["local_segment"])[col_name].min()
    if "avg" + features_type in features_to_compute:
        steps_features["avg" + features_type] = steps_data.groupby(["local_segment"])[col_name].mean()
    if "median" + features_type in features_to_compute:
        steps_features["median" + features_type] = steps_data.groupby(["local_segment"])[col_name].median()
    if "std" + features_type in features_to_compute:
        steps_features["std" + features_type] = steps_data.groupby(["local_segment"])[col_name].std()
    if (col_name == "steps") and ("firststeptime" in features_to_compute):
        steps_features["firststeptime"] = steps_data[steps_data["steps"].ne(0)].groupby(["local_segment"])["local_time"].first().apply(lambda x: (int(x.split(":")[0]) - reference_hour) * 60 + int(x.split(":")[1]) + (int(x.split(":")[2]) / 60))
    if (col_name == "steps") and ("laststeptime" in features_to_compute):
        steps_features["laststeptime"] = steps_data[steps_data["steps"].ne(0)].groupby(["local_segment"])["local_time"].last().apply(lambda x: (int(x.split(":")[0]) - reference_hour) * 60 + int(x.split(":")[1]) + (int(x.split(":")[2]) / 60))

    return steps_features

def activityFragmentationFeatures(steps_data, features_to_compute, steps_features, *args, **kwargs):
    if ("activetosedentarytransitionprobability" in features_to_compute):
        if not steps_data.empty:
            steps_features["activetosedentarytransitionprobability"] = steps_data.groupby(["local_segment"])["duration"].agg(lambda x: 1/x.mean())
        else:
            steps_features["activetosedentarytransitionprobability"] = np.nan
    if "sumdurationactivitylessthan5minutes" in features_to_compute:
        if not steps_data.empty:
            steps_features["sumdurationactivitylessthan5minutes"] = steps_data.groupby(["local_segment"]).apply(lambda x: x[x["duration"] < 5]["duration"].sum())
        else:
             steps_features["sumdurationactivitylessthan5minutes"] = np.nan
    if "sumdurationactivity5to105minutes" in features_to_compute:
        if not steps_data.empty:
            steps_features["sumdurationactivity5to105minutes"] = steps_data.groupby(["local_segment"]).apply(lambda x: x[x["duration"].between(5, 10, inclusive = "both")]["duration"].sum())
        else:
            steps_features["sumdurationactivity5to105minutes"] = np.nan
    if "sumdurationactivitygreaterthan10minutes" in features_to_compute:
        if not steps_data.empty:
            steps_features["sumdurationactivitygreaterthan10minutes"] = steps_data.groupby(["local_segment"]).apply(lambda x: x[x["duration"] > 10]["duration"].sum())
        else:
            steps_features["sumdurationactivitygreaterthan10minutes"] = np.nan
    if "ginicoefficient" in features_to_compute:
        if not steps_data.empty:
            steps_features["ginicoefficient"] = steps_data.groupby(["local_segment"])["duration"].agg(lambda x: calculateGiniCoefficient(x))
        else:
            steps_features["ginicoefficient"] = np.nan
            
    return steps_features

def walkingCadenceFeatures(steps_data, features_to_compute, steps_features, *args, **kwargs):
    threshold_device_nonwear_time = kwargs["threshold_device_nonwear_time"]

    peak_cadence_features = [x for x in features_to_compute if x.startswith("peak")]
    max_cadence_features = [x for x in features_to_compute if x.startswith("max")]
    cadence_band_features = [x for x in features_to_compute if (x.startswith("totalminutes") and not x.startswith("totalminutesabove"))]
    cadence_threshold_features = [x for x in features_to_compute if x.startswith("totalminutesabove")]

    if "meancadence" in features_to_compute:
        steps_features["meancadence"] = steps_data.groupby(["local_segment"])["steps"].mean()
    if "uncensoredmeancadence" in features_to_compute:
        weartimes = getDeviceWearTimes(steps_data, threshold_device_nonwear_time)
        steps_features["uncensoredmeancadence"] = weartimes.groupby(["local_segment"]).apply(lambda x: np.where(x["weartimeduration"].sum() > 0, x["sumsteps"].sum()/x["weartimeduration"].sum(), np.nan))
    if peak_cadence_features:
        for feature in peak_cadence_features:
            n = int(feature.split("peak")[1].split("minutecadence")[0])
            steps_features[feature] = steps_data.groupby(["local_segment"])["steps"].apply(lambda x: np.where(x.count() >= n, x.sort_values(ascending = False).head(n).mean(), np.nan))
    if max_cadence_features:
        for feature in max_cadence_features:
            n = int(feature.split("max")[1].split("minutecadence")[0])
            steps_features[feature] = steps_data.groupby(["local_segment"])["steps"].apply(lambda x: np.where(x.count() >= n, x.rolling(window = n, center = False).mean().max(), np.nan))
    if cadence_band_features:
        for feature in cadence_band_features:
            bounds = getCadenceBandBounds(feature)
            steps_features[feature] = steps_data.groupby(["local_segment"])["steps"].agg(lambda x: x.between(bounds[0], bounds[1], inclusive = "left").sum())
    if cadence_threshold_features:
        for feature in cadence_threshold_features:
            n = int(feature.split("totalminutesabove")[1].split("cadence")[0])
            steps_features[feature] = steps_data.groupby(["local_segment"])["steps"].agg(lambda x: (x > n).sum())

    return steps_features

def getBouts(steps_data):

    # put consecutive rows into the same group if they have the same values of "isactivebout", "local_timezone", and "local_segment"
    steps_data["group_idx"] =  (steps_data[["isactivebout", "local_timezone", "local_segment"]].shift() != steps_data[["isactivebout", "local_timezone", "local_segment"]]).any(axis=1).cumsum()
    
    # get bouts: duration column contains the number of minutes (rows) of sedentary and active activity for each episode
    grouped = steps_data.groupby("group_idx")
    bouts = grouped["local_segment"].agg(duration="count")
    bouts[["local_segment", "isactivebout"]] = grouped[["local_segment", "isactivebout"]].first()

    return bouts

def extractStepsFeaturesFromIntradayData(steps_intraday_data, reference_hour, threshold_active_bout, threshold_device_nonwear_time, intraday_features_to_compute_steps, intraday_features_to_compute_sedentarybout, intraday_features_to_compute_activebout, intraday_features_to_compute_activityfragmentation, intraday_features_to_compute_walkingcadence, steps_intraday_features):
    steps_intraday_features = pd.DataFrame()

    # statistics features of steps count
    steps_intraday_features = statsFeatures(steps_intraday_data, intraday_features_to_compute_steps, "steps", steps_intraday_features, reference_hour=reference_hour)

    # sedentary bout: less than THRESHOLD_ACTIVE_BOUT (default: 10) steps in a minute
    # active bout: greater or equal to THRESHOLD_ACTIVE_BOUT (default: 10) steps in a minute
    isactivebout = np.where(steps_intraday_data["steps"] < int(threshold_active_bout), 0, 1)
    steps_intraday_data = steps_intraday_data.assign(isactivebout = isactivebout)
    bouts = getBouts(steps_intraday_data)

    # statistics features of sedentary bout
    sedentary_bout = bouts[bouts["isactivebout"] == 0]
    steps_intraday_features = statsFeatures(sedentary_bout, intraday_features_to_compute_sedentarybout, "durationsedentarybout", steps_intraday_features)

    # statistics features of active bout
    active_bout = bouts[bouts["isactivebout"] == 1]
    steps_intraday_features = statsFeatures(active_bout, intraday_features_to_compute_activebout, "durationactivebout", steps_intraday_features)

    # activity fragmentation features
    steps_intraday_features = activityFragmentationFeatures(active_bout, intraday_features_to_compute_activityfragmentation, steps_intraday_features)

    # walking cadence features
    steps_intraday_features = walkingCadenceFeatures(steps_intraday_data, intraday_features_to_compute_walkingcadence, steps_intraday_features, threshold_device_nonwear_time=threshold_device_nonwear_time)

    steps_intraday_features.reset_index(inplace=True)

    return steps_intraday_features



def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    reference_hour = provider["REFERENCE_HOUR"]
    threshold_active_bout = provider["THRESHOLD_ACTIVE_BOUT"]
    threshold_device_nonwear_time = provider["THRESHOLD_DEVICE_NONWEAR_TIME"]
    threshold_minute_level_step_count = provider["THRESHOLD_MINUTE_LEVEL_STEP_COUNT"]
    include_zero_step_rows = provider["INCLUDE_ZERO_STEP_ROWS"]

    steps_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]

    requested_intraday_features_steps = [x + "steps" if x not in ["firststeptime", "laststeptime"] else x for x in requested_intraday_features["STEPS"]]
    requested_intraday_features_sedentarybout = [x + "sedentarybout" for x in requested_intraday_features["SEDENTARY_BOUT"]]
    requested_intraday_features_activebout = [x + "activebout" for x in requested_intraday_features["ACTIVE_BOUT"]]
    requested_intraday_features_activityfragmentation = requested_intraday_features["ACTIVITY_FRAGMENTATION"]
    requested_intraday_features_walkingcadence = [x + "cadence" for x in requested_intraday_features["WALKING_CADENCE"]]

    # name of the features this function can compute
    base_intraday_features_steps = ["sumsteps", "maxsteps", "minsteps", "avgsteps", "stdsteps", "firststeptime", "laststeptime"]
    base_intraday_features_sedentarybout = ["countepisodesedentarybout", "sumdurationsedentarybout", "maxdurationsedentarybout", "mindurationsedentarybout", "avgdurationsedentarybout", "stddurationsedentarybout"]
    base_intraday_features_activebout = ["countepisodeactivebout", "sumdurationactivebout", "maxdurationactivebout", "mindurationactivebout", "avgdurationactivebout", "stddurationactivebout"]
    base_intraday_features_activityfragmentation = ["activetosedentarytransitionprobability", "sumdurationactivitylessthan5minutes", "sumdurationactivity5to105minutes", "sumdurationactivitygreaterthan10minutes", "ginicoefficient"]
    base_intraday_features_walkingcadence = ["meancadence", "uncensoredmeancadence", "peak1minutecadence", "peak30minutecadence", "peak60minutecadence", "max5minutecadence", "max20minutecadence", "max30minutecadence", "max60minutecadence",
        "totalminutes0cadence", "totalminutes1to19cadence", "totalminutes20to39cadence", "totalminutes40to59cadence", "totalminutes60to79cadence", "totalminutes80to99cadence", "totalminutes100to119cadence", "totalminutes120pluscadence", 
        "totalminutesabove0cadence", "totalminutesabove19cadence", "totalminutesabove100cadence"]

    # the subset of requested features this function can compute
    intraday_features_to_compute_steps = list(set(requested_intraday_features_steps) & set(base_intraday_features_steps))
    intraday_features_to_compute_sedentarybout = list(set(requested_intraday_features_sedentarybout) & set(base_intraday_features_sedentarybout))
    intraday_features_to_compute_activebout = list(set(requested_intraday_features_activebout) & set(base_intraday_features_activebout))
    intraday_features_to_compute_activityfragmentation  = list(set(requested_intraday_features_activityfragmentation) & set(base_intraday_features_activityfragmentation))
    intraday_features_to_compute_walkingcadence  = list(set(requested_intraday_features_walkingcadence) & set(base_intraday_features_walkingcadence))

    intraday_features_to_compute = intraday_features_to_compute_steps + intraday_features_to_compute_sedentarybout + intraday_features_to_compute_activebout

    # exclude rows when the total step count is ZERO during the whole day
    if (not steps_intraday_data.empty) and (not include_zero_step_rows):
        dailycountstep = steps_intraday_data.groupby(["local_date"])[["steps"]].sum()
        zerocountdates = dailycountstep[dailycountstep["steps"] == 0].index.tolist()
        steps_intraday_data = steps_intraday_data[~steps_intraday_data["local_date"].isin(zerocountdates)]

    # exclude rows with minute-level step counts above the specified threshold (e.g., 1000 steps/minute)
    if (not steps_intraday_data.empty) and (threshold_minute_level_step_count > 0):
        highstepcountrows = steps_intraday_data[steps_intraday_data["steps"] > threshold_minute_level_step_count].index.tolist()
        steps_intraday_data = steps_intraday_data[~steps_intraday_data.index.isin(highstepcountrows)]

    # extract features from intraday features
    steps_intraday_features = pd.DataFrame(columns=["local_segment"] + intraday_features_to_compute)
    if not steps_intraday_data.empty:
        steps_intraday_data = filter_data_by_segment(steps_intraday_data, time_segment)

        if not steps_intraday_data.empty:
            steps_intraday_features = extractStepsFeaturesFromIntradayData(steps_intraday_data, reference_hour, threshold_active_bout, threshold_device_nonwear_time, intraday_features_to_compute_steps, intraday_features_to_compute_sedentarybout, intraday_features_to_compute_activebout, intraday_features_to_compute_activityfragmentation, intraday_features_to_compute_walkingcadence, steps_intraday_features)
    
    return steps_intraday_features
