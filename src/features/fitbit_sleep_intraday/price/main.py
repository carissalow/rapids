import pandas as pd
import itertools



def featuresFullNames(intraday_features_to_compute, sleep_levels_to_compute, day_types_to_compute, levels_include_all_groups):

    features_fullnames = ["local_segment"]

    sleep_level_with_group = []
    for sleep_level_group in sleep_levels_to_compute:
        for sleep_level in sleep_levels_to_compute[sleep_level_group]:
            sleep_level_with_group.append(sleep_level + sleep_level_group.lower())

    for feature in intraday_features_to_compute:
        if feature == "avgduration":
            features_fullnames.extend(["avgduration" + x[0] + "main" + x[1].lower() for x in itertools.product(sleep_level_with_group + (["all"] if levels_include_all_groups else []), day_types_to_compute)])
        elif feature == "avgratioduration":
            features_fullnames.extend(["avgratioduration" + x[0] + "withinmain" + x[1].lower() for x in itertools.product(sleep_level_with_group, day_types_to_compute)])
        elif feature in ["avgstarttimeofepisodemain", "avgendtimeofepisodemain", "avgmidpointofepisodemain", "stdstarttimeofepisodemain", "stdendtimeofepisodemain", "stdmidpointofepisodemain"]:
            features_fullnames.extend([feature + x.lower() for x in day_types_to_compute])
        else:
            features_fullnames.append(feature)
    
    return features_fullnames

def mergeSleepEpisodes(sleep_data, cols_for_groupby, base_sleep_levels):

    sleep_level_with_group = []
    for sleep_level_group in base_sleep_levels:
        for sleep_level in base_sleep_levels[sleep_level_group]:
            sleep_level_with_group.append(sleep_level + sleep_level_group.lower())
    
    sleep_episodes = pd.DataFrame(columns=["local_segment", "durationinbed", "start_timestamp", "end_timestamp", "local_start_date_time", "local_end_date_time"] + ["duration" + x for x in sleep_level_with_group])

    if cols_for_groupby and (not sleep_data.empty):
        sleep_data = sleep_data.groupby(by=cols_for_groupby, sort=False)
        sleep_episodes = sleep_data[["duration"]].sum().rename(columns={"duration": "durationinbed"})

        sleep_episodes["start_timestamp"] = sleep_data["start_timestamp"].first()
        sleep_episodes["end_timestamp"] = sleep_data["end_timestamp"].last()

        sleep_episodes["local_start_date_time"] = sleep_data["local_start_date_time"].first()
        sleep_episodes["local_end_date_time"] = sleep_data["local_end_date_time"].last()

        for sleep_level in sleep_level_with_group:
            sleep_episodes["duration" + sleep_level] = sleep_data.apply(lambda group: group[group["level"] == sleep_level.replace("classic", "").replace("stages", "").replace("unified", "")]["duration"].sum())

        sleep_episodes.reset_index(inplace=True, drop=False)
        del sleep_episodes["type_episode_id"]
    
    return sleep_episodes

def extractDailyFeatures(sleep_data):
    daily_grouped = sleep_data.groupby(["local_segment", "fake_date"])
    daily_features = daily_grouped[["start_minutes"]].first().rename(columns={"start_minutes": "starttimeofepisodemain"})
    daily_features["endtimeofepisodemain"] = daily_grouped["end_minutes"].last()
    daily_features["midpointofepisodemain"] = (daily_features["starttimeofepisodemain"] + daily_features["endtimeofepisodemain"]) / 2
    daily_features["durationinbedmain"] = daily_grouped["durationinbed"].sum()

    for col in sleep_data.columns:
        if col.startswith("duration") and col != "durationinbed":
            daily_features[col + "main"] = daily_grouped[col].sum().fillna(0)
            daily_features["ratio" + col + "withinmain"] = daily_features[col + "main"] / daily_features["durationinbedmain"]
    daily_features.reset_index(inplace=True)

    # Only keep one row per fake_date
    daily_features.drop_duplicates(subset=["fake_date"], keep="first", inplace=True)

    # The day of the week with Monday=0, Sunday=6. Set Friday and Saturday as Weekend, others as Weekday.
    daily_features["is_weekend"] = pd.to_datetime(daily_features["fake_date"]).dt.dayofweek.apply(lambda x: 1 if (x == 4 or x == 5) else 0)
    
    return daily_features

def statsOfDailyFeatures(daily_features, day_type, sleep_levels, intraday_features_to_compute, sleep_intraday_features, levels_include_all_groups):
    if day_type == "WEEKEND":
        daily_features = daily_features[daily_features["is_weekend"] == 1]
    elif day_type == "WEEK":
        daily_features = daily_features[daily_features["is_weekend"] == 0]
    elif day_type == "ALL":
        pass
    else:
        raise ValueError("Please make sure the [FITBIT_SLEEP_INTRADAY][PROVIDERS][PRICE][DAY_TYPES] parameter in config.yaml file only contains the subset of [WEEKEND, WEEK, ALL].")

    if daily_features.empty:
        return sleep_intraday_features
    
    if sleep_intraday_features.empty:
        sleep_intraday_features = pd.DataFrame()
    
    # Average of time related features
    if "avgstarttimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","starttimeofepisodemain"]].groupby("local_segment")["starttimeofepisodemain"].mean().to_frame().rename(columns={"starttimeofepisodemain": "avgstarttimeofepisodemain" + day_type.lower()})], axis=1)
    if "avgendtimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","endtimeofepisodemain"]].groupby("local_segment")["endtimeofepisodemain"].mean().to_frame().rename(columns={"endtimeofepisodemain": "avgendtimeofepisodemain" + day_type.lower()})], axis=1)
    if "avgmidpointofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","midpointofepisodemain"]].groupby("local_segment")["midpointofepisodemain"].mean().to_frame().rename(columns={"midpointofepisodemain": "avgmidpointofepisodemain" + day_type.lower()})], axis=1)
    
    # Std of time related features
    if "stdstarttimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","starttimeofepisodemain"]].groupby("local_segment")["starttimeofepisodemain"].std().to_frame().rename(columns={"starttimeofepisodemain": "stdstarttimeofepisodemain" + day_type.lower()})], axis=1)
    if "stdendtimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","endtimeofepisodemain"]].groupby("local_segment")["endtimeofepisodemain"].std().to_frame().rename(columns={"endtimeofepisodemain": "stdendtimeofepisodemain" + day_type.lower()})], axis=1)
    if "stdmidpointofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","midpointofepisodemain"]].groupby("local_segment")["midpointofepisodemain"].std().to_frame().rename(columns={"midpointofepisodemain": "stdmidpointofepisodemain" + day_type.lower()})], axis=1)

    # Duration & Ratio features
    for sleep_level_group in sleep_levels:
        for sleep_level in sleep_levels[sleep_level_group]:
            if "avgduration" in intraday_features_to_compute:
                col = "duration" + sleep_level + sleep_level_group.lower() + "main"
                sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment", col]].groupby("local_segment")[col].mean().to_frame().rename(columns={col: "avg" + col + day_type.lower()})], axis=1)
            if "avgratioduration" in intraday_features_to_compute:
                col = "ratioduration" + sleep_level + sleep_level_group.lower() + "withinmain"
                sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment", col]].groupby("local_segment")[col].mean().to_frame().rename(columns={col: "avg" + col + day_type.lower()})], axis=1)
    if levels_include_all_groups and ("avgduration" in intraday_features_to_compute):
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment", "durationinbedmain"]].groupby("local_segment")["durationinbedmain"].mean().to_frame().rename(columns={"durationinbedmain": "avgdurationallmain" + day_type.lower()})], axis=1)

    return sleep_intraday_features

def socialJetLagFeature(daily_features, sleep_intraday_features):
    daily_features_weekend = daily_features[daily_features["is_weekend"] == 1]
    sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features_weekend[["local_segment","midpointofepisodemain"]].groupby("local_segment")["midpointofepisodemain"].mean().to_frame().rename(columns={"midpointofepisodemain": "helper_weekend"})], axis=1)
    
    daily_features_weekday = daily_features[daily_features["is_weekend"] == 0]
    sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features_weekday[["local_segment","midpointofepisodemain"]].groupby("local_segment")["midpointofepisodemain"].mean().to_frame().rename(columns={"midpointofepisodemain": "helper_weekday"})], axis=1)

    sleep_intraday_features["socialjetlag"] = sleep_intraday_features["helper_weekend"] - sleep_intraday_features["helper_weekday"]

    for col in ["helper_weekend", "helper_weekday"]:
        del sleep_intraday_features[col]
    
    return sleep_intraday_features
    
def RMSSDFeatures(daily_features, intraday_features_to_compute, sleep_intraday_features):
    
    date_idx = pd.DataFrame(pd.date_range(start=daily_features["fake_date"].min(), end=daily_features["fake_date"].max(), freq="D"), columns=["fake_date"])
    date_idx["fake_date"] = date_idx["fake_date"].dt.date
    daily_features = daily_features.merge(date_idx, on="fake_date", how="right")

    for col in ["starttimeofepisodemain", "endtimeofepisodemain", "midpointofepisodemain"]:
        daily_features[col + "_diff"] = daily_features[col].diff().pow(2)

    if "rmssdmeanstarttimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","starttimeofepisodemain_diff"]].groupby("local_segment")["starttimeofepisodemain_diff"].mean().pow(0.5).to_frame().rename(columns={"starttimeofepisodemain_diff": "rmssdmeanstarttimeofepisodemain"})], axis=1)
    if "rmssdmeanendtimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","endtimeofepisodemain_diff"]].groupby("local_segment")["endtimeofepisodemain_diff"].mean().pow(0.5).to_frame().rename(columns={"endtimeofepisodemain_diff": "rmssdmeanendtimeofepisodemain"})], axis=1)
    if "rmssdmeanmidpointofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","midpointofepisodemain_diff"]].groupby("local_segment")["midpointofepisodemain_diff"].mean().pow(0.5).to_frame().rename(columns={"midpointofepisodemain_diff": "rmssdmeanmidpointofepisodemain"})], axis=1)

    if "rmssdmedianstarttimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","starttimeofepisodemain_diff"]].groupby("local_segment")["starttimeofepisodemain_diff"].median().pow(0.5).to_frame().rename(columns={"starttimeofepisodemain_diff": "rmssdmedianstarttimeofepisodemain"})], axis=1)
    if "rmssdmedianendtimeofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","endtimeofepisodemain_diff"]].groupby("local_segment")["endtimeofepisodemain_diff"].median().pow(0.5).to_frame().rename(columns={"endtimeofepisodemain_diff": "rmssdmedianendtimeofepisodemain"})], axis=1)
    if "rmssdmedianmidpointofepisodemain" in intraday_features_to_compute:
        sleep_intraday_features = pd.concat([sleep_intraday_features, daily_features[["local_segment","midpointofepisodemain_diff"]].groupby("local_segment")["midpointofepisodemain_diff"].median().pow(0.5).to_frame().rename(columns={"midpointofepisodemain_diff": "rmssdmedianmidpointofepisodemain"})], axis=1)

    return sleep_intraday_features




def price_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    last_night_end = provider["LAST_NIGHT_END"]

    sleep_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_intraday_features = provider["FEATURES"]
    levels_include_all_groups = provider["SLEEP_LEVELS"]["INCLUDE_ALL_GROUPS"]
    requested_sleep_levels = provider["SLEEP_LEVELS"]
    requested_day_types = provider["DAY_TYPES"]

    # Name of the features this function can compute
    base_intraday_features = ["avgduration", "avgratioduration", "avgstarttimeofepisodemain", "avgendtimeofepisodemain", "avgmidpointofepisodemain", "stdstarttimeofepisodemain", "stdendtimeofepisodemain", "stdmidpointofepisodemain", "socialjetlag", "rmssdmeanstarttimeofepisodemain", "rmssdmeanendtimeofepisodemain", "rmssdmeanmidpointofepisodemain", "rmssdmedianstarttimeofepisodemain", "rmssdmedianendtimeofepisodemain", "rmssdmedianmidpointofepisodemain"]
    base_sleep_levels = {"CLASSIC": ["awake", "restless", "asleep"],
                        "STAGES": ["wake", "deep", "light", "rem"],
                        "UNIFIED": ["awake", "asleep"]}
    base_day_types = ["WEEKEND", "WEEK", "ALL"]

    # The subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features))
    sleep_levels_to_compute = {key: list(set(requested_sleep_levels[key]) & set(base_sleep_levels[key])) for key in requested_sleep_levels if key in base_sleep_levels}
    day_types_to_compute = list(set(requested_day_types) & set(base_day_types))

    # Full names
    features_fullnames = featuresFullNames(intraday_features_to_compute, sleep_levels_to_compute, day_types_to_compute, levels_include_all_groups)
    sleep_intraday_features = pd.DataFrame(columns=features_fullnames)

    # Filter by segemnts and chunk episodes
    sleep_intraday_data = filter_data_by_segment(sleep_intraday_data, time_segment)

    if sleep_intraday_data.empty:
        return sleep_intraday_features

    # Discard segments shorter than one day
    sleep_intraday_data[["segment_start_datetime", "segment_end_datetime"]] = sleep_intraday_data["local_segment"].str.split("#", expand=True)[1].str.split(",", expand=True).astype("datetime64[ns]")
    sleep_intraday_data["segment_length"] = (sleep_intraday_data["segment_end_datetime"] - sleep_intraday_data["segment_start_datetime"]).dt.total_seconds()
    sleep_intraday_data = sleep_intraday_data[sleep_intraday_data["segment_length"] >= 24 * 60 * 60 - 1]
    for col in ["segment_start_datetime", "segment_end_datetime", "segment_length"]:
        del sleep_intraday_data[col]
    # Select main sleep records
    sleep_intraday_data = sleep_intraday_data[sleep_intraday_data["is_main_sleep"] == 1]

    if sleep_intraday_data.empty:
        return sleep_intraday_features

    # Merge rows to get sleep episodes
    main_sleep_episodes = mergeSleepEpisodes(sleep_intraday_data, ["local_segment", "type_episode_id"], base_sleep_levels)

    # Extract number of minutes after midnight as start time; add duration to get the end time
    main_sleep_episodes["start_minutes"] = main_sleep_episodes["local_start_date_time"].apply(lambda x: x.hour * 60 + x.minute + x.second / 60)
    main_sleep_episodes["end_minutes"] = main_sleep_episodes["start_minutes"] + main_sleep_episodes["durationinbed"]
    # Extract fake date
    """ The rule used for fake date extraction
    if start_minutes >= last_night_end
        assign today
    else:
        assign yesterday
    """
    main_sleep_episodes["fake_date_delta"] = main_sleep_episodes[["start_minutes"]].apply(lambda row: 0 if row["start_minutes"] >= last_night_end else -1, axis=1)
    main_sleep_episodes["fake_date"] = (main_sleep_episodes["local_start_date_time"] + pd.to_timedelta(main_sleep_episodes["fake_date_delta"], unit="d")).dt.date

    # Update "start_minutes" column based on START_TIME
    main_sleep_episodes["start_minutes"] = main_sleep_episodes[["start_minutes", "fake_date_delta"]].apply(lambda row: row["start_minutes"] - 24 * 60 * row["fake_date_delta"], axis=1)
    main_sleep_episodes["end_minutes"] = main_sleep_episodes["start_minutes"] + main_sleep_episodes["durationinbed"]
    
    # Sort main sleep episodes based on fake_date and start_minutes
    main_sleep_episodes = main_sleep_episodes.sort_values(["fake_date", "start_minutes"])
    # Extract daily features
    daily_features = extractDailyFeatures(main_sleep_episodes)
    
    # Extract features per segment based on daily features
    for day_type in day_types_to_compute:
        sleep_intraday_features = statsOfDailyFeatures(daily_features, day_type, sleep_levels_to_compute, intraday_features_to_compute, sleep_intraday_features, levels_include_all_groups)
    if "socialjetlag" in intraday_features_to_compute:
        sleep_intraday_features = socialJetLagFeature(daily_features, sleep_intraday_features)
    sleep_intraday_features = RMSSDFeatures(daily_features, intraday_features_to_compute, sleep_intraday_features)

    sleep_intraday_features.index.name = "local_segment"
    sleep_intraday_features.reset_index(inplace=True)
    

    return sleep_intraday_features
