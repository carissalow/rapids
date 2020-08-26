# Common.smk ##########################################################################################################

def infer_participant_platform(participant_file):
    with open(participant_file, encoding="ISO-8859-1") as external_file:
        external_file_content = external_file.readlines()
    platforms = external_file_content[1].strip().split(",")
    if platforms[0] == "multiple" or (len(platforms) > 1 and "android" in platforms and "ios" in platforms):
        platform = "android"
    else:
        platform = platforms[0]

    if platform not in ["android", "ios"]:
        raise ValueError("Platform (line 2) in a participant file should be 'android', 'ios', or 'multiple'. You typed '" + platforms + "'")

    return platform

# Preprocessing.smk ####################################################################################################

def optional_phone_sensed_bins_input(wildcards):
    platform = infer_participant_platform("data/external/"+wildcards.pid)
    
    if platform == "android":
        tables_platform = [table for table in config["PHONE_VALID_SENSED_BINS"]["DB_TABLES"] if table not in [config["CONVERSATION"]["DB_TABLE"]["IOS"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]]] # for android, discard any ios tables that may exist
    elif platform == "ios":
        tables_platform = [table for table in config["PHONE_VALID_SENSED_BINS"]["DB_TABLES"] if table not in [config["CONVERSATION"]["DB_TABLE"]["ANDROID"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"]]] # for ios, discard any android tables that may exist

    return expand("data/raw/{{pid}}/{table}_with_datetime.csv", table = tables_platform)

def find_day_segments_input_file(wildcards):
    for key, values in config.items():
        if  "DB_TABLE" in config[key] and config[key]["DB_TABLE"] == wildcards.sensor:
            if "DAY_SEGMENTS" in config[key]:
                return config[key]["DAY_SEGMENTS"]["FILE"]
            else:
                raise ValueError("{} should have a [DAY_SEGMENTS][FILE] parameter containing the path to its day segments file".format(wildcards.sensor))

def find_day_segments_input_type(wildcards):
    for key, values in config.items():
        if  "DB_TABLE" in config[key] and config[key]["DB_TABLE"] == wildcards.sensor:
            if "DAY_SEGMENTS" in config[key]:
                return config[key]["DAY_SEGMENTS"]["TYPE"]
            else:
                raise ValueError("{} should have a [DAY_SEGMENTS][TYPE] parameter containing INTERVAL, FREQUENCY, or EVENT".format(wildcards.sensor))

# Features.smk #########################################################################################################

def optional_ar_input(wildcards):
    platform = infer_participant_platform("data/external/"+wildcards.pid)

    if platform == "android": 
        return ["data/raw/{pid}/" + config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"] + "_with_datetime_unified.csv",
                "data/processed/{pid}/" + config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"] + "_deltas.csv"]
    elif platform == "ios":
        return ["data/raw/{pid}/"+config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]+"_with_datetime_unified.csv",
                "data/processed/{pid}/"+config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]+"_deltas.csv"]

def optional_conversation_input(wildcards):
    platform = infer_participant_platform("data/external/"+wildcards.pid)

    if platform == "android":
        return ["data/raw/{pid}/" + config["CONVERSATION"]["DB_TABLE"]["ANDROID"] + "_with_datetime_unified.csv"]
    elif platform == "ios":
        return ["data/raw/{pid}/" + config["CONVERSATION"]["DB_TABLE"]["IOS"] + "_with_datetime_unified.csv"]

def optional_location_barnett_input(wildcards):
    if config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED":
        return expand("data/raw/{{pid}}/{sensor}_resampled.csv", sensor=config["BARNETT_LOCATION"]["DB_TABLE"])
    else:
        return expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BARNETT_LOCATION"]["DB_TABLE"])

def optional_location_doryab_input(wildcards):
    if config["DORYAB_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED":
        return expand("data/raw/{{pid}}/{sensor}_resampled.csv", sensor=config["DORYAB_LOCATION"]["DB_TABLE"])
    else:
        return expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["DORYAB_LOCATION"]["DB_TABLE"])

def optional_steps_sleep_input(wildcards):
    if config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"] == True and config["STEP"]["EXCLUDE_SLEEP"]["TYPE"] == "FITBIT_BASED":
        return  "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv"
    else:
        return []

def optional_wifi_input(wildcards):
    if len(config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]) > 0 and len(config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]) == 0:
        return {"visible_access_points": expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"])}
    elif len(config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]) == 0 and len(config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]) > 0:
        return {"connected_access_points": expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"])}
    elif len(config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]) > 0 and len(config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]) > 0:
        return {"visible_access_points": expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]), "connected_access_points": expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"])}
    else:
        raise ValueError("If you are computing WIFI features you need to provide either VISIBLE_ACCESS_POINTS, CONNECTED_ACCESS_POINTS or both")


# Models.smk ###########################################################################################################

def input_merge_features_of_single_participant(wildcards):
    if wildcards.source == "phone_fitbit_features":
        return expand("data/processed/{pid}/{features}_{day_segment}.csv", pid=wildcards.pid, features=config["PARAMS_FOR_ANALYSIS"]["PHONE_FEATURES"] + config["PARAMS_FOR_ANALYSIS"]["FITBIT_FEATURES"], day_segment=wildcards.day_segment)
    else:
        return expand("data/processed/{pid}/{features}_{day_segment}.csv", pid=wildcards.pid, features=config["PARAMS_FOR_ANALYSIS"][wildcards.source.upper()], day_segment=wildcards.day_segment)

def optional_input_days_to_include(wildcards):
    if config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["ENABLED"]:
        # This input automatically trigers the rule days_to_analyse in mystudy.snakefile
        return ["data/interim/{pid}/days_to_analyse" + \
                    "_" + str(config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_BEFORE_SURGERY"]) + \
                    "_" + str(config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_IN_HOSPITAL"]) + \
                    "_" + str(config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_AFTER_DISCHARGE"]) + ".csv"]
    else:
        return []

def optional_input_valid_sensed_days(wildcards):
    if config["PARAMS_FOR_ANALYSIS"]["DROP_VALID_SENSED_DAYS"]["ENABLED"]:
        # This input automatically trigers the rule phone_valid_sensed_days in preprocessing.snakefile
        return ["data/interim/{pid}/phone_valid_sensed_days_{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins.csv"]
    else:
        return []

# Reports.smk ###########################################################################################################

def optional_heatmap_days_by_sensors_input(wildcards):
    platform = infer_participant_platform("data/external/"+wildcards.pid)
    
    if platform == "android":
        tables_platform = [table for table in config["HEATMAP_DAYS_BY_SENSORS"]["DB_TABLES"] if table not in [config["CONVERSATION"]["DB_TABLE"]["IOS"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]]] # for android, discard any ios tables that may exist
    elif platform == "ios":
        tables_platform = [table for table in config["HEATMAP_DAYS_BY_SENSORS"]["DB_TABLES"] if table not in [config["CONVERSATION"]["DB_TABLE"]["ANDROID"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"]]] # for ios, discard any android tables that may exist

    return expand("data/raw/{{pid}}/{table}_with_datetime.csv", table = tables_platform)
