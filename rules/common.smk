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

# Features.smk #########################################################################################################
def find_features_files(wildcards):
    feature_files = []
    for provider_key, provider in config[(wildcards.sensor_key).upper()]["PROVIDERS"].items():
        if provider["COMPUTE"]:
            feature_files.extend(expand("data/interim/{{pid}}/{sensor_key}_features/{sensor_key}_{language}_{provider_key}.csv", sensor_key=wildcards.sensor_key.lower(), language=provider["SRC_LANGUAGE"].lower(), provider_key=provider_key.lower()))
    return(feature_files)

def optional_steps_sleep_input(wildcards):
    if config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"] == True and config["STEP"]["EXCLUDE_SLEEP"]["TYPE"] == "FITBIT_BASED":
        return  "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv"
    else:
        return []

def input_merge_sensor_features_for_individual_participants(wildcards):
    feature_files = []
    for config_key in config.keys():
        if config_key.startswith(("PHONE", "FITBIT")) and "PROVIDERS" in config[config_key]:
            for provider_key, provider in config[config_key]["PROVIDERS"].items():
                if "COMPUTE" in provider.keys() and provider["COMPUTE"]:
                    feature_files.append("data/processed/features/{pid}/" + config_key.lower() + ".csv")
                    break
    return feature_files

# Reports.smk ###########################################################################################################

def optional_heatmap_days_by_sensors_input(wildcards):
    platform = infer_participant_platform("data/external/"+wildcards.pid)
    
    if platform == "android":
        tables_platform = [table for table in config["HEATMAP_DAYS_BY_SENSORS"]["DB_TABLES"] if table not in [config["CONVERSATION"]["DB_TABLE"]["IOS"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]]] # for android, discard any ios tables that may exist
    elif platform == "ios":
        tables_platform = [table for table in config["HEATMAP_DAYS_BY_SENSORS"]["DB_TABLES"] if table not in [config["CONVERSATION"]["DB_TABLE"]["ANDROID"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"]]] # for ios, discard any android tables that may exist

    return expand("data/raw/{{pid}}/{table}_with_datetime.csv", table = tables_platform)
