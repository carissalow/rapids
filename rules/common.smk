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
        if config_key.startswith(("PHONE", "FITBIT", "EMPATICA")) and "PROVIDERS" in config[config_key]:
            for provider_key, provider in config[config_key]["PROVIDERS"].items():
                if "COMPUTE" in provider.keys() and provider["COMPUTE"]:
                    feature_files.append("data/processed/features/{pid}/" + config_key.lower() + ".csv")
                    break
    return feature_files

