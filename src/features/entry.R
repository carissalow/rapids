source("renv/activate.R")
source("src/features/utils/utils.R")
library("dplyr",warn.conflicts = F)
library("tidyr")

sensor_data_files <- snakemake@input

provider <- snakemake@params["provider"][["provider"]]
provider_key <- snakemake@params["provider_key"]
sensor_key <- snakemake@params["sensor_key"]

if("time_segments_labels" %in% names(sensor_data_files)){
    # Extract sensor features
    sensor_data_files$time_segments_labels <- NULL
    time_segments_file <- snakemake@input[["time_segments_labels"]]
    sensor_features <- fetch_provider_features(provider, provider_key, sensor_key, sensor_data_files, time_segments_file)
}else{
    # Data cleaning
    sensor_features = run_provider_cleaning_script(provider, provider_key, sensor_key, sensor_data_files)
}

write.csv(sensor_features, snakemake@output[[1]], row.names = FALSE)