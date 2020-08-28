source("renv/activate.R")
source("src/features/utils/utils.R")
library("dplyr")
library("tidyr")

sensor_data_file <-  snakemake@input[["sensor_data"]]
day_segments_file <-  snakemake@input[["day_segments_labels"]]
provider <- snakemake@params["provider"][["provider"]]
provider_key <- snakemake@params["provider_key"]

sensor_features <- fetch_provider_features(provider, provider_key, "calls", sensor_data_file, day_segments_file)

write.csv(sensor_features, snakemake@output[[1]], row.names = FALSE)
