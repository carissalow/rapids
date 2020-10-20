source("renv/activate.R")
source("src/data/unify_utils.R")
library(yaml)

sensor_data <- read.csv(snakemake@input[["sensor_data"]], stringsAsFactors = FALSE)
participant_info <- snakemake@input[["participant_info"]]
sensor <- snakemake@params[["sensor"]]

participant <- read_yaml(participant_info)
platforms <- participant$PHONE$PLATFORMS
platform <- ifelse(platforms[1] == "multiple" | (length(platforms) > 1 & "android" %in% platforms & "ios" %in% platforms), "android", platforms[1])

sensor_data <- unify_data(sensor_data, sensor, platform)

write.csv(sensor_data, snakemake@output[[1]], row.names = FALSE)
