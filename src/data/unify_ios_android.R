source("renv/activate.R")
source("src/data/unify_utils.R")

sensor_data <- read.csv(snakemake@input[["sensor_data"]], stringsAsFactors = FALSE)
participant_info <- snakemake@input[["participant_info"]]
sensor <- snakemake@params[["sensor"]]
unifiable_sensors = snakemake@params[["unifiable_sensors"]]

platforms <- strsplit(readLines(participant_info, n=2)[[2]], ",")[[1]]
platform <- ifelse(platforms[1] == "multiple" | (length(platforms) > 1 & "android" %in% platforms & "ios" %in% platforms), "android", platforms[1])

sensor_data <- unify_data(sensor_data, sensor, platform, unifiable_sensors)

write.csv(sensor_data, snakemake@output[[1]], row.names = FALSE)
