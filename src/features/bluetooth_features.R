source("renv/activate.R")
source("src/features/bluetooth/bluetooth_base.R")
library(dplyr)
library(tidyr)

bluetooth_data <- read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
day_segment <- snakemake@params[["day_segment"]]
requested_features <-  snakemake@params[["features"]]
features = data.frame(local_date = character(), stringsAsFactors = FALSE)

# Compute base bluetooth features
features <- merge(features, base_bluetooth_features(bluetooth_data, day_segment, requested_features), by="local_date", all = TRUE)

if(ncol(features) != length(requested_features) + 1)
  stop(paste0("The number of features in the output dataframe (=", ncol(features),") does not match the expected value (=", length(requested_features)," + 1). Verify your bluetooth feature extraction functions"))


write.csv(features, snakemake@output[[1]], row.names = FALSE)