source("renv/activate.R")
source("src/features/wifi/wifi_base.R")
library(dplyr)

wifi_data <- read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
day_segment <- snakemake@params[["day_segment"]]
requested_features <-  snakemake@params[["features"]]
features = data.frame(local_date = character(), stringsAsFactors = FALSE)

# Compute base wifi features
features <- merge(features, base_wifi_features(wifi_data, day_segment, requested_features), by="local_date", all = TRUE)

if(ncol(features) != length(requested_features) + 1)
  stop(paste0("The number of features in the output dataframe (=", ncol(features),") does not match the expected value (=", length(requested_features)," + 1). Verify your wifi feature extraction functions"))

write.csv(features, snakemake@output[[1]], row.names = FALSE)
