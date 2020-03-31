# If you want to implement extra features, source(..) a new file and duplicate the line "features  <- merge(...)", then
# swap base_sms_features(...) for your own function

source("packrat/init.R")
source("src/features/sms/sms_base.R")
library(dplyr, warn.conflicts = FALSE)

sms <-  read.csv(snakemake@input[[1]])
day_segment <- snakemake@params[["day_segment"]]
metrics <-  snakemake@params[["metrics"]]
sms_type <-  snakemake@params[["sms_type"]]
features  <-  data.frame(local_date = character(), stringsAsFactors = FALSE)

# Compute base SMS features
features <- merge(features, base_sms_features(sms, sms_type, day_segment, metrics), by="local_date", all = TRUE)

if(ncol(features) != length(metrics) + 1)
  stop(paste0("The number of features in the output dataframe (=", ncol(features),") does not match the expected value (=", length(metrics)," + 1). Verify your SMS feature extraction functions"))

write.csv(features, snakemake@output[[1]], row.names = FALSE)
