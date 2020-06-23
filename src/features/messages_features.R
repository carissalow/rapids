# If you want to implement extra features, source(..) a new file and duplicate the line "features  <- merge(...)", then
# swap base_sms_features(...) for your own function

source("renv/activate.R")
source("src/features/messages/messages_base.R")
library(dplyr, warn.conflicts = FALSE)

sms <-  read.csv(snakemake@input[[1]])
day_segment <- snakemake@params[["day_segment"]]
requested_features <-  snakemake@params[["features"]]
sms_type <-  snakemake@params[["messages_type"]]
features  <-  data.frame(local_date = character(), stringsAsFactors = FALSE)

# Compute base SMS features
features <- merge(features, base_sms_features(sms, sms_type, day_segment, requested_features), by="local_date", all = TRUE)

if(ncol(features) != length(requested_features) + 1)
  stop(paste0("The number of features in the output dataframe (=", ncol(features),") does not match the expected value (=", length(requested_features)," + 1). Verify your SMS feature extraction functions"))

write.csv(features, snakemake@output[[1]], row.names = FALSE)
