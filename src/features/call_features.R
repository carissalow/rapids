source("renv/activate.R")
source("src/features/call/call_base.R")
library(dplyr)

calls <-  read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
day_segments_labels <- read.csv(snakemake@input[["day_segments_labels"]])
requested_features <-  snakemake@params[["features"]]
call_type <- snakemake@params[["call_type"]]
features = data.frame(local_segment = character(), stringsAsFactors = FALSE)

day_segments <- day_segments_labels %>% pull(label)
for (day_segment in day_segments)
  features <- merge(features, base_call_features(calls, call_type, day_segment, requested_features),  all = TRUE)

if(ncol(features) != length(requested_features) + 1)
  stop(paste0("The number of features in the output dataframe (=", ncol(features),") does not match the expected value (=", length(requested_features)," + 1). Verify your Call feature extraction functions"))

features <- features %>% separate(col = local_segment, 
                                  into = c("segment", "local_start_date", "local_start_time", "local_end_date", "local_end_time"),
                                  sep = "#", 
                                  remove = FALSE)

write.csv(features, snakemake@output[[1]], row.names = FALSE)
