source("renv/activate.R")

library("tidyr")
library("dplyr", warn.conflicts = F)

location_features_files <- snakemake@input[["sensor_features"]]
location_features <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("local_segment"))


for(location_features_file in location_features_files){
    location_features <- merge(location_features, read.csv(location_features_file), all = TRUE)
}

write.csv(location_features %>% arrange(local_segment), snakemake@output[[1]], row.names = FALSE)