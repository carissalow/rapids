source("renv/activate.R")
source("src/features/utils/utils.R")
library("dplyr")
library("stringr")
library("tidyr")

location_data <-  read.csv(snakemake@input[["location_data"]], stringsAsFactors = FALSE)
day_segments_labels <-  read.csv(snakemake@input[["day_segments_labels"]], stringsAsFactors = FALSE)
provider <- snakemake@params["provider"][["provider"]]
provider_key <- snakemake@params["provider_key"]

location_features  <-  data.frame(local_segment = character(), stringsAsFactors = FALSE)

if(!"FEATURES" %in% names(provider))
        stop(paste0("Provider config[LOCATION][PROVIDERS][", provider_key,"] is missing a FEATURES attribute in config.yaml"))

if(provider[["COMPUTE"]] == TRUE){
  code_path <- paste0("src/features/location/", provider[["SRC_FOLDER"]], "/main.R")  
  source(code_path)
  features_function <- match.fun(paste0(provider[["SRC_FOLDER"]], "_location_features"))
  day_segments <- day_segments_labels %>% pull(label)
  for (day_segment in day_segments){
    print(paste(rapids_log_tag,"Processing", provider_key, day_segment))

    features <- features_function(location_data, day_segment, provider)

    # Check all features names contain the provider key so they are unique
    features_names <- colnames(features %>% select(-local_segment))
    if(any(!grepl(paste0(".*(",str_to_lower(provider_key),").*"), features_names)))
      stop(paste("The name of all location features of", provider_key," must contain its name in lower case but the following don't [", paste(features_names[!grepl(paste0(".*(",str_to_lower(provider_key),").*"), features_names)], collapse = ", "), "]"))

    location_features <- merge(location_features, features, all = TRUE)
  }
} else {
  for(feature in provider[["FEATURES"]])
    location_features[,feature] <- NA
}

location_features <- location_features %>% separate(col = local_segment, 
                                  into = c("local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"),
                                  sep = "#", 
                                  remove = FALSE)

write.csv(location_features, snakemake@output[[1]], row.names = FALSE)
