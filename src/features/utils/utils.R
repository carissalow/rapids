library("stringr")

rapids_log_tag <- "RAPIDS:"

filter_data_by_segment <- function(data, day_segment){
    # Filter the rows that belong to day_segment, and put the segment full name in a new column for grouping
    date_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2}"
    hour_regex = "[0-9]{2}:[0-9]{2}:[0-9]{2}"
    data <- data %>% 
        filter(grepl(paste0("\\[", day_segment, "#"), assigned_segments)) %>% 
        mutate(local_segment = str_extract(assigned_segments, paste0("\\[", day_segment, "#", date_regex, "#", hour_regex, "#", date_regex, "#", hour_regex, "\\]")),
                local_segment = str_sub(local_segment, 2, -2)) # get rid of first and last character([])
    return(data)
}

fetch_provider_features <- function(provider, provider_key, config_key, sensor_data_file, day_segments_file){
    sensor_features  <-  data.frame(local_segment = character(), stringsAsFactors = FALSE)

    sensor_data <-  read.csv(sensor_data_file, stringsAsFactors = FALSE)
    day_segments_labels <-  read.csv(day_segments_file, stringsAsFactors = FALSE)

    if(!"FEATURES" %in% names(provider))
        stop(paste0("Provider config[", config_key,"][PROVIDERS][", provider_key,"] is missing a FEATURES attribute in config.yaml"))

    if(provider[["COMPUTE"]] == TRUE){
        code_path <- paste0("src/features/", config_key,"/", provider[["SRC_FOLDER"]], "/main.R")  
        source(code_path)
        features_function <- match.fun(paste0(provider[["SRC_FOLDER"]], "_features"))
        day_segments <- day_segments_labels %>% pull(label)
        for (day_segment in day_segments){
            print(paste(rapids_log_tag,"Processing", config_key, provider_key, day_segment))

            features <- features_function(sensor_data, day_segment, provider)

            # Check all features names contain the provider key so they are unique
            features_names <- colnames(features %>% select(-local_segment))
            if(any(!grepl(paste0(".*(",str_to_lower(provider_key),").*"), features_names)))
            stop(paste("The name of all calls features of", provider_key," must contain its name in lower case but the following don't [", paste(features_names[!grepl(paste0(".*(",str_to_lower(provider_key),").*"), features_names)], collapse = ", "), "]"))

            sensor_features <- merge(sensor_features, features, all = TRUE)
        }
    } else {
    for(feature in provider[["FEATURES"]])
        sensor_features[,feature] <- NA
    }

    sensor_features <- sensor_features %>% separate(col = local_segment, 
                                    into = c("local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"),
                                    sep = "#", 
                                    remove = FALSE)
    return(sensor_features)
}