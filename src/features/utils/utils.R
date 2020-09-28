library("stringr")

rapids_log_tag <- "RAPIDS:"

filter_data_by_segment <- function(data, day_segment){
  # Filter the rows that belong to day_segment, and put the segment full name in a new column for grouping
  datetime_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}"
  timestamp_regex = "[0-9]{13}"
  data <- data %>% 
    filter(grepl(paste0("\\[", day_segment, "#"), assigned_segments)) %>% 
    mutate(local_segment = str_extract(assigned_segments, paste0("\\[", day_segment, "#", datetime_regex, ",", datetime_regex, ";", timestamp_regex, ",", timestamp_regex, "\\]"))) %>% 
    extract(local_segment, into = c("local_segment", "timestamps_segment"), paste0("\\[(", day_segment, "#", datetime_regex, ",", datetime_regex, ");(", timestamp_regex, ",", timestamp_regex, ")\\]")) %>% 
    select(-assigned_segments)
  return(data)
}

chunk_episodes <- function(sensor_episodes){
  columns_to_drop <- c("timestamp", "duration", "utc_date_time", "local_date_time", "local_date", "local_time", "local_hour", "local_minute", "segment_start", "segment_end", 'timestamp_plus_duration' )
  
  chunked_episodes <- sensor_episodes %>% separate(col = local_segment, 
                                                   into = c("local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"),
                                                   sep = "#", 
                                                   remove = FALSE) %>% 
    unite(col = "segment_start", "local_start_date", "local_start_time", sep = " ",remove = TRUE) %>% 
    unite(col = "segment_end", "local_end_date", "local_end_time", sep = " ",remove = TRUE) %>% 
    mutate(local_segment_label = NULL,
           timestamp_plus_duration = timestamp + (duration * 1000 * 60)) %>% 
    group_by(local_timezone) %>% 
    nest() %>% 
    mutate(
      data = map(data, ~.x %>% mutate(segment_start = as.numeric(lubridate::ymd_hms(segment_start, tz = local_timezone)) * 1000,
                             segment_end = as.numeric(lubridate::ymd_hms(segment_end, tz = local_timezone)) * 1000)),
      # We group by episode_id and those variables from the original episodes we want to keep once we summarise
      data = map(data, ~.x %>% group_by_at(vars(c("episode_id", setdiff(colnames(.x), columns_to_drop)  ))) %>%
                   summarize(chunked_start = max(first(timestamp), first(segment_start)), 
                             chunked_end = min(last(timestamp_plus_duration), last(segment_end)),
                             duration = (chunked_end - chunked_start) / (1000 * 60 ),
                             chunked_start = format(lubridate::as_datetime(chunked_start / 1000, tz = local_timezone),  "%Y-%m-%d %H:%M:%S"), 
                             chunked_end = format(lubridate::as_datetime(chunked_end  / 1000, tz = local_timezone),  "%Y-%m-%d %H:%M:%S")))
        ) %>%
    unnest(data)
    
  return(chunked_episodes)
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
    } else { # This is redundant, if COMPUTE is FALSE this script will be never executed
    for(feature in provider[["FEATURES"]])
        sensor_features[,feature] <- NA
    }

    sensor_features <- sensor_features %>% extract(col = local_segment, 
                                                into = c("local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"),
                                                "(.*)#(.*),(.*)", 
                                                remove = FALSE)
    return(sensor_features)
}