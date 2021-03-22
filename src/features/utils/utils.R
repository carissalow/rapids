library("stringr")

rapids_log_tag <- "RAPIDS:"

filter_data_by_segment <- function(data, time_segment){
  # Filter the rows that belong to time_segment, and put the segment full name in a new column for grouping
  datetime_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}"
  timestamp_regex = "[0-9]{13}"
  data <- data %>% 
    filter(grepl(paste0("\\[", time_segment, "#"), assigned_segments)) %>% 
    mutate(local_segment = str_extract(assigned_segments, paste0("\\[", time_segment, "#", datetime_regex, ",", datetime_regex, ";", timestamp_regex, ",", timestamp_regex, "\\]"))) %>% 
    extract(local_segment, into = c("local_segment", "timestamps_segment"), paste0("\\[(", time_segment, "#", datetime_regex, ",", datetime_regex, ");(", timestamp_regex, ",", timestamp_regex, ")\\]")) %>% 
    select(-assigned_segments)
  return(data)
}

chunk_episodes <- function(sensor_episodes){
  columns_to_drop <- c("^timestamp$","local_date_time", "local_date", "local_time", "local_hour", "local_minute", "segment_start", "segment_end" )

  chunked_episodes <- sensor_episodes %>% 
    separate(col = timestamps_segment,
             into = c("segment_start_timestamp", "segment_end_timestamp"),
             sep = ",", convert = TRUE, remove = TRUE) %>% 
    group_by(local_timezone) %>% 
    nest() %>% 
    mutate(data = map(data, ~.x %>% 
                   distinct(start_timestamp, end_timestamp, local_segment, .keep_all = TRUE) %>% 
                   mutate(start_timestamp = pmax(start_timestamp, segment_start_timestamp),
                          end_timestamp = pmin(end_timestamp, segment_end_timestamp),
                          duration = (end_timestamp - start_timestamp) / (1000 * 60)) %>% 
                   select(-matches(columns_to_drop)) %>% 
                   group_by_at(vars(setdiff(colnames(.), c("start_timestamp", "end_timestamp", "duration")))) %>%
                   summarize(start_timestamp = first(start_timestamp),
                             end_timestamp = last(end_timestamp),
                             duration = sum(duration)) %>% 
                   mutate(local_start_date_time = format(lubridate::as_datetime(start_timestamp / 1000, tz = local_timezone),  "%Y-%m-%d %H:%M:%S"),
                          local_end_date_time = format(lubridate::as_datetime(end_timestamp  / 1000, tz = local_timezone),  "%Y-%m-%d %H:%M:%S")) %>% 
                   ungroup())
           ) %>%
    unnest(data) %>% 
    ungroup() %>% 
    select(-local_timezone)
  
  return(chunked_episodes)
}

fetch_provider_features <- function(provider, provider_key, sensor_key, sensor_data_files, time_segments_file){
    sensor_features  <-  data.frame(local_segment = character(), stringsAsFactors = FALSE)

    time_segments_labels <-  read.csv(time_segments_file, stringsAsFactors = FALSE)

    if(!"FEATURES" %in% names(provider))
        stop(paste0("Provider config[", sensor_key,"][PROVIDERS][", provider_key,"] is missing a FEATURES attribute in config.yaml"))

    if(provider[["COMPUTE"]] == TRUE){
        source(provider[["SRC_SCRIPT"]])
        features_function <- match.fun(paste0(tolower(provider_key), "_features"))
        time_segments <- time_segments_labels %>% pull(label)
        for (time_segment in time_segments){
            print(paste(rapids_log_tag,"Processing", sensor_key, provider_key, time_segment))

            features <- features_function(sensor_data_files, time_segment, provider)
            if(!"local_segment" %in% colnames(features))
              stop(paste0("The dataframe returned by the ",sensor_key," provider '", provider_key,"' is missing the 'local_segment' column added by the 'filter_data_by_segment()' function. Check the provider script is using such function and is not removing 'local_segment' by accident (", provider[["SRC_SCRIPT"]],")\n  The 'local_segment' column is used to index a provider's features (each row corresponds to a different time segment instance (e.g. 2020-01-01, 2020-01-02, 2020-01-03, etc.)"))
            features <- features %>% rename_at(vars(!matches("local_segment")), ~ paste(sensor_key, provider_key, ., sep = "_"))
            sensor_features <- merge(sensor_features, features, all = TRUE)
        }
    } else { # This is redundant, if COMPUTE is FALSE this script will be never executed
    for(feature in provider[["FEATURES"]])
        sensor_features[,feature] <- NA
    }
    sensor_features <- sensor_features %>% mutate(local_segment = str_remove(local_segment, "_RR\\d+SS")) %>% 
                                            extract(col = local_segment, 
                                                into = c("local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"),
                                                "(.*)#(.*),(.*)", 
                                                remove = FALSE)
    return(sensor_features)
}