source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library("stringr")

# Load Ian Barnett's code. Taken from https://scholar.harvard.edu/ibarnett/software/gpsmobility
file.sources = list.files(c("src/features/phone_locations/barnett/library"), pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

create_empty_file <- function(requested_features){
  return(data.frame(local_segment= character(), 
                        hometime= numeric(), 
                        disttravelled= numeric(), 
                        rog= numeric(), 
                        maxdiam= numeric(), 
                        maxhomedist= numeric(), 
                        siglocsvisited= numeric(), 
                        avgflightlen= numeric(), 
                        stdflightlen= numeric(), 
                        avgflightdur= numeric(), 
                        stdflightdur= numeric(), 
                        probpause= numeric(), 
                        siglocentropy= numeric(), 
                        minsmissing= numeric(), 
                        circdnrtn= numeric(), 
                        wkenddayrtn= numeric(),
                        minutes_data_used= numeric()
                      ) %>% select(all_of(requested_features)))
}

barnett_features <- function(sensor_data_files, time_segment, params){
  
  location_data <-  read.csv(sensor_data_files[["sensor_data"]], stringsAsFactors = FALSE)
  location_features <- NULL

  location <- location_data
  accuracy_limit <- params[["ACCURACY_LIMIT"]]
  minutes_data_used <- params[["MINUTES_DATA_USED"]]

  # Compute what features were requested
  available_features <- c("hometime","disttravelled","rog","maxdiam", "maxhomedist","siglocsvisited","avgflightlen", "stdflightlen",
                          "avgflightdur","stdflightdur", "probpause","siglocentropy","minsmissing", "circdnrtn","wkenddayrtn")
  requested_features <- intersect(unlist(params["FEATURES"], use.names = F), available_features)
  requested_features <- c("local_segment", requested_features)
  if(minutes_data_used)
    requested_features <- c(requested_features, "minutes_data_used")

  # Excludes datasets with less than 24 hours of data
  if(max(location$timestamp) - min(location$timestamp) < 86400000)
    location <- head(location, 0)

  if (nrow(location) > 1){
    # Filter by segment and skipping any non-daily segment
    location <- location %>% filter_data_by_segment(time_segment)
    
    datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
    datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"
    location <- location %>% mutate(is_daily = str_detect(local_segment, paste0(time_segment, "#", datetime_start_regex, ",", datetime_end_regex))) 

    if(!all(location$is_daily)){
      message(paste("Barnett's location features cannot be computed for time segmentes that are not daily (cover 00:00:00 to 23:59:59 of every day). Skipping ", time_segment))
      location_features <- create_empty_file(requested_features)  
    } else {
      # Count how many minutes of data we use to get location features
      # Some minutes have multiple fused rows
      location_minutes_used <- location %>% 
        group_by(local_date, local_hour) %>% 
        summarise(n_minutes = n_distinct(local_minute), .groups = 'drop_last') %>% 
        group_by(local_date) %>% 
        summarise(minutes_data_used = sum(n_minutes), .groups = 'drop_last') %>% 
        select(local_date, minutes_data_used)

      # Save time segment to attach it later
      location_dates_segments <- location %>% select(local_date, local_segment) %>% distinct(local_date, .keep_all = TRUE)

      # Select only the columns that the algorithm needs
      all_timezones <- table(location %>% pull(local_timezone))
      location <- location %>% select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)
      if(nrow(location %>% filter(accuracy < accuracy_limit)) > 1){
        timezone <- names(all_timezones)[as.vector(all_timezones)==max(all_timezones)]
        outputMobility <- MobilityFeatures(location, ACCURACY_LIM = accuracy_limit, tz = timezone)
      } else {
        print(paste("Cannot compute Barnett location features because there are no rows with an accuracy value lower than ACCURACY_LIMIT", accuracy_limit))
        outputMobility <- NULL
      }

      if(is.null(outputMobility)){
        location_features <- create_empty_file(requested_features)
      } else{
        # Copy index (dates) as a column 
        features <- cbind(rownames(outputMobility$featavg), outputMobility$featavg)
        features <- as.data.frame(features)
        features[-1] <- lapply(lapply(features[-1], as.character), as.numeric)
        colnames(features)=c("local_date",tolower(colnames(outputMobility$featavg)))
        # Add the minute count column
        features <- left_join(features, location_minutes_used, by = "local_date")
        # Add the time segment column for consistency
        features <- left_join(features, location_dates_segments, by = "local_date")
        location_features <- features %>% select(all_of(requested_features))
      }
    } 
  } else {
        location_features <- create_empty_file(requested_features)  
  }

  if(ncol(location_features) != length(requested_features))
      stop(paste0("The number of features in the output dataframe (=", ncol(location_features),") does not match the expected value (=", length(requested_features),"). Verify your barnett location features"))
  return(location_features)
}