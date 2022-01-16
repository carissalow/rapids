source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library("stringr")
library("lubridate")
library("purrr")

# Load Ian Barnett's code. From https://scholar.harvard.edu/ibarnett/software/gpsmobility
file.sources = list.files(c("src/features/phone_locations/barnett/library"), pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
output_apply <- sapply(file.sources,source,.GlobalEnv)

create_empty_file <- function(){
  return(data.frame(local_date= character(), hometime= numeric(), disttravelled= numeric(), rog= numeric(), maxdiam= numeric(), 
                    maxhomedist= numeric(), siglocsvisited= numeric(), avgflightlen= numeric(), stdflightlen= numeric(), 
                    avgflightdur= numeric(), stdflightdur= numeric(), probpause= numeric(), siglocentropy= numeric(), minsmissing= numeric(), 
                    circdnrtn= numeric(), wkenddayrtn= numeric(), minutes_data_used= numeric()
  ))
}

barnett_daily_features <- function(snakemake){
  location_features <- NULL
  location <- read.csv(snakemake@input[["sensor_data"]], stringsAsFactors = FALSE)
  segment_labels <- read.csv(snakemake@input[["time_segments_labels"]], stringsAsFactors = FALSE)
  accuracy_limit = 999999999 # We filter rows based on accuracy in src/data/process_location_types.R script
  datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
  datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"
  location <- location %>% 
    mutate(is_daily = str_detect(assigned_segments, paste0(".*#", datetime_start_regex, ",", datetime_end_regex, ".*")))

  does_not_span = nrow(segment_labels) == 0 || nrow(location) == 0 || all(location$is_daily == FALSE) || (max(location$timestamp) - min(location$timestamp) < 86400000)

  if(is.na(does_not_span) || does_not_span){
      warning("Barnett's location features cannot be computed for data or time segments that do not span one or more entire days (00:00:00 to 23:59:59). Values below point to the problem:",
            "\nLocation data rows within a daily time segment: ", nrow(filter(location, is_daily)),
            "\nLocation data time span in days: ", round((max(location$timestamp) - min(location$timestamp)) / 86400000, 2)
            )
    location_features <- create_empty_file()  
  } else{
    # Count how many minutes of data we use to get location features. Some minutes have multiple fused rows
    location_minutes_used <- location %>% 
      group_by(local_date, local_hour) %>% 
      summarise(n_minutes = n_distinct(local_minute), .groups = 'drop_last') %>% 
      group_by(local_date) %>% 
      summarise(minutes_data_used = sum(n_minutes), .groups = 'drop_last') %>% 
      select(local_date, minutes_data_used)
    
    # Select only the columns that the algorithm needs
    all_timezones <- table(location %>% pull(local_timezone))
    location <- location %>% select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)
    timezone <- names(all_timezones)[as.vector(all_timezones)==max(all_timezones)]
    outputMobility <- MobilityFeatures(location, ACCURACY_LIM = accuracy_limit, tz = timezone)
    
    if(is.null(outputMobility)){
      location_features <- create_empty_file()
    } else {
      # Copy index (dates) as a column 
      features <- cbind(rownames(outputMobility$featavg), outputMobility$featavg)
      features <- as.data.frame(features)
      features[-1] <- lapply(lapply(features[-1], as.character), as.numeric)
      colnames(features)=c("local_date",tolower(colnames(outputMobility$featavg)))
      location_features <- left_join(features, location_minutes_used, by = "local_date")
    }
    
  }
  write.csv(location_features, snakemake@output[[1]], row.names =FALSE)
}

barnett_daily_features(snakemake)