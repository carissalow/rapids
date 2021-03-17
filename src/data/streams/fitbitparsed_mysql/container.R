# if you need a new package, you should add it with renv::install(package) so your renv venv is updated
library(RMariaDB)
library(yaml)

#' @description
#' Auxiliary function to parse the connection credentials from a specifc group in ./credentials.yaml
#' You can reause most of this function if you are connection to a DB or Web API.
#' It's OK to delete this function if you don't need credentials, e.g., you are pulling data from a CSV for example.
#' @param group the yaml key containing the credentials to connect to a database
#' @preturn dbEngine a database engine (connection) ready to perform queries
get_db_engine <- function(group){
  # The working dir is aways RAPIDS root folder, so your credentials file is always /credentials.yaml
  credentials <- read_yaml("./credentials.yaml")
  if(!group %in% names(credentials))
    stop(paste("The credentials group",group, "does not exist in ./credentials.yaml. The only groups that exist in that file are:", paste(names(credentials), collapse = ","), ". Did you forget to set the group in [PHONE_DATA_STREAMS][fitbitparsed_mysql][DATABASE_GROUP] in config.yaml?"))
  dbEngine <- dbConnect(MariaDB(), db = credentials[[group]][["database"]],
                                    username = credentials[[group]][["user"]],
                                    password = credentials[[group]][["password"]],
                                    host = credentials[[group]][["host"]],
                                    port = credentials[[group]][["port"]])
  return(dbEngine)
}

#' @description
#' Gets the sensor data for a specific device id from a database table, file or whatever source you want to query
#' 
#' @param stream_parameters The PHONE_STREAM_PARAMETERS key in config.yaml. If you need specific parameters add them there.
#' @param device A device ID string
#' @param sensor_container database table or file containing the sensor data for all participants. This is the FITBIT_SENSOR[CONTAINER] key in config.yaml
#' @param columns the columns needed from this sensor (we recommend to only return these columns instead of every column in sensor_container)
#' @return A dataframe with the sensor data for device

pull_data <- function(stream_parameters, device, sensor, sensor_container, columns){
  dbEngine <- get_db_engine(stream_parameters$DATABASE_GROUP)
  query <- paste0("SELECT ", paste(columns, collapse = ",")," FROM ", sensor_container, " WHERE ",columns$DEVICE_ID," = '", device,"'")
  # Letting the user know what we are doing
  message(paste0("Executing the following query to download data: ", query)) 
  sensor_data <- dbGetQuery(dbEngine, query)
  
  dbDisconnect(dbEngine)
  
  if(nrow(sensor_data) == 0)
    warning(paste("The device '", device,"' did not have data in ", sensor_container))

  return(sensor_data)
}

