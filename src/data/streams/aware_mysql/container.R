# if you need a new package, you should add it with renv::install(package) so your renv venv is updated
library(RMariaDB)

# This file gets executed for each PHONE_SENSOR of each participant
# If you are connecting to a database the env file containing its credentials is available at "./.env"
# If you are reading a CSV file instead of a DB table, the @param sensor_container wil contain the file path as set in config.yaml
# You are not bound to databases or files, you can query a web API or whatever data source you need.

#' @description
#' RAPIDS allows users to use the keyword "infer" (previously "multiple") to automatically infer the mobile Operative System a device was running.
#' If you have a way to infer the OS of a device ID, implement this function. For example, for AWARE data we use the "aware_device" table.
#'  
#' If you don't have a way to infer the OS, call stop("Error Message") so other users know they can't use "infer" or the inference failed, 
#' and they have to assign the OS manually in the participant file
#' 
#' @param data_configuration The PHONE_DATA_CONFIGURATION key in config.yaml. If you need specific parameters add them there.
#' @param device A device ID string
#' @return The OS the device ran, "android" or "ios"

infer_device_os <- function(data_configuration, device){
  group <- data_configuration$SOURCE$DATABASE_GROUP # specified DB credentials group in config.yaml
  
  dbEngine <- dbConnect(MariaDB(), default.file = "./.env", group = group)
  query <- paste0("SELECT device_id,brand FROM aware_device WHERE device_id = '", device, "'")
  message(paste0("Executing the following query to infer phone OS: ", query)) 
  os <- dbGetQuery(dbEngine, query)
  dbDisconnect(dbEngine)
  
  if(nrow(os) > 0)
    return(os %>% mutate(os = ifelse(brand == "iPhone", "ios", "android")) %>% pull(os))
  else
    stop(paste("We cannot infer the OS of the following device id because it does not exist in the aware_device table:", device))
  
  return(os)
}

#' @description
#' Gets the sensor data for a specific device id from a database table, file or whatever source you want to query
#' 
#' @param data_configuration The PHONE_DATA_CONFIGURATION key in config.yaml. If you need specific parameters add them there.
#' @param device A device ID string
#' @param sensor_container database table or file containing the sensor data for all participants. This is the PHONE_SENSOR[TABLE] key in config.yaml
#' @param columns the columns needed from this sensor (we recommend to only return these columns instead of every column in sensor_container)
#' @return A dataframe with the sensor data for device

download_data <- function(data_configuration, device, sensor_container, columns){
  group <- data_configuration$SOURCE$DATABASE_GROUP
  dbEngine <- dbConnect(MariaDB(), default.file = "./.env", group = group)
  
  
  query <- paste0("SELECT ", paste(columns, collapse = ",")," FROM ", sensor_container, " WHERE device_id = '", device,"'")
  # Letting the user know what we are doing
  message(paste0("Executing the following query to download data: ", query)) 
  sensor_data <- dbGetQuery(dbEngine, query)
  
  dbDisconnect(dbEngine)
  
  if(nrow(sensor_data) == 0)
    warning(paste("The device '", device,"' did not have data in ", sensor_container))

  return(sensor_data)
}

