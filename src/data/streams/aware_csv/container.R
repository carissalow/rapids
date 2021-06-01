# if you need a new package, you should add it with renv::install(package) so your renv venv is updated
library(readr)

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
#' @param stream_parameters The PHONE_STREAM_PARAMETERS key in config.yaml. If you need specific parameters add them there.
#' @param device A device ID string
#' @return The OS the device ran, "android" or "ios"

infer_device_os <- function(stream_parameters, device){
  if(!dir.exists(stream_parameters$FOLDER))
    stop("[PHONE_DATA_STREAMS][aware_csv][FOLDER] does not exist: ", stream_parameters$FOLDER)
  data_file <- file.path(stream_parameters$FOLDER, "aware_device.csv")

  if(!file.exists(data_file))
    stop("If you want to infer the OS of a smartphone using aware_csv, you need to have a CSV file called aware_device.csv with a 'device_id' and 'brand' columns, but this file does not exist: ", data_file)

  # Letting the user know what we are doing
  message(paste0("Reading this CSV file: ", data_file))

  os <- read_delim_chunked(data_file, escape_backslash = TRUE, delim = ",", escape_double = FALSE, quote = "\"",
    callback = DataFrameCallback$new(function(x, pos) x[x[["device_id"]] == device, c("brand", "device_id")] ), progress = T, chunk_size = 50000)

  if(is.null(os)) # emtpy file
    os <- read.csv(data_file)

  if(nrow(os) > 0)
    return(os %>% mutate(os = ifelse(brand == "iPhone", "ios", "android")) %>% pull(os))
  else
    stop(paste("We cannot infer the OS of the following device id because it does not exist in the aware_device table:", device))
  return(os)
}


#' @description
#' Gets the sensor data for a specific device id from a database table, file or whatever source you want to query
#' 
#' @param stream_parameters The PHONE_STREAM_PARAMETERS key in config.yaml. If you need specific parameters add them there.
#' @param device A device ID string
#' @param sensor_container database table or file containing the sensor data for all participants. This is the PHONE_SENSOR[CONTAINER] key in config.yaml
#' @param columns the columns needed from this sensor (we recommend to only return these columns instead of every column in sensor_container)
#' @return A dataframe with the sensor data for device

pull_data <- function(stream_parameters, device, sensor, sensor_container, columns){
  if(!dir.exists(stream_parameters$FOLDER))
    stop("[PHONE_DATA_STREAMS][aware_csv][FOLDER] does not exist: ", stream_parameters$FOLDER)
  data_file <- file.path(stream_parameters$FOLDER, sensor_container)

  if(!file.exists(data_file))
    stop("The CSV file with ",sensor," data does not exist: '", data_file, "'. In config.yaml, configure [",sensor,"][CONTAINER] with the name of your CSV file including its '.csv' extension (you set the folder name in [PHONE_DATA_STREAMS][aware_csv]")

  if(!endsWith(data_file, ".csv"))
    stop("The data container should be a CSV file: ", data_file)

  # Letting the user know what we are doing
  message(paste0("Reading this CSV file: ", data_file))

  sensor_data <- read_delim_chunked(data_file, escape_backslash = TRUE, delim = ",", escape_double = FALSE, quote = "\"",
    callback = DataFrameCallback$new(function(x, pos) x[x[[columns$DEVICE_ID]] == device, unlist(columns, use.names = FALSE)] ), progress = T, chunk_size = 50000)
    # callback = DataFrameCallback$new(function(x, pos) subset(x,x[[columns$DEVICE_ID]] == device, select = unlist(columns))), progress = T, chunk_size = 50000)
  if(is.null(sensor_data)) # emtpy file
    sensor_data <- read.csv(data_file) %>% select(unlist(columns, use.names = FALSE))

  if(nrow(sensor_data) == 0)
    warning("The device '", device,"' did not have data in ", sensor_container)
  return(sensor_data)
}
