# if you need a new package, you should add it with renv::install(package) so your renv venv is updated
library(influxdbr)
library(tidyverse)
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
    stop(paste("The credentials group",group, "does not exist in ./credentials.yaml. The only groups that exist in that file are:", paste(names(credentials), collapse = ",")))
  
  #replace with credentials values
  conn_object <- influx_connection(host=credentials[[group]][["host"]], 
                                   user=credentials[[group]][["user"]],
                                   pass=credentials[[group]][["password"]],
                                   port= credentials[[group]][["port"]])
  
  return(conn_object)
}

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
  dbEngine <- get_db_engine(stream_parameters$DATABASE_GROUP)
  #need to re-fetch the YAML for the DB name
  credentials <- read_yaml("./credentials.yaml")
  message(paste0("Utilizing the Influx query for: ", device)) 
  #execute query string
  query_object <- influx_select(dbEngine,
                                db = credentials[[stream_parameters$DATABASE_GROUP]][["database"]],
                                field_keys="device_id,brand",
                                measurement="aware_device",
                                where= paste0("device_id = '",device,"'"),
                                return_xts = FALSE)
  
  #fetches the table from the query_object, filtering rows with ALL n/a
  #a behavior of influxdbr is that one all NA row will be returned with no matches
  columns = c("brand", "device_id")
  if(! all(columns %in% colnames( query_object[[1]])))
    os <- data.frame(matrix(ncol=length(columns),nrow=0, dimnames=list(NULL, columns)))
  else
    os <- query_object[[1]] %>% filter_all(any_vars(!is.na(.))) %>% select(columns)
  
  
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
  dbEngine <- get_db_engine(stream_parameters$DATABASE_GROUP)
  #need to re-fetch the YAML for the DB name
  credentials <- read_yaml("./credentials.yaml")
  

  # Letting the user know what we are doing
  message(paste0("Executing an Influx query for: ", device, " ", sensor, ". Extracting ", columns, " from ", sensor_container)) 
  #execute query string
  query_object <- influx_select(dbEngine,
                                db = credentials[[stream_parameters$DATABASE_GROUP]][["database"]],
                                field_keys=paste(columns, collapse = ","),
                                measurement=sensor_container,
                                where= paste0(columns$DEVICE_ID, " = '",device,"'"),
                                return_xts=FALSE)  

  columns = unlist(columns, use.names = FALSE)
  if(! all(columns %in% colnames( query_object[[1]])))
    sensor_data <- data.frame(matrix(ncol=length(columns),nrow=0, dimnames=list(NULL, columns)))
  else
    sensor_data <- query_object[[1]] %>% filter_all(any_vars(!is.na(.))) %>% select(columns)
  
  if(nrow(sensor_data) == 0)
    warning(paste("The device '", device,"' did not have data in ", sensor_container))

  return(sensor_data)
}

