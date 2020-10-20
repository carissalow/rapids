source("renv/activate.R")
library(RMySQL)
library(dplyr)
library(readr)
library(stringr)
library(yaml)


participant_file <- snakemake@input[[1]]
source <- snakemake@params[["source"]]
sensor <- snakemake@params[["sensor"]]
table <- snakemake@params[["table"]]
sensor_file <- snakemake@output[[1]]

participant <- read_yaml(participant_file)
if(! "FITBIT" %in% names(participant)){
  stop(paste("The following participant file does not have a FITBIT section, create one manually or automatically (see the docs):", participant_file))
}
device_ids <- participant$FITBIT$DEVICE_IDS
unified_device_id <- tail(device_ids, 1)
# As opposed to phone data, we dont' filter by date here because data can still be in JSON format, we need to parse it first

if(source$TYPE == "DATABASE"){
  dbEngine <- dbConnect(MySQL(), default.file = "./.env", group = source$DATABASE_GROUP)
  query <- paste0("SELECT * FROM ", table, " WHERE ",source$DEVICE_ID_COLUMN," IN ('", paste0(device_ids, collapse = "','"), "')")
  sensor_data <- dbGetQuery(dbEngine, query)
  dbDisconnect(dbEngine)
  sensor_data <- sensor_data %>%
    rename(device_id = source$DEVICE_ID_COLUMN) %>% 
    mutate(device_id = unified_device_id) # Unify device_id
  
  if(FALSE) # For MoSHI use, we didn't split fitbit sensors into different tables
    sensor_data <- sensor_data %>% filter(fitbit_data_type == str_split(sensor, "_", simplify = TRUE)[[2]])

  # Droping duplicates on all columns except for _id or id
  sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), c("_id", "id"))))

  write_csv(sensor_data, sensor_file)

}
