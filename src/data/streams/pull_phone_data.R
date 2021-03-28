source("renv/activate.R")

library(yaml)
library(dplyr)
library(readr)
# we use reticulate but only load it if we are going to use it to minimize the case when old RAPIDS deployments need to update ther renv

validate_deviceid_platforms <- function(device_ids, platforms, participant){
  if(length(device_ids) > 1 && length(platforms) == 1){
    if(platforms != "android" && platforms != "ios" && platforms != "infer")
      stop(paste0("If you have more than 1 device_id, platform should be 'android', 'ios' OR 'infer' but you typed: '", paste0(platforms, collapse = "s,"), "'. Participant file: ", participant))
  } else if(length(device_ids) > 1 && length(platforms) > 1){
    if(length(device_ids) != length(platforms))
      stop(paste0("The number of device_ids should match the number of platforms. Participant file:", participant))
    if(all(intersect(c("android", "ios"), unique(platforms)) != c("android", "ios")))
      stop(paste0("If you have more than 1 device_id and more than 1 platform, the platforms should be a mix of 'android' AND 'ios' but you typed: '", paste0(platforms, collapse = ","), "'. Participant file: ", participant))
  }
}

validate_inferred_os <- function(stream_container, participant_file, device, device_os){
  if(!is.na(device_os) && device_os != "android" && device_os != "ios")
    stop(paste0("We tried to infer the OS for ", device, " but 'infer_device_os' function inside '",stream_container,"' returned '",device_os,"' instead of 'android' or 'ios'. You can assign the OS manually in the participant file or report this bug on GitHub.\nParticipant file ", participant_file))
}

mutate_data <- function(scripts, data, data_configuration){
  for(script in scripts){
    if(grepl("\\.(R)$", script)){
      myEnv <- new.env()    
      source(script, local=myEnv)
      attach(myEnv, name="sourced_scripts_rapids")
      if(exists("main", myEnv)){
        message(paste("Applying mutation script", script))
        data <- main(data, data_configuration)
      } else{
        stop(paste0("The following mutation script does not have main function: ", script))
      }
      # rm(list = ls(envir = myEnv), envir = myEnv, inherits = FALSE)
      detach("sourced_scripts_rapids")
    } else{ # python
      library(reticulate)
      module <- gsub(pattern = "\\.py$", "", basename(script))
      script_functions <- import_from_path(module, path = dirname(script))
      if(py_has_attr(script_functions, "main")){
        message(paste("Applying mutation script", script))
        data <- script_functions$main(data, data_configuration)
      } else{
        stop(paste0("The following mutation script does not have a main function: ", script))
      }
    }
  }

  return(data)
}

rename_columns <- function(name_maps, data){
  for(name in names(name_maps))
    data <- data %>% rename(!!tolower(name) := name_maps[[name]])
  return(data)
}

validate_expected_columns_mapping <- function(schema, rapids_schema, sensor, rapids_schema_file, stream_format){
  rapids_columns <- rapids_schema[[sensor]]
  if(is.null(rapids_columns))
    stop(paste(sensor, " columns are not listed in RAPIDS' column specification. If you are adding support for a new phone sensor, add any mandatory columns in ", rapids_schema_file))
  
  if("ANDROID" %in% schema[[sensor]]){
    android_columns <- names(schema[[sensor]][["ANDROID"]][["RAPIDS_COLUMN_MAPPINGS"]])
    if(length(setdiff(rapids_columns, android_columns)) > 0)
      stop(paste(sensor," mappings are missing one or more mandatory columns for ANDROID. The missing column mappings are for ", paste(setdiff(rapids_columns, android_columns), collapse=","),"in", stream_format, " (the mappings are case sensitive)"))
    if(length(setdiff(android_columns, rapids_columns)) > 0)
      stop(paste(sensor," mappings have one or more columns than required for ANDROID. The extra column mappings are for ", paste(setdiff(android_columns, rapids_columns), collapse=","),"in", stream_format, " (the mappings are case sensitive)"))
  }

  if("IOS" %in% schema[[sensor]]){
    ios_columns <- names(schema[[sensor]][["IOS"]][["RAPIDS_COLUMN_MAPPINGS"]])
    if(length(setdiff(rapids_columns, ios_columns)) > 0)
      stop(paste(sensor," mappings are missing one or more mandatory columns for IOS. The missing column mappings are for ", paste(setdiff(rapids_columns, ios_columns), collapse=","),"in", stream_format, " (the mappings are case sensitive)"))
    if(length(setdiff(ios_columns, rapids_columns)) > 0)
      stop(paste(sensor," mappings have one or more columns than required for IOS. The extra column mappings are for ", paste(setdiff(ios_columns, rapids_columns), collapse=","),"in", stream_format, " (the mappings are case sensitive)"))
  }
}

load_container_script <- function(stream_container){
  language <- if_else(endsWith(tolower(stream_container), "py"), "python", "r")
  if(language == "python"){
    library(reticulate)
    container <- import_from_path(gsub(pattern = "\\.py$", "", basename(stream_container)), path = dirname(stream_container))
    if(!py_has_attr(container, "pull_data"))
      stop(paste0("The following container.py script does not have a pull_data function: ", stream_container))
    if(!py_has_attr(container, "infer_device_os"))
      stop(paste0("The following container.py script does not have a infer_device_os function: ", stream_container))
    return(list("infer_device_os" = container$infer_device_os, "pull_data" = container$pull_data))
  } else if(language == "r"){
    source(stream_container)
    if(!exists("pull_data"))
      stop(paste0("The following container.R script does not have a pull_data function: ", stream_container))
    if(!exists("infer_device_os"))
      stop(paste0("The following container.R script does not have a infer_device_os function: ", stream_container))
    return(list("infer_device_os" = infer_device_os, "pull_data" = pull_data))
  }
}

get_devices_ids <- function(participant_data){
  devices_ids = c()
  for(device in participant_data)
    for(attribute in names(device))
      if(attribute == "DEVICE_IDS")
        devices_ids <- c(devices_ids, device[[attribute]])
  return(devices_ids)
}

validate_participant_file_without_device_ids <- function(participant_file){
  participant_data <- yaml::read_yaml(participant_file)
  participant_devices <- get_devices_ids(participant_data)
  if(length(participant_devices) == 0)
    stop("There are no device ids in this participant file for smartphones or wearables: ", participant_file)
}

pull_phone_data <- function(){
  participant_file <- snakemake@input[["participant_file"]]
  stream_format <- snakemake@input[["stream_format"]]
  rapids_schema_file <- snakemake@input[["rapids_schema_file"]]
  stream_container <- snakemake@input[["stream_container"]]
  data_configuration <- snakemake@params[["data_configuration"]]
  tables <- snakemake@params[["tables"]]
  sensor <- toupper(snakemake@params[["sensor"]])
  device_type <- "phone"
  output_data_file <- snakemake@output[[1]]

  validate_participant_file_without_device_ids(participant_file)
  participant_data <- read_yaml(participant_file)
  stream_schema <- read_yaml(stream_format)
  rapids_schema <- read_yaml(rapids_schema_file)
  devices <- participant_data$PHONE$DEVICE_IDS
  device_oss <- participant_data$PHONE$PLATFORMS
  device_oss <- replace(device_oss, device_oss == "multiple", "infer") # support multiple for retro compatibility
  validate_deviceid_platforms(devices, device_oss, participant_file)

  if(length(device_oss) == 1)
    device_oss <- rep(device_oss, length(devices))

  validate_expected_columns_mapping(stream_schema, rapids_schema, sensor, rapids_schema_file, stream_format)
  # ANDROID or IOS RAPIDS_COLUMN_MAPPINGS are guaranteed to be the same at this point (see validate_expected_columns_mapping function)
  expected_columns <- tolower(rapids_schema[[sensor]])
  participant_data <- setNames(data.frame(matrix(ncol = length(expected_columns), nrow = 0)), expected_columns)

  if(length(devices) == 0){
    warning("There were no PHONE device ids in this participant file:", participant_file)
    write_csv(participant_data, output_data_file)
    return()
  }

  container_functions <- load_container_script(stream_container)
  infer_device_os_container <- container_functions$infer_device_os
  pull_data_container <- container_functions$pull_data

  for(idx in seq_along(devices)){ #TODO remove length
    
    device <- devices[idx]
    message(paste0("\nProcessing ", sensor, " for ", device))
    device_os <- ifelse(device_oss[idx] == "infer", infer_device_os_container(data_configuration, device), device_oss[idx])
    validate_inferred_os(basename(stream_container), participant_file, device, device_os)
    
    if(!toupper(device_os) %in% names(stream_schema[[sensor]])){ # the current sensor is only available in a single OS (like PHONE_MESSAGES)
      warning(sensor, " data is not available for ", device_os, ". No data to download for ", device)
      next
    }

    os_table <- ifelse(length(tables) > 1, tables[[toupper(device_os)]], tables) # some sensor tables have a different name for android and ios    

    columns_to_download <- c(stream_schema[[sensor]][[toupper(device_os)]][["RAPIDS_COLUMN_MAPPINGS"]], stream_schema[[sensor]][[toupper(device_os)]][["MUTATION"]][["COLUMN_MAPPINGS"]])
    columns_to_download <- columns_to_download[(columns_to_download != "FLAG_TO_MUTATE")]
    data <- pull_data_container(data_configuration, device, sensor, os_table, columns_to_download)
    
    if(!setequal(columns_to_download, colnames(data)))
      stop(paste0("The pulled data for ", device, " does not have the expected columns (including [RAPIDS_COLUMN_MAPPINGS] and [MUTATE][COLUMN_MAPPINGS]). The container script returned [", paste(colnames(data), collapse=","),"] but the format mappings expected [",paste(columns_to_download, collapse=","), "]. The conainer script is: ", stream_container))
    
    renamed_data <- rename_columns(columns_to_download, data)
    
    mutation_scripts <- stream_schema[[sensor]][[toupper(device_os)]][["MUTATION"]][["SCRIPTS"]]
    mutated_data <- mutate_data(mutation_scripts, renamed_data, data_configuration)

    if(!setequal(expected_columns, colnames(mutated_data)))
      stop(paste0("The mutated data for ", device, " does not have the columns RAPIDS expects. The mutation script returned [", paste(colnames(mutated_data), collapse=","),"] but RAPIDS expected [",paste(expected_columns, collapse=","), "]. One ore more mutation scripts in [", sensor,"][MUTATION][SCRIPTS] are adding extra columns or removing or not adding the ones expected"))
    participant_data <- rbind(participant_data, mutated_data %>% distinct())
      
  }
  participant_data <- participant_data %>% arrange(timestamp)
  write_csv(participant_data, output_data_file)
}

pull_phone_data()