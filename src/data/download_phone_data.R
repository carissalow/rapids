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

validate_inferred_os <- function(source_download_file, participant_file, device, device_os){
  if(!is.na(device_os) && device_os != "android" && device_os != "ios")
    stop(paste0("We tried to infer the OS for ", device, " but 'infer_device_os' function inside '",source_download_file,"' returned '",device_os,"' instead of 'android' or 'ios'. You can assign the OS manually in the participant file or report this bug on GitHub.\nParticipant file ", participant_file))
}

mutate_data <- function(scripts, data){
  for(script in scripts){
    if(grepl("\\.(R)$", script)){
      myEnv <- new.env()    
      source(script, local=myEnv)
      attach(myEnv, name="sourced_scripts_rapids")
      if(exists("main", myEnv)){
        message(paste("Applying mutation script", script))
        data <- main(data)
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
        data <- script_functions$main(data)
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

validate_expected_columns_mapping <- function(schema, rapids_schema, sensor, rapids_schema_file){
  android_columns <- names(schema[[sensor]][["ANDROID"]][["COLUMN_MAPPINGS"]])
  android_columns <- android_columns[(android_columns != "FLAG_AS_EXTRA")]

  ios_columns <- names(schema[[sensor]][["IOS"]][["COLUMN_MAPPINGS"]])
  ios_columns <- ios_columns[(ios_columns != "FLAG_AS_EXTRA")]
  rapids_columns <- rapids_schema[[sensor]]

  if(is.null(rapids_columns))
    stop(paste(sensor, " columns are not listed in RAPIDS' column specification. If you are adding support for a new phone sensor, add any mandatory columns in ", rapids_schema_file))
  if(length(setdiff(rapids_columns, android_columns)) > 0)
    stop(paste(sensor," mappings are missing one or more mandatory columns for ANDROID. The missing column mappings are for ", paste(setdiff(rapids_columns, android_columns), collapse=","),"in", rapids_schema_file))
  if(length(setdiff(rapids_columns, ios_columns)) > 0)
    stop(paste(sensor," mappings are missing one or more mandatory columns for IOS. The missing column mappings are for ", paste(setdiff(rapids_columns, ios_columns), collapse=","),"in", rapids_schema_file))
  if(length(setdiff(android_columns, rapids_columns)) > 0)
    stop(paste(sensor," mappings have one or more columns than required for ANDROID, add them as FLAG_AS_EXTRA instead. The extra column mappings are for ", paste(setdiff(android_columns, rapids_columns), collapse=","),"in", rapids_schema_file))
  if(length(setdiff(ios_columns, rapids_columns)) > 0)
    stop(paste(sensor," mappings have one or more columns than required for IOS, add them as FLAG_AS_EXTRA instead. The extra column mappings are for ", paste(setdiff(ios_columns, rapids_columns), collapse=","),"in", rapids_schema_file))
}

download_phone_data <- function(){
  participant_file <- snakemake@input[["participant_file"]]
  source_schema_file <- snakemake@input[["source_schema_file"]]
  rapids_schema_file <- snakemake@input[["rapids_schema_file"]]
  source_download_file <- snakemake@input[["source_download_file"]]
  data_configuration <- snakemake@params[["data_configuration"]]
  tables <- snakemake@params[["tables"]]
  sensor <- toupper(snakemake@params[["sensor"]])
  output_data_file <- snakemake@output[[1]]

  source(source_download_file)

  participant_data <- read_yaml(participant_file)
  schema <- read_yaml(source_schema_file)
  rapids_schema <- read_yaml(rapids_schema_file)
  devices <- participant_data$PHONE$DEVICE_IDS
  device_oss <- participant_data$PHONE$PLATFORMS
  device_oss <- replace(device_oss, device_oss == "multiple", "infer") # support multiple for retro compatibility
  validate_deviceid_platforms(devices, device_oss, participant_file)

  if(length(device_oss) == 1)
    device_oss <- rep(device_oss, length(devices))

  validate_expected_columns_mapping(schema, rapids_schema, sensor, rapids_schema_file)
  # ANDROID or IOS COLUMN_MAPPINGS are guaranteed to be the same at this point (see validate_expected_columns_mapping function)
  expected_columns <- tolower(names(schema[[sensor]][["ANDROID"]][["COLUMN_MAPPINGS"]]))
  expected_columns <- expected_columns[(expected_columns != "flag_extra")]
  participant_data <- setNames(data.frame(matrix(ncol = length(expected_columns), nrow = 0)), expected_columns)

  for(idx in seq_along(devices)){ #TODO remove length
    
    device <- devices[idx]
    message(paste0("\nProcessing ", sensor, " for ", device))
    device_os <- ifelse(device_oss[idx] == "infer", infer_device_os(data_configuration, device), device_oss[idx])
    validate_inferred_os(basename(source_download_file), participant_file, device, device_os)
    os_table <- ifelse(length(tables) > 1, tables[[toupper(device_os)]], tables) # some sensor tables have a different name for android and ios    

    columns_to_download <- schema[[sensor]][[toupper(device_os)]][["COLUMN_MAPPINGS"]]
    columns_to_download <- columns_to_download[(columns_to_download != "FLAG_TO_MUTATE")]
    data <- download_data(data_configuration, device, os_table, columns_to_download)
    
    # Rename all COLUMN_MAPPINGS except those mapped as FLAG_AS_EXTRA or FLAG_TO_MUTATE
    columns_to_rename <- schema[[sensor]][[toupper(device_os)]][["COLUMN_MAPPINGS"]]
    columns_to_rename <- (columns_to_rename[(columns_to_rename != "FLAG_TO_MUTATE" & names(columns_to_rename) != "FLAG_AS_EXTRA")])
    renamed_data <- rename_columns(columns_to_rename, data)
    
    mutation_scripts <- schema[[sensor]][[toupper(device_os)]][["MUTATION_SCRIPTS"]]
    mutated_data <- mutate_data(mutation_scripts, renamed_data)

    if(length(setdiff(expected_columns, colnames(mutated_data))) > 0)
      stop(paste("The mutated data for", device, "is missing these columns expected by RAPIDS: [", paste(setdiff(expected_columns, colnames(mutated_data)), collapse=","),"]. One ore more mutation scripts in [", sensor,"][",toupper(device_os), "]","[MUTATION_SCRIPTS] are removing or not adding these columns"))
    participant_data <- rbind(participant_data, mutated_data)
      
  }

  write_csv(participant_data, output_data_file)
}

download_phone_data()