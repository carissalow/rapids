source("renv/activate.R")

library(RMariaDB)
library(stringr)
library(purrr)
library(readr)
library("dplyr", warn.conflicts = F)

config <- snakemake@params[["config"]]
group <- config$SOURCE$DATABASE_GROUP
timezone <- config$SOURCE$TIMEZONE
phone_device_id_column = "device_id"
fitbit_device_id_column = "fitbit_id"
empatica_device_id_column = "empatica_id"
add_phone_section = config$PHONE_SECTION$ADD
add_fitbit_section = config$FITBIT_SECTION$ADD
add_empatica_section = config$EMPATICA_SECTION$ADD
phone_ignored = config$PHONE_SECTION$IGNORED_DEVICE_IDS
fitbit_ignored = config$FITBIT_SECTION$IGNORED_DEVICE_IDS
empatica_ignored = config$EMPATICA_SECTION$IGNORED_DEVICE_IDS

rmysql.settingsfile <- "./.env"

participants <- read_csv(config$CSV_FILE_PATH, col_types=cols_only(device_id="c",pid="c",label="c",platform="c",
                          start_date=col_datetime(),end_date=col_datetime(),fitbit_id="c",empatica_id="c")) %>% 
                          mutate(start_date = as.character(start_date), end_date = as.character(end_date)) # we read as date to validate format
participants <- participants %>% 
mutate(!!phone_device_id_column := str_replace_all(!!rlang::sym(phone_device_id_column), ";",","),
        !!fitbit_device_id_column := str_replace_all(!!rlang::sym(fitbit_device_id_column), ";",","),
        !!empatica_device_id_column := str_replace_all(!!rlang::sym(empatica_device_id_column), ";",","),
        platform = str_replace_all(platform, ";",","),
        !!phone_device_id_column := if_else(!!rlang::sym(phone_device_id_column) %in% phone_ignored, NA_character_, !!rlang::sym(phone_device_id_column)),
        !!empatica_device_id_column := if_else(!!rlang::sym(empatica_device_id_column) %in% empatica_ignored, NA_character_, !!rlang::sym(empatica_device_id_column)),
        !!fitbit_device_id_column := if_else(!!rlang::sym(fitbit_device_id_column) %in% fitbit_ignored, NA_character_, !!rlang::sym(fitbit_device_id_column)))

dir.create(file.path("./data/external/participant_files/"))

participants %>%
  pwalk(function(add_phone_section, add_fitbit_section, phone_device_id_column, fitbit_device_id_column, ...) {
    empty_phone <- c("PHONE:", "  DEVICE_IDS:", "  PLATFORMS:","  LABEL:", "  START_DATE:", "  END_DATE:")
    empty_fitbit <- c("FITBIT:", "  DEVICE_IDS:", "  LABEL:", "  START_DATE:", "  END_DATE:")
    empty_empatica <- c("EMPATICA:", "  LABEL:", "  START_DATE:", "  END_DATE:")
    row <- tibble(...)
    lines <- c()
    start_date = if_else(is.na(row$start_date), "", row$start_date)
    end_date = if_else(is.na(row$end_date), "", row$end_date)

    if(add_phone_section == TRUE && !is.na(row[phone_device_id_column])){
      lines <- append(lines, c("PHONE:", paste0("  DEVICE_IDS: [",row[phone_device_id_column],"]"), paste0("  PLATFORMS: [",row$platform,"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", start_date), paste("  END_DATE:", end_date)))
    }else
      lines <- append(lines, empty_phone)
    
    if(add_fitbit_section == TRUE && !is.na(row[fitbit_device_id_column])){
      lines <- append(lines, c("FITBIT:", paste0("  DEVICE_IDS: [",row[fitbit_device_id_column],"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", start_date), paste("  END_DATE:", end_date)))
    } else
      lines <- append(lines, empty_fitbit)

    if(add_empatica_section == TRUE && !is.na(row[empatica_device_id_column])){
      lines <- append(lines, c("EMPATICA:", paste0("  DEVICE_IDS: [",row[empatica_device_id_column],"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", start_date), paste("  END_DATE:", end_date)))
    } else
      lines <- append(lines, empty_empatica)
    lines <- append(lines, "\n")

    file_connection <- file(paste0("./data/external/participant_files/", row$pid, ".yaml"))
    writeLines(lines, file_connection)
    close(file_connection)

  }, add_phone_section, add_fitbit_section, phone_device_id_column, fitbit_device_id_column, empatica_device_id_column)

file_lines <-readLines("./config.yaml")
for (i in 1:length(file_lines)){
  if(startsWith(file_lines[i], "PIDS:")){
    file_lines[i] <- paste0("PIDS: ['", paste(participants$pid, collapse = "', '"), "']")
  }
}
writeLines(file_lines, con = "./config.yaml") 