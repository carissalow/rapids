source("renv/activate.R")
library(tidyr)
library("dplyr", warn.conflicts = F)
library(tidyverse)
library(caret)
library(corrr)

rapids_cleaning <- function(sensor_data_files, provider){

    clean_features <- read.csv(sensor_data_files[["sensor_data"]], stringsAsFactors = FALSE)
    platforms <- read.csv(sensor_data_files[["platforms"]], stringsAsFactors = FALSE)
    impute_selected_event_features <- provider[["IMPUTE_SELECTED_EVENT_FEATURES"]]
    cols_nan_threshold <- as.numeric(provider[["COLS_NAN_THRESHOLD"]])
    drop_zero_variance_columns <- as.logical(provider[["COLS_VAR_THRESHOLD"]])
    rows_nan_threshold <- as.numeric(provider[["ROWS_NAN_THRESHOLD"]])
    data_yield_unit <- tolower(str_split_fixed(provider[["DATA_YIELD_FEATURE"]], "_", 4)[[4]])
    data_yield_column <- paste0("phone_data_yield_rapids_ratiovalidyielded", data_yield_unit)
    data_yield_ratio_threshold <- as.numeric(provider[["DATA_YIELD_RATIO_THRESHOLD"]])
    drop_highly_correlated_features <- provider[["DROP_HIGHLY_CORRELATED_FEATURES"]]

    # Impute selected event features
    if(as.logical(impute_selected_event_features$COMPUTE)){
        if(!"phone_data_yield_rapids_ratiovalidyieldedminutes" %in% colnames(clean_features)){
            stop("Error: RAPIDS provider needs to impute the selected event features based on phone_data_yield_rapids_ratiovalidyieldedminutes column, please set config[PHONE_DATA_YIELD][PROVIDERS][RAPIDS][COMPUTE] to True and include 'ratiovalidyieldedminutes' in [FEATURES].")
        }
        column_names <- colnames(clean_features)
        # Features which can be extracted from both Android and iOS devices
        selected_battery_features <- column_names[grepl("^phone_battery_rapids_", column_names)]
        selected_calls_features <- column_names[grepl("^phone_calls_rapids_.*_(count|distinctcontacts|sumduration|minduration|maxduration|meanduration|modeduration)", column_names)]
        selected_screen_features <- column_names[grepl("^phone_screen_rapids_(sumduration|maxduration|minduration|avgduration|countepisode)", column_names)]
        selected_wifi_connected_features <- column_names[grepl("^phone_wifi_connected_rapids_", column_names)]
        selected_androidios_columns <- c(selected_battery_features, selected_calls_features, selected_screen_features, selected_wifi_connected_features)
        clean_features[selected_androidios_columns][is.na(clean_features[selected_androidios_columns]) & (clean_features$phone_data_yield_rapids_ratiovalidyieldedminutes > impute_selected_event_features$MIN_DATA_YIELDED_MINUTES_TO_IMPUTE)] <- 0

        # Features which can only be extracted from Android devices
        selected_apps_features <- column_names[grepl("^phone_applications_foreground_rapids_(countevent|countepisode|minduration|maxduration|meanduration|sumduration)", column_names)]        
        selected_keyboard_features <- column_names[grepl("^phone_keyboard_rapids_(sessioncount|averagesessionlength|changeintextlengthlessthanminusone|changeintextlengthequaltominusone|changeintextlengthequaltoone|changeintextlengthmorethanone|maxtextlength|totalkeyboardtouches)", column_names)]
        selected_messages_features <- column_names[grepl("^phone_messages_rapids_.*_(count|distinctcontacts)", column_names)]
        selected_wifi_visible_features <- column_names[grepl("^phone_wifi_visible_rapids_", column_names)]
        selected_android_columns <- c(selected_apps_features, selected_keyboard_features, selected_messages_features, selected_wifi_visible_features)
        if(nrow(platforms) == 1){ # Single platform
            if(platforms$os == "android"){
                clean_features[selected_android_columns][is.na(clean_features[selected_android_columns]) & (clean_features$phone_data_yield_rapids_ratiovalidyieldedminutes > impute_selected_event_features$MIN_DATA_YIELDED_MINUTES_TO_IMPUTE)] <- 0
            }
        }else{ # Multiple platforms
            for(idx in 1:nrow(platforms)){
                if(platforms[idx, "os"] == "android"){
                    if(idx < nrow(platforms)){
                        selected_android_rows <- (clean_features$local_segment_start_datetime >= platforms[idx, "local_date_time"]) & (clean_features$local_segment_end_datetime <= platforms[idx + 1, "local_date_time"])
                    }else{
                        selected_android_rows <- clean_features$local_segment_start_datetime >= platforms[idx, "local_date_time"]
                    }
                    clean_features[selected_android_rows, selected_android_columns][is.na(clean_features[selected_android_rows, selected_android_columns]) & (clean_features[selected_android_rows, ]$phone_data_yield_rapids_ratiovalidyieldedminutes > impute_selected_event_features$MIN_DATA_YIELDED_MINUTES_TO_IMPUTE)] <- 0
                }
            }
        }
    }
    
    # Drop rows with the value of data_yield_column less than data_yield_ratio_threshold
    if(!data_yield_column %in% colnames(clean_features)){
        stop(paste0("Error: RAPIDS provider needs to clean data based on ", data_yield_column, " column, please set config[PHONE_DATA_YIELD][PROVIDERS][RAPIDS][COMPUTE] to True and include 'ratiovalidyielded", data_yield_unit, "' in [FEATURES]."))
    }
    clean_features <- clean_features %>% 
        filter(.[[data_yield_column]] >= data_yield_ratio_threshold)

    # Drop columns with a percentage of NA values above cols_nan_threshold
    if(nrow(clean_features))
        clean_features <- clean_features %>% select_if(~ sum(is.na(.)) / length(.) <= cols_nan_threshold )

    # Drop columns with zero variance
    if(drop_zero_variance_columns)
    clean_features <- clean_features %>% select_if(grepl("pid|local_segment|local_segment_label|local_segment_start_datetime|local_segment_end_datetime",names(.)) | sapply(., n_distinct, na.rm = T) > 1)

    # Drop highly correlated features
    if(as.logical(drop_highly_correlated_features$COMPUTE)){
        
        min_overlap_for_corr_threshold <- as.numeric(drop_highly_correlated_features$MIN_OVERLAP_FOR_CORR_THRESHOLD)
        corr_threshold <- as.numeric(drop_highly_correlated_features$CORR_THRESHOLD)

        features_for_corr <- clean_features %>% 
            select_if(is.numeric) %>% 
            select_if(sapply(., n_distinct, na.rm = T) > 1)

        valid_pairs <- crossprod(!is.na(features_for_corr)) >= min_overlap_for_corr_threshold * nrow(features_for_corr)

        if((nrow(features_for_corr) != 0) & (ncol(features_for_corr) != 0)){

            highly_correlated_features <- features_for_corr %>% 
                correlate(use = "pairwise.complete.obs", method = "spearman") %>% 
                column_to_rownames(., var = "term") %>% 
                as.matrix() %>% 
                replace(!valid_pairs | is.na(.), 0) %>% 
                findCorrelation(., cutoff = corr_threshold, verbose = F, names = T)

            clean_features <- clean_features[, !names(clean_features) %in% highly_correlated_features]
        
        }
    }

    # Drop rows with a percentage of NA values above rows_nan_threshold
    clean_features <- clean_features %>% 
        mutate(percentage_na =  rowSums(is.na(.)) / ncol(.)) %>% 
        filter(percentage_na <= rows_nan_threshold) %>% 
        select(-percentage_na)

    return(clean_features)
}

