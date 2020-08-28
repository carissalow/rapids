library("stringr")
filter_data_by_segment <- function(data, day_segment){
    # Filter the rows that belong to day_segment, and put the segment full name in a new column for grouping
    date_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2}"
    hour_regex = "[0-9]{2}:[0-9]{2}:[0-9]{2}"
    data <- data %>% 
        filter(grepl(paste0("\\[", day_segment, "#"), assigned_segments)) %>% 
        mutate(local_segment = str_extract(assigned_segments, paste0("\\[", day_segment, "#", date_regex, "#", hour_regex, "#", date_regex, "#", hour_regex, "\\]")),
                local_segment = str_sub(local_segment, 2, -2)) # get rid of first and last character([])
    return(data)
}
rapids_log_tag <- "RAPIDS:"