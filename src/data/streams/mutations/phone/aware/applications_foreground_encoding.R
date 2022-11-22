source("renv/activate.R")
library("dplyr", warn.conflicts = FALSE)

convert_encoding_to_utf8 <- function(phone_applications_foreground) {
    phone_applications_foreground <- phone_applications_foreground %>%
        mutate(application_name = stringi::stri_encode(application_name, to = "UTF-8"))

    return(phone_applications_foreground)
}

main <- function(data, stream_parameters) {
    return(convert_encoding_to_utf8(data))
}