source("renv/activate.R")
library("dplyr", warn.conflicts = FALSE)

 convert_encoding_to_utf8 <- function(phone_bluetooth) {
     phone_bluetooth <- phone_bluetooth %>%
         mutate(bt_name = stringi::stri_encode(bt_name, to = "UTF-8"))

     return(phone_bluetooth)
 }

 main <- function(data, stream_parameters) {
     return(convert_encoding_to_utf8(data))
 }