source("renv/activate.R")
library("dplyr", warn.conflicts = FALSE)

 convert_encoding_to_utf8 <- function(phone_keyboard) {
     phone_keyboard <- phone_keyboard %>%
         mutate(
            package_name = stringi::stri_encode(package_name, to = "UTF-8"),
            before_text = stringi::stri_encode(before_text, to = "UTF-8"),
            current_text = stringi::stri_encode(current_text, to = "UTF-8")
        )

     return(phone_keyboard)
 }

 main <- function(data, stream_parameters) {
     return(convert_encoding_to_utf8(data))
 }