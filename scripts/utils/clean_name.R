clean_name <- function(x) {
  x %>%
    # Convert to ASCII, replacing special chars (including �) with _
    iconv(to = "ASCII//TRANSLIT", sub = "_") %>%  
    # Remove all non-alphanumeric characters except underscores
    str_replace_all("[^[:alnum:]_]", "") %>%  
    # Remove trailing underscores
    str_replace_all("_+$", "")
}