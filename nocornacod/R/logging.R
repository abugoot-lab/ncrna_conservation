get_time_str <- function() {
    return(paste0(
        "[", format(Sys.time(), format="%Y-%m-%dT%H:%M:%S%z"), "]"
    ))
}

msgf <- function(...) {
    return(paste0(get_time_str(), " ", ...))
}
