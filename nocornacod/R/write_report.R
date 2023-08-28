# write_report.R

write_report <- function(x, output) {
    report <- list(
        coordinates=x$coordinates
    )

    if (!is.null(x$cpt5)) {
        report[["changepoint5"]] <- x$cpt5
        report[["mean_ncrna5"]] <- mean(x$coordinates$ncrna5_len)
    }

    if (!is.null(x$cpt3)) {
        report[["changepoint3"]] <- x$cpt3
        report[["mean_ncrna3"]] <- mean(x$coordinates$ncrna3_len)
    }

    if (!is.null(x$fold5)) {
        report[["energy5"]] <- as.list(x$fold5$energies)
    }
    if (!is.null(x$fold3)) {
        report[["energy3"]] <- as.list(x$fold3$energies)
    }

    jsonlite::write_json(report, output, pretty=TRUE)
}
