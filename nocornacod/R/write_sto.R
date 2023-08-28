write_sto <- function(msa, ss, output, seq_width=80L) {
    con <- file(output, "wt")
    on.exit(close(con), add=TRUE)

    m <- regexec("^([^[:space:]]+)[[:space:]]*(.*)$", names(msa))
    parts <- regmatches(names(msa), m)
    ids <- vapply(parts, `[[`, character(1), 2)
    descs <- vapply(
        parts, \(x) ifelse(length(x) > 2, x[[3]], ""), character(1)
    )

    id_width <- max(nchar(ids))
    cat("# STOCKHOLM 1.0\n\n", file=con)
    has_desc <- nzchar(descs)
    if (any(has_desc)) {
        fmt <- paste0("#=GS %-", id_width, "s DE %s")
        cat(
            sprintf(fmt, ids[has_desc], descs[has_desc]),
            file=con, sep="\n"
        )
        cat("\n", file=con)
    }

    msa_width <- Biostrings::width(msa)[[1]]
    i <- 1L
    fmt <- paste0("%-", id_width, "s %s")
    while (i <= msa_width) {
        j <- i + seq_width - 1L
        j <- ifelse(j <= msa_width, j, msa_width)
        block <- Biostrings::subseq(msa, start=i, end=j) |>
            as.character()
        ss_block <- substr(ss, i, j)
        cat(sprintf(fmt, ids, block), file=con, sep="\n")
        cat(sprintf(fmt, "#=GC SS_cons", ss_block), file=con, sep="\n")
        cat("\n", file=con)
        i <- i + seq_width
    }

    cat("//\n", file=con)
}
