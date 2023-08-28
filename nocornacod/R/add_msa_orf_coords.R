add_msa_relative_coords <- function(x) {
    midx <- Biostrings::vmatchPattern("N", x$msa5, fixed=FALSE)
    x$coordinates$context5_mr <- get_nth_base(midx, x$coordinates$context5_cr)
    midx <- Biostrings::vmatchPattern("N", x$msa3, fixed=FALSE)
    x$coordinates$context3_mr <- get_nth_base(midx, x$coordinates$context3_cr)

    return(x)
}

get_nth_base <- function(midx, n) {
    x <- mapply(
        \(i, j) if (length(i) >= j) IRanges::start(i[j]) else NA,
        midx, n
    )
    return(x)
}
