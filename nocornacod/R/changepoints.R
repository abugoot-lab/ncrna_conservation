cpts_pelt <- function(x) {
    results <- tryCatch(
        warning=function(cnd) NA_real_,
        error=function(cnd) NA_real_,
        changepoint::cpt.mean(
            x, method="PELT", penalty="CROPS", pen.value=c(4, 1500)
        )
    )

    if (suppressWarnings(is.na(results))) {
        return(NULL)
    }

    # The first row should always have the most change points
    cpts <- changepoint::cpts.full(results)[1, ]
    if (all(is.na(cpts))) {
        return(NULL)
    }

    return(as.integer(cpts[!is.na(cpts)]))
}

# The change point detection algorithm used does not report an optimal number
# of change points so we pick the smallest range that fits the criteria of
# 1. Minimum ncRNA length
# 2. Minimum average conservation score
predict_cpts <- function(x, min_ncrna_len, min_avg_score) {
    co <- x$coordinates
    msa5 <- x$msa5
    msa3 <- x$msa3
    msa5_scores <- x$msa5_scores
    msa3_scores <- x$msa3_scores
    sink(file=nullfile())
    cpt5 <- cpts_pelt(msa5_scores)
    cpt3 <- cpts_pelt(msa3_scores)
    sink(file=NULL)

    # Only keep change points that would result in ncRNA lengths greater than
    # the minimum
    cpt5_pass <- FALSE
    cpt3_pass <- FALSE
    if (!is.null(cpt5)) {
        cpt5_tmp <- Filter(\(x) all(x < co$context5_mr), cpt5)
        if (length(cpt5_tmp) > 0) {
            ncrna5_len_tmp <- lapply(
                cpt5_tmp, \(x) cnt_bases(msa5, x, co$context5_mr)
            )
            # The score on the 5' end is always from the changepoint to the
            # end of the scores
            ncrna5_score_pass <- vapply(
                cpt5_tmp,
                \(x) mean(msa5_scores[-seq_len(x - 1)]) >= min_avg_score,
                logical(1)
            )
            ncrna5_len_pass <- vapply(
                ncrna5_len_tmp, \(x) all(x >= min_ncrna_len), logical(1)
            )
            ncrna5_pass <- ncrna5_score_pass & ncrna5_len_pass
            if (any(ncrna5_pass)) {
                cpt5_pass <- TRUE
                i5 <- which.max(cpt5_tmp[ncrna5_pass])
                cpt5 <- cpt5_tmp[ncrna5_pass][[i5]]
                ncrna5_len <- ncrna5_len_tmp[ncrna5_pass][[i5]]
            }
        }
    }

    if (!is.null(cpt3)) {
        cpt3_tmp <- Filter(\(x) all(x > co$context3_mr), cpt3)
        if (length(cpt3_tmp) > 0) {
            ncrna3_len_tmp <- lapply(
                cpt3_tmp, \(x) cnt_bases(msa3, co$context3_mr, x)
            )
            # The score on the 3' end is always from the start of the
            # scores to the change point
            ncrna3_score_pass <- vapply(
                cpt3_tmp,
                \(x) mean(msa3_scores[1:x]) >= min_avg_score,
                logical(1)
            )
            ncrna3_len_pass <- vapply(
                ncrna3_len_tmp, \(x) all(x >= min_ncrna_len), logical(1)
            )
            ncrna3_pass <- ncrna3_score_pass & ncrna3_len_pass
            if (any(ncrna3_pass)) {
                cpt3_pass <- TRUE
                i3 <- which.min(cpt3_tmp[ncrna3_pass])
                cpt3 <- cpt3_tmp[ncrna3_pass][[i3]]
                ncrna3_len <- ncrna3_len_tmp[ncrna3_pass][[i3]]
            }
        }
    }

    if (cpt5_pass) {
        co$ncrna5_len <- ncrna5_len
        x$cpt5 <- cpt5
    } else {
        co$ncrna5_len <- NA_integer_
    }

    if (cpt3_pass) {
        co$ncrna3_len <- ncrna3_len
        x$cpt3 <- cpt3
    } else {
        co$ncrna3_len <- NA_integer_
    }
    x$coordinates <- co

    return(x)
}

# Count number of bases between two points in a sequence, inclusive
cnt_bases <- function(x, start, end) {
    s <- Biostrings::subseq(x, start, end)
    m <- Biostrings::vmatchPattern("N", s, fixed=FALSE)
    return(lengths(m))
}
