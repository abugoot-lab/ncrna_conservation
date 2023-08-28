DNA_ALPHABET <- c("A", "C", "G", "T", "-")

score_conservation <- function(x, smooth) {
    msa5 <- as.matrix(x$msa5)
    msa3 <- as.matrix(x$msa3)
    assert(
        smooth < ncol(msa5), smooth < ncol(msa3),
        msg="Smoothing window is too large"
    )
    x$msa5_scores <- valdar01(msa5) |> smooth_scores(smooth)
    x$msa3_scores <- valdar01(msa3) |> smooth_scores(smooth)

    return(x)
}

# The scoring function is one adapted from Valdar and Thornton (2001)
# which describes a scoring function for protein multiple sequence
# alignments. This function only scores alignments of DNA. The function
# calculates a score between 0 and 1 for each column in the MSA, reflecting
# the conservation of bases in that column, where 0 is no conservation
# (random) and 1 is complete conservation. Only the characters, A, T, C, G,
# N, and - are permitted in the alignment. All other bases will be converted
# to N and treated as not matching any other base.
#
# Valdar, WSJ. and Thornton, JM. (2001), Protein-protein interfaces: Analysis
# of amino acid conservation in homodimers. Proteins, 42:108-124.
# 
# Valdar WSJ. (2002), Scoring residue conservation. Proteins, 48:227-241.
#
# https://www.ebi.ac.uk/thornton-srv/databases/valdarprograms/scorecons_server_help.html
valdar01 <- function(msa) {
    msa[!(msa %in% DNA_ALPHABET)] <- "N"
    weights <- valdar_weights(msa)
    # If all weights are 0, then the sequences are all identical. We could,
    # in this case, weight all the sequences equally, but in the context of
    # the pipeline, it would not make sense to proceed because change point
    # detection would fail. Also, the computation of lambda would result in
    # NaN.
    assert(!all(weights == 0), msg="Sequences are identical. Can't proceed")
    lambda <- valdar_lambda(weights)
    scores <- apply(msa, 2, valdar_col, lambda=lambda, w=weights)
    attr(scores, "weights") <- weights
    attr(scores, "lambda") <- lambda

    return(scores)
}

# The rolling mean is symmetric and will result in NA values at the tails of
# resulting vector. In order to make it easier to calculate change points
# and align the MSA coordinates with the score coordinates, we pad the tails
# by repeating the scores at each end.
# TODO investigate issues with a smaller smoothing window
smooth_scores <- function(x, k) {
    a <- attributes(x)
    if (k > length(x)) {
        warning("Rounding smoothing window to width of MSA")
        k <- length(x)
    }

    scores <- rolling_avg(x, k)
    if (length(x) > k) {
        no_na <- which(!is.na(scores))
        start_tail <- scores[no_na[[1]]]
        end_tail <- scores[no_na[[length(no_na)]]]
        scores[seq_len(no_na[[1]] - 1)] <- start_tail
        scores[seq.int(no_na[[length(no_na)]], length(scores))] <- end_tail
    }
    attributes(scores) <- a

    return(scores)
}

valdar_col <- function(col, lambda, w) {
    x <- 0
    for (i in seq.int(1, length(col) - 1, 1)) {
        for (j in seq.int(i + 1, length(col), 1)) {
            if (col[i] %in% c("-", "N") | col[j] %in% c("-", "N")) {
                next
            }
            if (col[i] == col[j]) {
                x <- x + (w[i] * w[j]) 
            }
        }
    }

    return(x * lambda)
}

# The similarity is given as the number of matching non-N bases, divided by the
# number of alignment positions where there is a base in at least one of the
# sequences. The distance is then one minus the similarity. The original
# distance measure was designed for protein sequences and used a transformed
# scoring matrix to measure the similarity between two residues. With DNA
# alignments, there is no scoring matrix, only matches and mismatches.
valdar_dist <- function(x, y) {
    no_gaps <- x != "-" | y != "-"
    no_n <- x != "N" & y != "N"
    scores <- x[no_gaps & no_n] == y[no_gaps & no_n]
    return(1 - (sum(scores) / length(no_gaps)))
}

# The weight for a sequence is given as the average distance between that
# sequence and all other sequences in the alignment.
valdar_weights <- function(x) {
    n <- nrow(x)
    w <- double(n)
    for (i in seq_len(n)) {
        total <- 0
        for (j in seq_len(n)) {
            if (i == j) {
                next
            }
            total <- total + valdar_dist(x[i, ], x[j, ])
        }
        w[i] <- total / (n - 1)
    }

    return(w)
}

valdar_lambda <- function(x) {
    lambda <- 0
    for (i in seq.int(1, length(x) - 1, 1)) {
        for (j in seq.int(i + 1, length(x), 1)) {
            lambda <- lambda + (x[i] * x[j])
        }
    }

    return(1 / lambda)
}

rolling_avg <- function(x, k) {
    # Don't quite understand how this works
    avgs <- stats::filter(x, rep_len(1/k, k), sides=2)
    return(as.double(avgs))
}
