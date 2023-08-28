add_context_boundaries <- function(x, context) {
    co <- x$coordinates
    assert(
        co$end - co$start + 1 >= context * 2,
        msg=paste("CDS is not large enough for ", context, " bps of context")
    )
    co$context5_start <- 1
    co$context5_end <- co$start + context - 1
    co$context3_start <- co$end - context + 1
    co$context3_end <- co$seqlength
    # CDS boundaries relative to context
    co$context5_cr <- co$start
    co$context3_cr <- context

    x$coordinates <- co

    return(x)
}
