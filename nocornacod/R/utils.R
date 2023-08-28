is_string <- function(x) {
    return(is.character(x) && length(x) == 1)
}

is_scalar <- function(x) {
    if (anyNA(x)) {
        return(FALSE)
    }

    if (!(typeof(x) %in% c("integer", "double"))) {
        return(FALSE)
    }

    if (any(is.infinite(x))) {
        return(FALSE)
    }

    if (length(x) != 1) {
        return(FALSE)
    }

    return(TRUE)
}

is_flag <- function(x) {
    if (anyNA(x)) {
        return(FALSE)
    }

    if (!is.logical(x)) {
        return(FALSE)
    }

    if (length(x) != 1) {
        return(FALSE)
    }

    return(TRUE)
}

is_integerish <- function(x) {
    if (anyNA(x)) {
        return(FALSE)
    }

    if (!(typeof(x) %in% c("integer", "double"))) {
        return(FALSE)
    }

    if (!is.vector(x)) {
        return(FALSE)
    }

    if (any(is.infinite(x))) {
        return(FALSE)
    }

    if (!all(x == trunc(x))) {
        return(FALSE)
    }

    return(TRUE)
}

assert <- function(..., msg=NULL) {
    exprs <- match.call(expand.dots=FALSE)$...
    for (exp in exprs) {
        result <- eval.parent(exp)
        if (!is.logical(result)) {
            stop("All arguments must evaluate to a logical")
        }

        if (!all(result)) {
            if (is.null(msg)) {
                stop(deparse(substitute(exp)), " is FALSE", call.=FALSE)
            } else {
                stop(msg, call.=FALSE)
            }
        }
    }

    return(invisible(TRUE))
}

prepare_outdir <- function(path, force) {
    if (dir.exists(path)) {
        if (force) {
            unlink(path, recursive=TRUE)
        } else {
            stop("Output directory exists", call.=FALSE)
        }
    }

    made_dir <- dir.create(path)
    if (!made_dir) {
        stop("Failed to create output directory", call.=FALSE)
    }

    return(invisible(TRUE))
}

cleanup <- function(...) {
    paths <- Filter(file.exists, list(...)) |> unlist()
    unlink(paths, recursive=TRUE)
}
