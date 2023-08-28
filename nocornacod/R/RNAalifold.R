#' Run RNAalifold
#' @noRd
rnaalifold <- function(file, tmpdir) {
    out <- tempfile(pattern="RNAalifold", tmpdir=tmpdir)
    on.exit(try(file.remove(out), silent=TRUE), add=TRUE)

    rtn <- system2(
        "RNAalifold",
        args=c(
            "--quiet", "--jobs=1", "--noPS", "--input-format=F", "--noLP",
            file
        ), stdout=out
    )

    if (rtn != 0) {
        stop("No output from RNAalifold", call.=FALSE)
    }

    ss <- parse_rnaalifold_output(out)

    return(ss)
}

#' Parse output of RNAalifold
#'
#' This only works if the program was run with a single alignment because the
#' assumption is that there is only the output for a single fold prediction.
#' In this case, the output is two lines: the first is the consensus sequence
#' of the alignment and the second is the secondary structure followed by the
#' folding energies.
#' @noRd
parse_rnaalifold_output <- function(file) {
    output <- readLines(file)
    stopifnot(length(output) == 2)

    consensus <- output[[1]]
    m <- regexec("^(\\S+)\\s+\\((.*)\\)$", output[[2]])
    parts <- regmatches(output[[2]], m)
    # Should have the original line, the fold and the energies
    stopifnot(length(parts[[1]]) == 3)

    ss <- parts[[1]][[2]]
    energies <- strsplit(parts[[1]][[3]], split="( *= *)|( *\\+ *)") |>
        (\(x) x[[1]])() |>
        as.double() |>
        as.list()
    names(energies) <- c("total", "fold", "cov")

    rtn <- list("consenus"=consensus, "structure"=ss, "energies"=energies)
    return(rtn)
}
