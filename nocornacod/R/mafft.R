msa <- function(x, tmpdir, threads) {
    co <- x$coordinates
    s <- x$sequences
    fivep <- Biostrings::subseq(
        s, start=co$context5_start, end=co$context5_end
    )
    threep <- Biostrings::subseq(
        s, start=co$context3_start, end=co$context3_end
    )

    x$msa5 <- mafft(fivep, tmpdir, threads)
    x$msa3 <- mafft(threep, tmpdir, threads)

    return(x)
}

mafft <- function(seqs, tmpdir, threads) {
    tmp_in <- tempfile("mafft_in_", tmpdir=tmpdir, fileext=".fna")
    tmp_out <- tempfile("mafft_out_", tmpdir=tmpdir, fileext=".fna")
    on.exit(suppressWarnings(file.remove(c(tmp_in, tmp_out))), add=TRUE)

    Biostrings::writeXStringSet(seqs, tmp_in)
    rtn <- system2(
        "mafft",
        args=c(
            "--auto", "--thread", threads, "--nuc", "--quiet",
            tmp_in
          ),
        stdout=tmp_out
    )

    if (rtn != 0) {
        stop("mafft returned a non-zero exit code", call.=FALSE)
    }

    return(Biostrings::readDNAStringSet(tmp_out))
}
