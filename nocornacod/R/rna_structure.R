# Sample n sequences from s with weights w
sample_seqs <- function(s, w, n) {
    top <- order(w, decreasing=TRUE)
    x <- top[1:n]
    s <- s[x]
    return(s)
}

predict_rna_struct <- function(x, pad, max_fold, outdir_path) {
    msa5 <- NULL
    msa3 <- NULL
    if (!is.null(x$cpt5) && all(x$cpt5 - pad > 0)) {
        msa5 <- Biostrings::subseq(
            x$msa5, start=x$cpt5 - pad, end=Biostrings::width(x$msa5)
        )
        if (max_fold < length(msa5)) {
            msa5 <- sample_seqs(msa5, attr(x$msa5_scores, "weights"), max_fold) 
        }
    }

    if (!is.null(x$cpt3) && all(x$cpt3 + pad <= Biostrings::width(x$msa3))) {
        msa3 <- Biostrings::subseq(
            x$msa3, start=1, end=x$cpt3 + pad
        )
        if (max_fold < length(msa3)) {
            msa3 <- sample_seqs(msa3, attr(x$msa3_scores, "weights"), max_fold)
        }
    }

    if (!is.null(msa5)) {
        message(msgf("Folding 5'"))
        out <- file.path(outdir_path, "R-scape_out_5")
        x$fold5 <- do_predict(msa5, out)
    }

    if (!is.null(msa3)) {
        message(msgf("Folding 3'"))
        out <- file.path(outdir_path, "R-scape_out_3")
        x$fold3 <- do_predict(msa3, out)
    }

    return(x)
}

do_predict <- function(msa, outdir) {
    dir.create(outdir, showWarnings=FALSE)
    msa_fa <- tempfile("RNAalifold_in", tmpdir=outdir, fileext=".mfa")
    
    Biostrings::writeXStringSet(msa, msa_fa)
    rnaalifold_results <- tryCatch(
        error=function(cnd) {
            message(conditionMessage(cnd))
            unlink(outdir, recursive=TRUE)
            return(NULL)
        },
        {
            rnaalifold(msa_fa, outdir)
        }
    )
    file.remove(msa_fa)
    
    if (is.null(rnaalifold_results)) {
        return(NULL)
    }

    msa_sto <- tempfile("RNAalifold_out", tmpdir=outdir, fileext=".sto")
    write_sto(msa, rnaalifold_results$structure, msa_sto)
    rscape_done <- tryCatch(
        error=function(cnd) {
            message(conditionMessage(cnd))
            unlink(outdir, recursive=TRUE)
            return(FALSE)
        },
        {
            rscape(msa_sto, outdir)        
            TRUE
        }
    )
    file.remove(msa_sto)

    if (!rscape_done) {
        return(NULL)
    }

    return(list("rscape"=outdir, "energies"=rnaalifold_results$energies))
}

rscape <- function(msa, outdir) {
    rtn <- system2(
        "R-scape",
        args=c(
            "-s", "--outdir", outdir, "--outname", "results", msa
        ),
        stdout=NULL, stderr=NULL
    )

    if (rtn != 0 || !file.exists(outdir)) {
        stop("R-scape failed", call.=FALSE)
    }

    return(outdir)
}
