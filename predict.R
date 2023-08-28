#!/usr/bin/env Rscript

################################################################################
# predict.R
# This script implements a pipeline to predict ncRNA regions surrounding a
# protein CDS by scoring the conservation of an alignment of the regions. The
# alignment is scored by computing a per-column base conservation score using
# a sum-of-pairs strategy. If a conserved non-coding region is found, the
# pipeline also predicts the consensus secondary structure.
################################################################################

suppressPackageStartupMessages(library(optparse))
parser <- optparse::OptionParser(
  usage="usage: %prog [options] <FASTA> <COORDINATES> <OUTPUT>",
  description="Predict non-coding RNA" 
)
parser <- optparse::add_option(
  parser, c("-t", "--threads"), action="store", type="integer",
  dest="threads", default=1L, help="Number of threads to use [1]",
  metavar="N"
)
parser <- optparse::add_option(
  parser, c("--inner_pad"), action="store", type="integer",
  dest="inner_pad", default=100L, metavar="N",
  help=paste0(
    "Number of bases into the ORF to pad alignment region when predicting",
    " structure [100]"
  )
)
parser <- optparse::add_option(
  parser, c("--outer_pad"), action="store", type="integer",
  dest="outer_pad", default=11L, metavar="N",
  help=paste0(
    "Number of bases outside the conserved region to pad alignment region",
    " when predicting structure [11]"
  )
)
parser <- optparse::add_option(
  parser, c("--min_ncrna_len"), action="store", type="integer",
  dest="min_ncrna_len", default=13L, metavar="N",
  help=paste0(
    "Number of bases outside CDS that must be conserved",
    " to be considered ncNRA [13]"
  )
) 
parser <- optparse::add_option(
    parser, c("--min_score"), action="store", type="double",
    dest="min_score", default=0.5, metavar="R",
    help="Minimum score of ncRNA in MSA [0.5]"
)
parser <- optparse::add_option(
  parser, c("--smooth_win"), action="store", type="integer",
  dest="smooth_win", default=27L, metavar="N",
  help="Smoothing window for conservation score [27]"
)
parser <- optparse::add_option(
  parser, c("--max_cons_fold"), action="store", type="integer",
  dest="max_cons_fold", default=79L, metavar="N",
  help="Maximum number of sequences to use in consensus folding [79]"
)
parser <- optparse::add_option(
  parser, c("--force"), action="store_true", type="logical",
  dest="force", default=FALSE, help="Force overwrite of output directory"
)

argv <- suppressWarnings(
        optparse::parse_args(parser, positional_arguments=TRUE)
)
if (length(argv$args) != 3L) {
        message("Incorrect number of arguments")
        optparse::print_help(parser)
        q(save="no", status=64)
}

library(nocornacod)
completed <- predict_ncrna(
  argv$args[[1]], argv$args[[2]], argv$args[[3]], argv$options$inner_pad,
  argv$options$outer_pad, argv$options$min_ncrna_len, argv$options$min_score,
  argv$options$smooth_win, argv$options$max_cons_fold, argv$options$threads,
  argv$options$force
)
if (!completed) {
        q(save="no", status=65)
}
