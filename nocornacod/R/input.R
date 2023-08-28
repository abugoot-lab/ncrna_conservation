MIN_CLUSTER_SIZE <- 2L
# There needs to be enough bases outside of the CDS for the ncRNA prediction to
# work, otherwise we would just be predicting conservation of the CDS.
MIN_CONTEXT_LEN <- 500L

read_user_input <- function(seqs_path, coords_path) {
    seqs <- read_seqs(seqs_path)
    coords <- read_coords(coords_path)
    tokeep <- intersect(coords$contig, names(seqs))

    assert(
        length(tokeep) > 0,
        msg="No shared IDs between sequences and coordinates"
    )

    coords <- coords[coords$contig %in% tokeep, ]
    coords$seqlength <- Biostrings::width(seqs[coords$contig])
    coords <- coords[
        coords$start > MIN_CONTEXT_LEN &
        coords$seqlength - coords$end > MIN_CONTEXT_LEN,
    ]

    assert(
        nrow(coords) >= MIN_CLUSTER_SIZE,
        msg=paste0(
            "Only ", nrow(coords), " sequences after filtering. Need ",
            MIN_CLUSTER_SIZE
        )
    )

    # This is important because some other functions assume the coordinates and
    # sequences are in the same order.
    seqs <- seqs[coords$contig]
    input <- list(coordinates=coords, sequences=seqs) |>
        convert2plus()

    return(input)
}

read_seqs <- function(path) {
    seqs <- tryCatch(
        error=function(cnd) stop("Failed to read sequence file", call.=FALSE),
        Biostrings::readDNAStringSet(path)
    )

    ids <- trimws(names(seqs), which="both") |>
        strsplit(split="[[:space:]]") |>
        vapply(FUN=`[[`, FUN.VALUE=character(1L), 1L)
    dups <- duplicated(ids)
    if (any(dups)) {
        Map(
            \(x) warning(x, call.=FALSE),
            sprintf("Ignoring sequence with duplicate ID %s", ids[dups])
        )
    }
    seqs <- seqs[!dups]
    names(seqs) <- ids[!dups]

    useqs <- get_unique_seqs(seqs)
    message(msgf(sprintf(
        "%d / %d unique sequences", length(useqs), length(seqs)
    )))

    return(useqs)
}

read_coords <- function(path) {
    tbl <- tryCatch(
        error=function(cnd) {
            stop("Failed to parse coordinates", call.=FALSE)
        },
        utils::read.table(
            path, col.names=c("contig", "start", "end", "strand"), sep="\t",
            colClasses=c("character", "integer", "integer", "character")
        )
    )

    ids <- trimws(tbl$contig, which="both") |>
        strsplit(split="[[:space:]]") |>
        vapply(FUN=`[[`, FUN.VALUE=character(1L), 1L)

    dups <- duplicated(ids)
    if (any(dups)) {
        Map(
            \(x) warning(x, call.=FALSE),
            sprintf("Ignoring coordiantes with duplicate contig %s", ids[dups])
        )
    }
    tbl <- tbl[!dups, ]
    tbl$contig <- ids[!dups]

    return(tbl)
}

convert2plus <- function(x) {
    co <- x$coordinates
    s <- x$sequences
    if (all(co$strand == "+")) {
        return(x)
    }
    w <- co$end - co$start
    strd <- co$strand
    co$start <- ifelse(co$strand == "+", co$start, co$seqlength - co$end + 1)
    co$end <- co$start + w 
    co$strand <- "+"
    s[strd == "-"] <- Biostrings::reverseComplement(s[strd == "-"])

    x$coordinates <- co
    x$sequences <- s

    return(x)
}

get_unique_seqs <- function(x) {
    hashes <- vapply(
        methods::as(x, "character"),
        digest::digest, character(1),
        algo="xxhash64", serialize=FALSE
    )
    uniq <- x[names(hashes)[!duplicated(hashes)]]

    return(uniq)
}
