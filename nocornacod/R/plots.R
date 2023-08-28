R2R_PAD <- 5

make_plots <- function(x, outdir_path) {
    co <- x$coordinates
    cons_plot5 <- make_conservation_plot(x$msa5_scores, x$cpt5)
    gene_plot5 <- make_genes_plot(
        x$msa5_scores, x$cpt5, co$contig, co$context5_mr,
        length(x$msa5_scores), "-"
    )
    cons_plot3 <- make_conservation_plot(x$msa3_scores, x$cpt3)
    gene_plot3 <- make_genes_plot(
        x$msa3_scores, x$cpt3, co$contig, 1, co$context3_mr, "+"
    )

    cpt_temp <- tempfile("conservation", outdir_path, ".svg")
    on.exit(cleanup(cpt_temp), add=TRUE)

    p <- patchwork::wrap_plots(
        cons_plot5, gene_plot5, cons_plot3, gene_plot3, ncol=1, nrow=4
    )
    svglite::svglite(cpt_temp)
    print(p)
    invisible(grDevices::dev.off())
    output <- file.path(outdir_path, "conservation.svg")
    five_prime <- NULL
    three_prime <- NULL
    if (!is.null(x$fold5)) {
        five_prime <- list(
            svg=file.path(
                x$fold5$rscape, "results.R2R.sto.svg"
            ),
            energy=x$fold5$energies[["total"]]
        )
    }

    if (!is.null(x$fold3)) {
        three_prime <- list(
            svg=file.path(
                x$fold3$rscape, "results.R2R.sto.svg"
            ),
            energy=x$fold3$energies[["total"]]
        )
    }
    compose_plots(output, cpt_temp, five_prime, three_prime)
}

make_gene_tbl <- function(contig, start, end, strand, max_genes) {
    tbl <- data.frame(
        molecule=contig,
        gene=paste("gene", seq_along(contig), sep=""),
        start=start,
        end=end,
        orientation=ifelse(strand == "+", TRUE, FALSE)
    )

    if (nrow(tbl) > max_genes) {
        tbl <- tbl[sample(seq_len(nrow(tbl)), max_genes), ]
    }

    return(tbl)
}

make_conservation_plot <- function(scores, cpts) {
    column <- NULL
    score <- NULL

    scores_tbl <- data.frame(column=seq_along(scores), score=scores)

    p <- ggplot2::ggplot(scores_tbl, ggplot2::aes(x=column, y=score)) +
        ggplot2::geom_line()

    if (!is.null(cpts)) {
        p <- p + ggplot2::geom_vline(xintercept=cpts, colour="red")
    }
    p <- p +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(y="Conservation Score") +
        ggplot2::theme(
            axis.text.x=ggplot2::element_blank(),
            axis.title.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank(),
            axis.title=ggplot2::element_text(size=8)
        )

    return(p)
}

make_genes_plot <- function(
    scores, cpts, contig, start, end, strand, max_genes=9L
) {
    molecule <- NULL
    gene <- NULL
    orientation <- NULL

    gene_tbl <- make_gene_tbl(contig, start, end, strand, max_genes)
    p <- ggplot2::ggplot(gene_tbl,
        ggplot2::aes(
            xmin=start, xmax=end, y=molecule, label=gene, forward=orientation
        )
    ) +
      gggenes::geom_gene_arrow(
        arrow_body_height=grid::unit(4, "mm"),
        arrowhead_height=grid::unit(4, "mm"),
        arrowhead_width=grid::unit(1, "mm"),
        fill="coral"
      ) +
      gggenes::geom_gene_label()

    if (!is.null(cpts)) {
        p <- p + ggplot2::geom_vline(xintercept=cpts, colour="red")
    }
    p <- p +
        ggplot2::xlim(1, length(scores)) +
        ggplot2::labs(x="Column") +
        gggenes::theme_genes() +
        ggplot2::theme(
            axis.text.y=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank(),
            axis.title=ggplot2::element_text(size=8)
        )

    return(p)
}

compose_plots <- function(output, cpt, fp, tp) {
    cpt <- xml2::read_xml(cpt)
    cpt_width <- xml2::xml_attr(cpt, "width") |>
        (\(x) sub("pt$", "", x))() |>
        as.double()
    cpt_height <- xml2::xml_attr(cpt, "height") |>
        (\(x) sub("pt$", "", x))() |>
        as.double()
    xml2::xml_name(cpt) <- "g"
    xml2::xml_find_first(cpt, "d1:rect[@height='100%' and @width='100%']") |>
        xml2::xml_remove()

    if (!is.null(fp)) {
        r2r1 <- xml2::read_xml(fp$svg)
        r2r1_width <- xml2::xml_attr(r2r1, "width") |>
            as.double() |>
            trunc()
        r2r1_height <- xml2::xml_attr(r2r1, "height") |>
            as.double() |>
            trunc()
        xml2::xml_find_first(r2r1, "/svg:svg/svg:text/svg:tspan") |>
            xml2::xml_set_text(sprintf("5' Fold (%.2f)", fp$energy)) |>
            invisible()
        xml2::xml_name(r2r1) <- "g"
    }

    if (!is.null(tp)) {
        r2r2 <- xml2::read_xml(tp$svg)
        r2r2_width <- xml2::xml_attr(r2r2, "width") |>
            as.double() |>
            trunc()
        r2r2_height <- xml2::xml_attr(r2r2, "height") |>
            as.double() |>
            trunc()
        xml2::xml_find_first(r2r2, "/svg:svg/svg:text/svg:tspan") |>
            xml2::xml_set_text(sprintf("3' Fold (%.2f)", tp$energy)) |>
            invisible()
        xml2::xml_name(r2r2) <- "g"
    }

    if (!is.null(fp) && is.null(tp)) {
        r2r_width <- r2r1_width
        r2r_height <- r2r1_height * 2
        r2r2_width <- r2r1_width
        r2r2_height <- r2r1_height
    } else if (is.null(fp) && !is.null(tp)) {
        r2r_width <- r2r2_width
        r2r_height <- r2r2_height * 2
        r2r1_width <- r2r2_width
        r2r1_height <- r2r2_height
    } else {
        r2r_width <- max(c(r2r1_width, r2r2_width))
        r2r_height <- r2r1_height + r2r2_height
    }

    p_width <- cpt_width + r2r_width
    p_height <- max(c(cpt_height, r2r_height))
    r2r_h_adj <- ifelse(
        p_height > r2r_height,
        (p_height - r2r_height) %/% 2,
        0
    )
    cpt_coords <- c(
        0,
        ifelse(cpt_height < p_height, (p_height - cpt_height) %/% 2, 0)
    )

    if (r2r1_width < r2r2_width) {
        r2r1_coords <- c(cpt_width + (r2r2_width - r2r1_width) %/% 2, r2r_h_adj)
        r2r2_coords <- c(cpt_width, r2r1_height + r2r_h_adj)
    } else {
        r2r1_coords <- c(cpt_width, r2r_h_adj)
        r2r2_coords <- c(
            cpt_width + (r2r1_width - r2r2_width) %/% 2,
            r2r1_height + r2r_h_adj
        )
    }

    root <- xml2::xml_new_root("svg")
    xml2::xml_attr(root, "xmlns") <- "http://www.w3.org/2000/svg"
    xml2::xml_attr(root, "xmlns:xlink") <- "http://www.w3.org/1999/xlink"
    xml2::xml_attr(root, "version") <- "1.0"
    xml2::xml_attr(root, "width") <- p_width
    xml2::xml_attr(root, "height") <- p_height + R2R_PAD
    xml2::xml_attr(root, "viewBox") <- sprintf("0 0 %d %d", p_width, p_height + R2R_PAD)
    rec <- xml2::xml_new_root("rect")
    xml2::xml_attrs(rec) <- c(
        x="0", y="0", width="100%", height="100%", fill="#FFFFFF",
        stroke="#FFFFFF"
    )

    xml2::xml_attrs(cpt) <- c(
        "transform"=sprintf(
            "translate(%d,%d)",
            cpt_coords[1],
            cpt_coords[2] + R2R_PAD
        ),
        "class"="svglite"
    )
    xml2::xml_add_child(root, rec)
    xml2::xml_add_child(root, cpt)

    if (!is.null(fp)) {
        xml2::xml_attrs(r2r1) <- c(
            "transform"=sprintf(
                "translate(%d,%d)",
                r2r1_coords[1],
                r2r1_coords[2] + R2R_PAD
            )
        )
        xml2::xml_add_child(root, r2r1)
    }

    if (!is.null(tp)) {
        xml2::xml_attrs(r2r2) <- c(
            "transform"=sprintf(
                "translate(%d,%d)",
                r2r2_coords[1],
                r2r2_coords[2] + R2R_PAD
            )
        )
        xml2::xml_add_child(root, r2r2)
    }

    xml2::write_xml(root, output)
}
