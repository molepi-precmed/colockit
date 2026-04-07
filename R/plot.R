#' Generate locus plots and colocalisation visualisations
#'
#' Dual-mode plotting function:
#' * **Preview mode** (default, `coloc = NULL`): plots both GWAS across the
#'   locus, highlighting SNPs common to both in red.
#' * **Colocalisation mode** (`coloc` and `ldmat` supplied): colours SNPs by
#'   their posterior probability of H4, labels the top SNP, and appends LD
#'   heatmaps beneath each locus plot.
#'
#' @param src_region A `data.table` of GWAS 1 summary statistics. Required
#'   columns: `snp`, `dbsnpid`, `chrom`, `pos`, `pvalue`, `beta`,
#'   `trait.name`.
#' @param trgt_region A `data.table` of GWAS 2 summary statistics (same
#'   requirements as `src_region`).
#' @param locus A single-row `data.table` with columns `chrom`, `startpos`,
#'   `endpos`.
#' @param plot_file Full path to the output PNG file.
#' @param ldmat Optional LD correlation matrix (required for colocalisation
#'   mode). Default `NULL`.
#' @param coloc Optional `data.table` or `data.frame` with columns `snp` and
#'   `SNP.PP.H4` (per-SNP posterior probability of H4). Default `NULL`.
#'
#' @return A named list with entries `src` and `trgt` containing only the
#'   SNPs common to both GWAS within the locus.
#'   Returns `NULL` if there are no common SNPs.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_color_gradient
#'   theme element_text labs unit ggsave
#' @importFrom ggrepel geom_label_repel
#' @importFrom gridExtra grid.arrange
#' @importFrom LDheatmap LDheatmap LDheatmap.addScatterplot
#' @keywords internal
plot_regions <- function(src_region, trgt_region, locus, plot_file,
                          ldmat = NULL, coloc = NULL) {
    check_missing_cols(names(src_region),  c("snp", "dbsnpid", "chrom",
                                              "pos", "pvalue", "beta",
                                              "trait.name"))
    check_missing_cols(names(trgt_region), c("snp", "dbsnpid", "chrom",
                                              "pos", "pvalue", "beta",
                                              "trait.name"))

    src_region[,  `:=`(pos    = as.numeric(pos),
                        beta   = as.numeric(beta),
                        pvalue = as.numeric(pvalue))]
    trgt_region[, `:=`(pos    = as.numeric(pos),
                        beta   = as.numeric(beta),
                        pvalue = as.numeric(pvalue))]

    ## ---- helpers -------------------------------------------------------
    .title    <- function(name)
        sprintf("Locus plot for %s", gsub("_", " ", name))

    .subtitle <- function(chrom, xmin, xmax)
        sprintf("Chromosome %s, %.2f-%.2f Mb",
                chrom, xmin / 1e6, xmax / 1e6)

    .base_plot <- function(region, chrom, xmin, xmax) {
        ggplot2::ggplot(region) +
            ggplot2::geom_point(ggplot2::aes(x = pos, y = -log10(pvalue)),
                                shape = 1) +
            ggplot2::scale_x_continuous(
                labels = function(x) paste(round(x / 1e6, 2), "Mb"),
                limits = c(xmin, xmax)) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5)) +
            ggplot2::labs(
                x        = "Position",
                y        = bquote(-log[10](italic(P))),
                title    = .title(region[1, trait.name]),
                subtitle = .subtitle(chrom, xmin, xmax))
    }

    .color_h4 <- function(p, region, coloc) {
        if (!"data.table" %in% class(coloc)) data.table::setDT(coloc)
        dt       <- region[coloc, on = "snp", nomatch = NULL]
        top_snps <- dt[SNP.PP.H4 == max(SNP.PP.H4), snp]
        p <- p +
            ggplot2::geom_point(
                data  = dt,
                ggplot2::aes(x = pos, y = -log10(pvalue),
                             color = SNP.PP.H4),
                size = 2, alpha = 0.5) +
            ggplot2::scale_color_gradient(low = "darkblue", high = "red") +
            ggplot2::geom_point(
                data  = region[snp %in% top_snps],
                ggplot2::aes(x = pos, y = -log10(pvalue)),
                shape = 23, color = "red", fill = "red", size = 2) +
            ggrepel::geom_label_repel(
                data              = region[snp %in% top_snps],
                ggplot2::aes(x = pos, y = -log10(pvalue), label = snp),
                size              = 2,
                force             = 4,
                box.padding       = ggplot2::unit(1.2, "lines"),
                point.padding     = ggplot2::unit(0.3, "lines"),
                min.segment.length = ggplot2::unit(0.1, "lines"),
                max.overlaps      = 30)
        p
    }

    .ld_map <- function(ldmat, pos, pvalue) {
        map <- LDheatmap::LDheatmap(
            ldmat, genetic.distances = pos,
            add.map = TRUE, title = NULL,
            color = grDevices::heat.colors(20), flip = TRUE)
        LDheatmap::LDheatmap.addScatterplot(
            map, -log10(pvalue),
            height = 0.4,
            ylab   = bquote(-log[10](italic(p))))
    }

    ## ---- common SNPs ---------------------------------------------------
    common_snps <- intersect(src_region$dbsnpid, trgt_region$dbsnpid)
    if (length(common_snps) == 0) {
        cat("No common SNPs detected in the region.\n")
        return(NULL)
    }
    src_common  <- src_region[dbsnpid  %in% common_snps]
    trgt_common <- trgt_region[dbsnpid %in% common_snps]

    ## ---- preview mode --------------------------------------------------
    if (is.null(coloc)) {
        xmin <- min(src_region$pos,  trgt_region$pos)
        xmax <- max(src_region$pos, trgt_region$pos)

        p1 <- .base_plot(src_region,  locus$chrom, xmin, xmax) +
            ggplot2::geom_point(data = src_common,
                                ggplot2::aes(x = pos, y = -log10(pvalue)),
                                color = "red")
        p2 <- .base_plot(trgt_region, locus$chrom, xmin, xmax) +
            ggplot2::geom_point(data = trgt_common,
                                ggplot2::aes(x = pos, y = -log10(pvalue)),
                                color = "red")

        plot <- gridExtra::grid.arrange(p1, p2, nrow = 2)
        ggplot2::ggsave(plot_file, plot = plot,
                        width = 2600, height = 1800, units = "px", dpi = 300)
        cat("Preview locus plot saved to:", plot_file, "\n")
    }

    ## ---- colocalisation mode -------------------------------------------
    if (!is.null(coloc) && !is.null(ldmat)) {
        xmin <- min(src_common$pos,  trgt_common$pos)
        xmax <- max(src_common$pos, trgt_common$pos)

        p1 <- .base_plot(src_common,  locus$chrom, xmin, xmax)
        p1 <- .color_h4(p1, src_common, coloc)

        p2 <- .base_plot(trgt_common, locus$chrom, xmin, xmax)
        p2 <- .color_h4(p2, trgt_common, coloc)

        map1 <- .ld_map(ldmat, src_common$pos,  src_common$pvalue)
        map2 <- .ld_map(ldmat, trgt_common$pos, trgt_common$pvalue)

        plot <- gridExtra::grid.arrange(
            p1, map1$LDheatmapGrob,
            p2, map2$LDheatmapGrob,
            layout_matrix = rbind(c(1, 1, 2), c(3, 3, 4)))
        ggplot2::ggsave(plot_file, plot = plot,
                        width = 3000, height = 2400, units = "px", dpi = 300)
        cat("Colocalisation plot saved to:", plot_file, "\n")
    }

    list(src = src_common, trgt = trgt_common)
}
