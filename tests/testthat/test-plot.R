## Tests for plot_regions()

## ---- helpers -------------------------------------------------------------

## Two GWAS with n_common SNPs in common plus n_extra unique to each
make_plot_pair <- function(n_common = 10L, n_extra = 3L) {
    set.seed(5L)

    ## common SNPs share the same snp name AND dbsnpid in both GWAS
    common_snps   <- paste0("rs", seq_len(n_common))
    common_ids    <- seq_len(n_common)

    ## extra SNPs are unique to each GWAS — use widely spaced IDs to avoid overlap
    extra_src_snps  <- paste0("rs", 1000L + seq_len(n_extra))
    extra_src_ids   <- 1000L + seq_len(n_extra)
    extra_trgt_snps <- paste0("rs", 2000L + seq_len(n_extra))
    extra_trgt_ids  <- 2000L + seq_len(n_extra)

    make_gwas <- function(snps, ids, trait) {
        n <- length(snps)
        data.table::data.table(
            snp        = snps,
            dbsnpid    = as.integer(ids),
            chrom      = 1L,
            pos        = as.integer(seq(1e6, by = 1000L, length.out = n)),
            allele     = rep("A", n),
            oth.allele = rep("G", n),
            beta       = rnorm(n, sd = 0.1),
            pvalue     = runif(n, 1e-5, 0.5),
            n_samples  = 1000L,
            freq       = runif(n, 0.05, 0.45),
            trait.name = trait
        )
    }

    src  <- make_gwas(c(common_snps, extra_src_snps),
                      c(common_ids,  extra_src_ids),  "src_trait")
    trgt <- make_gwas(c(common_snps, extra_trgt_snps),
                      c(common_ids,  extra_trgt_ids), "trgt_trait")
    locus <- make_test_locus(endpos = 1100000L)

    list(src = src, trgt = trgt, locus = locus)
}

## ---- preview mode -------------------------------------------------------

test_that("plot_regions preview mode returns a list with src and trgt", {
    d       <- make_plot_pair()
    out_dir <- tempfile(); dir.create(out_dir)
    plot_f  <- file.path(out_dir, "preview.png")

    res <- plot_regions(d$src, d$trgt, d$locus, plot_f)

    expect_type(res, "list")
    expect_named(res, c("src", "trgt"))
})

test_that("plot_regions preview mode returns only common SNPs", {
    d       <- make_plot_pair(n_common = 10L, n_extra = 3L)
    out_dir <- tempfile(); dir.create(out_dir)
    plot_f  <- file.path(out_dir, "preview.png")

    res <- plot_regions(d$src, d$trgt, d$locus, plot_f)

    ## Only the 10 common dbsnpids should be in both returned datasets
    expect_equal(nrow(res$src),  10L)
    expect_equal(nrow(res$trgt), 10L)
    expect_identical(sort(res$src$dbsnpid), sort(res$trgt$dbsnpid))
})

test_that("plot_regions preview mode writes a PNG file", {
    d       <- make_plot_pair()
    out_dir <- tempfile(); dir.create(out_dir)
    plot_f  <- file.path(out_dir, "preview.png")

    plot_regions(d$src, d$trgt, d$locus, plot_f)

    expect_true(file.exists(plot_f))
    expect_gt(file.info(plot_f)$size, 0L)
})

test_that("plot_regions returns NULL when there are no common SNPs", {
    set.seed(9L)
    make_disjoint <- function(trait, id_offset) {
        data.table::data.table(
            snp        = paste0("rs", id_offset + 1:5),
            dbsnpid    = as.integer(id_offset + 1:5),
            chrom      = 1L,
            pos        = as.integer(seq(1e6, by = 1000L, length.out = 5L)),
            allele     = rep("A", 5L),
            oth.allele = rep("G", 5L),
            beta       = rnorm(5L, sd = 0.1),
            pvalue     = runif(5L, 1e-5, 0.5),
            n_samples  = 1000L,
            freq       = runif(5L, 0.1, 0.4),
            trait.name = trait
        )
    }
    src   <- make_disjoint("src_trait",  0L)
    trgt  <- make_disjoint("trgt_trait", 100L)
    locus <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- plot_regions(src, trgt, locus,
                        file.path(out_dir, "no_common.png"))
    expect_null(res)
})

## ---- colocalisation mode ------------------------------------------------

test_that("plot_regions coloc mode writes a PNG file when ldmat and coloc supplied", {
    d       <- make_plot_pair(n_common = 5L, n_extra = 0L)
    out_dir <- tempfile(); dir.create(out_dir)
    plot_f  <- file.path(out_dir, "coloc.png")

    ## minimal coloc results table (snp + SNP.PP.H4)
    common_snps <- d$src$snp
    coloc_res   <- data.table::data.table(
        snp       = common_snps,
        SNP.PP.H4 = runif(length(common_snps), 0, 1)
    )

    n     <- length(common_snps)
    ldmat <- diag(n)
    rownames(ldmat) <- common_snps
    colnames(ldmat) <- common_snps

    plot_regions(d$src, d$trgt, d$locus, plot_f,
                 ldmat = ldmat, coloc = coloc_res)

    expect_true(file.exists(plot_f))
    expect_gt(file.info(plot_f)$size, 0L)
})

## ---- input validation ---------------------------------------------------

test_that("plot_regions errors when required columns are missing from src", {
    d       <- make_plot_pair()
    d$src[, pvalue := NULL]
    out_dir <- tempfile(); dir.create(out_dir)

    expect_error(
        plot_regions(d$src, d$trgt, d$locus,
                     file.path(out_dir, "plot.png")),
        regexp = "pvalue"
    )
})
