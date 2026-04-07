## Tests for format_region()

## Build a minimal region with all columns required by all methods
make_format_region <- function(n = 5L) {
    set.seed(7L)
    snps <- paste0("rs", seq_len(n))
    dt <- data.table::data.table(
        snp        = snps,
        dbsnpid    = as.integer(seq_len(n)),
        pos        = as.integer(seq(1e6, by = 1000L, length.out = n)),
        allele     = rep("A", n),
        oth.allele = rep("G", n),
        beta       = rnorm(n, sd = 0.2),
        se         = rep(0.05, n),
        varbeta    = rep(0.0025, n),
        pvalue     = runif(n, 1e-6, 0.1),
        n_samples  = 5000L,
        freq       = runif(n, 0.05, 0.45)
    )
    ldmat <- diag(n)
    rownames(ldmat) <- snps
    colnames(ldmat) <- snps
    list(region = dt, ldmat = ldmat)
}

## ---- ABF ---------------------------------------------------------------

test_that("format_region ABF returns a named list with correct elements", {
    d   <- make_format_region()
    out <- format_region(d$region, "abf", d$ldmat)

    expect_type(out, "list")
    ## LD is appended after the core fields when ldmat is supplied
    expect_named(out, c("snp", "position", "beta", "varbeta",
                        "pvalues", "MAF", "N", "type", "LD"))
})

test_that("format_region ABF values are correctly mapped from region columns", {
    d   <- make_format_region()
    out <- format_region(d$region, "abf", d$ldmat)

    expect_equal(out$snp,      d$region$snp)
    expect_equal(out$position, d$region$pos)
    expect_equal(out$beta,     d$region$beta)
    expect_equal(out$varbeta,  d$region$varbeta)
    expect_equal(out$pvalues,  d$region$pvalue)
    expect_equal(out$MAF,      d$region$freq)
    expect_equal(out$N,        5000L)
    expect_equal(out$type,     "quant")
})

test_that("format_region ABF passes coloc::check_dataset without error", {
    d <- make_format_region()
    expect_no_error(format_region(d$region, "abf", d$ldmat))
})

test_that("format_region ABF works with NULL ldmat (LD omitted from list)", {
    d   <- make_format_region()
    ## When ldmat is NULL, LD should not appear in the output list at all
    ## (coloc::check_dataset errors if LD=NULL is present)
    out <- suppressMessages(format_region(d$region, "abf", NULL))
    expect_false("LD" %in% names(out))
})

## ---- SuSiE -------------------------------------------------------------

test_that("format_region SuSiE returns same structure as ABF", {
    d    <- make_format_region()
    abf  <- format_region(d$region, "abf",   d$ldmat)
    sus  <- format_region(d$region, "susie",  d$ldmat)

    expect_named(sus, names(abf))
    expect_equal(sus$snp,     abf$snp)
    expect_equal(sus$beta,    abf$beta)
    expect_equal(sus$varbeta, abf$varbeta)
    expect_equal(sus$LD,      abf$LD)
})

## ---- PWCoCo ------------------------------------------------------------

test_that("format_region PWCoCo returns a data.table with correct columns", {
    d   <- make_format_region()
    out <- format_region(d$region, "pwcoco", d$ldmat)

    expect_true(data.table::is.data.table(out))
    expect_named(out, c("snp", "eff.allele", "oth.allele",
                        "freq", "beta", "se", "pvalue", "n_samples"))
})

test_that("format_region PWCoCo maps allele column to eff.allele", {
    d   <- make_format_region()
    out <- format_region(d$region, "pwcoco", d$ldmat)

    expect_equal(out$eff.allele, d$region$allele)
    expect_equal(out$oth.allele, d$region$oth.allele)
    expect_equal(out$freq,       d$region$freq)
})

## ---- Error handling ----------------------------------------------------

test_that("format_region errors on unknown analysis_type", {
    d <- make_format_region()
    expect_error(format_region(d$region, "unknown_method", d$ldmat),
                 regexp = "Unknown analysis_type")
})

test_that("format_region errors when required columns are missing", {
    d <- make_format_region()
    d$region[, beta := NULL]
    expect_error(format_region(d$region, "abf", d$ldmat),
                 regexp = "beta")
})
