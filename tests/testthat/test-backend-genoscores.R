## Tests for genoscores backend helpers
##
## Tier 0 — no external dependencies: pure functions testable with synthetic data.
## Tier DB — require genoscores package + live DB connection (skip otherwise).
##
## The live-DB tests assume a connection object named `scorecon` exists in the
## global environment (the conventional variable used by genoscores scripts).
## If it does not exist the tests are skipped gracefully.

has_genoscores <- requireNamespace("genoscores", quietly = TRUE)
has_scorecon   <- exists("scorecon", envir = .GlobalEnv, inherits = FALSE)

## ---- compare_freqs (no DB required) -------------------------------------

test_that("compare_freqs removes SNPs where |freq - eff.freq| > tolerance", {
    region <- data.table::data.table(
        snp      = c("rs1", "rs2", "rs3"),
        dbsnpid  = 1:3L,
        freq     = c(0.30, 0.30, 0.30),   # dbSNP freq
        eff.freq = c(0.30, 0.45, 0.10),   # SNPTEST freq  (rs2 and rs3 deviate)
        oth.freq = c(0.70, 0.55, 0.90)
    )
    out <- colockit:::compare_freqs(region, tolerance = 0.1)

    ## rs2 differs by 0.15, rs3 by 0.20 — both removed
    expect_equal(nrow(out), 1L)
    expect_equal(out$snp, "rs1")
})

test_that("compare_freqs replaces freq with eff.freq and drops helper columns", {
    region <- data.table::data.table(
        snp      = "rs1",
        dbsnpid  = 1L,
        freq     = 0.28,
        eff.freq = 0.30,
        oth.freq = 0.70
    )
    out <- colockit:::compare_freqs(region, tolerance = 0.1)

    expect_equal(out$freq, 0.30)
    expect_false("eff.freq" %in% names(out))
    expect_false("oth.freq" %in% names(out))
})

test_that("compare_freqs keeps all SNPs when all frequencies agree", {
    region <- data.table::data.table(
        snp      = paste0("rs", 1:5),
        dbsnpid  = 1:5L,
        freq     = rep(0.3, 5L),
        eff.freq = rep(0.3, 5L),
        oth.freq = rep(0.7, 5L)
    )
    out <- colockit:::compare_freqs(region, tolerance = 0.1)
    expect_equal(nrow(out), 5L)
})

## ---- process_private_gwas (no DB required) --------------------------------

make_raw_snptest <- function(n = 5L, seed = 1L) {
    set.seed(seed)
    data.table::data.table(
        snp        = paste0("rs", seq_len(n)),
        dbsnpid    = as.integer(seq_len(n)),
        chrom      = 1L,
        pos        = as.integer(seq(1e6, by = 1000L, length.out = n)),
        oth.allele = rep("A", n),
        allele     = rep("G", n),
        beta       = rnorm(n, sd = 0.1),
        se         = rep(0.05, n),
        pvalue     = runif(n, 1e-5, 0.2),
        n_AA       = 200L,
        n_AB       = 400L,
        n_BB       = 400L,
        maf        = rep(0.4, n),
        n_samples  = 1000L
    )
}

test_that("process_private_gwas errors when non-rsid SNPs are present", {
    gwas <- make_raw_snptest()
    gwas[1L, snp := "chr1:1000000"]   # non-rs identifier
    expect_error(colockit:::process_private_gwas(gwas, "phenotype"),
                 regexp = "rsid")
})

test_that("process_private_gwas computes eff.freq and renames to eaf", {
    gwas <- make_raw_snptest()
    out  <- colockit:::process_private_gwas(gwas, "phenotype")
    expect_true("eaf" %in% names(out))
})

test_that("process_private_gwas drops n_AA, n_AB, n_BB, oth.freq, minor", {
    gwas <- make_raw_snptest()
    out  <- colockit:::process_private_gwas(gwas, "phenotype")
    dropped <- c("n_AA", "n_AB", "n_BB", "oth.freq", "minor")
    expect_false(any(dropped %in% names(out)))
})

test_that("process_private_gwas sets gwasid to 0 and trait.name to outcome", {
    gwas <- make_raw_snptest()
    out  <- colockit:::process_private_gwas(gwas, "my_phenotype")
    expect_equal(unique(out$gwasid), 0L)
    expect_equal(unique(out$trait.name), "my_phenotype")
})

test_that("process_private_gwas drops rows with NA pvalue", {
    gwas <- make_raw_snptest(n = 6L)
    gwas[3L, pvalue := NA_real_]
    out <- colockit:::process_private_gwas(gwas, "phenotype")
    expect_equal(nrow(out), 5L)
})

## ---- read_snptest (no DB required, uses a temp file) --------------------

test_that("read_snptest reads a SNPTEST-like file and renames columns", {
    ## build a minimal file mimicking SNPTEST output (header line prepended)
    header  <- paste(
        "alternate_ids",
        "rsid", "chromosome", "position",
        "alleleA", "alleleB",
        "frequentist_add_beta_1", "frequentist_add_se_1",
        "frequentist_add_pvalue",
        "all_AA", "all_AB", "all_BB", "all_maf", "all_total",
        sep = " ")
    data_row <- "rs1_alt rs1 1 1000000 A G 0.12 0.05 0.02 50 100 50 0.25 200"

    tmp <- tempfile(fileext = ".snptest")
    writeLines(c(header, data_row), tmp)

    out <- colockit:::read_snptest(tmp)
    unlink(tmp)

    expect_equal(names(out), c("snp", "chrom", "pos", "oth.allele", "allele",
                                "beta", "se", "pvalue",
                                "n_AA", "n_AB", "n_BB", "maf", "n_samples"))
    expect_equal(out$snp, "rs1")
    expect_equal(out$beta, 0.12)
})

## ---- DB-dependent tests ---------------------------------------------------

test_that("get_published_gwas returns a data.table with expected columns", {
    skip_if_not(has_genoscores,
                "genoscores package is not available")
    skip_if_not(has_scorecon,
                "scorecon DB connection not found in global env")

    ## gwasid 1 is assumed to exist in the test/development database;
    ## adjust if your instance uses different ids
    con <- get("scorecon", envir = .GlobalEnv)
    expect_error(
        out <- colockit:::get_published_gwas(1L, con),
        NA)  # no error expected
    expect_true(data.table::is.data.table(out))
    expect_true(all(c("snp", "dbsnpid", "gwasid", "trait.name",
                       "n_samples") %in% names(out)))
})

test_that("get_freqs_gs returns freq column in (0, 1)", {
    skip_if_not(has_genoscores,
                "genoscores package is not available")
    skip_if_not(has_scorecon,
                "scorecon DB connection not found in global env")

    con  <- get("scorecon", envir = .GlobalEnv)
    gwas <- data.table::data.table(dbsnpid = 997567L, allele = "T")
    out  <- colockit:::get_freqs_gs(gwas, con)

    expect_true(data.table::is.data.table(out))
    expect_true(all(out$freq > 0 & out$freq < 1))
})

test_that("find_missing_alleles_gs resolves alleles from dbSNP", {
    skip_if_not(has_genoscores,
                "genoscores package is not available")
    skip_if_not(has_scorecon,
                "scorecon DB connection not found in global env")

    con     <- get("scorecon", envir = .GlobalEnv)
    missing <- data.table::data.table(
        snp      = "rs997567",
        dbsnpid  = 997567L,
        allele   = "T",
        chrom    = 22L,
        pos      = 29000000L
    )
    out <- colockit:::find_missing_alleles_gs(missing, con)

    expect_true("oth.allele" %in% names(out))
    expect_equal(nrow(out), 1L)
    expect_false(is.na(out$oth.allele))
})

## ---- Top-level dispatch test ---------------------------------------------

test_that("prepare_regions dispatches to genoscores backend and returns list", {
    skip_if_not(has_genoscores,
                "genoscores package is not available")
    skip_if_not(has_scorecon,
                "scorecon DB connection not found in global env")
    skip_if_not(
        file.exists("/Users/andrico/Sites/tests-data/kg.2020.hg38.eur.bim"),
        "Reference PLINK panel not available")

    con   <- get("scorecon", envir = .GlobalEnv)
    locus <- make_test_locus(chrom = 22L, startpos = 28900000L,
                              endpos = 29100000L)

    ## gwasid 1 and 2 must exist in the development database
    res <- prepare_regions(
        locus        = locus,
        gwas1        = 1L,
        gwas2        = 2L,
        output_dir   = tempdir(),
        refplinkfile = "/Users/andrico/Sites/tests-data/kg.2020.hg38.eur",
        backend      = "genoscores",
        gwas_type    = "published",
        con          = con)

    ## result may be NULL if no common SNPs — that is also acceptable
    expect_true(is.null(res) || is.list(res))
    if (!is.null(res)) {
        expect_true(all(c("src", "trgt", "ldmat") %in% names(res)))
    }
})
