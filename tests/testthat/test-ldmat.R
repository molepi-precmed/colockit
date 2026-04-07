## Tests for get_ldmat()
##
## Three tiers:
##   Tier 1 — input validation only, no reference panel needed (always run)
##   Tier 2 — structural + mathematical properties (requires reference panel)
##   Tier 3 — allele-sign correction correctness (requires reference panel)
##
## Test SNPs confirmed polymorphic in kg.2020.hg38.eur (EUR 1000G):
##
##   rs997567   dbsnpid=997567   chr22 pos=29002010  a1=T  a2=C  MAF=0.34
##   rs7289239  dbsnpid=7289239  chr22 pos=29003256  a1=C  a2=T  MAF=0.45
##   rs8142788  dbsnpid=8142788  chr22 pos=29004527  a1=A  a2=G  MAF=0.17
##
## Sign-flip test uses rs7289239: allele=C matches a1 (no flip),
##                                allele=T matches a2 (flip needed).

refplink_test  <- "/Users/andrico/Sites/tests-data/kg.2020.hg38.eur"
panel_available <- file.exists(paste0(refplink_test, ".bim"))

## ======================================================================
##  Tier 1 — Input validation (no reference panel needed)
## ======================================================================

test_that("get_ldmat errors when required columns are missing", {
    region <- data.table::data.table(
        dbsnpid = c(997567L, 7289239L),
        allele  = c("T", "C")
    )
    expect_error(get_ldmat(region, refplink_test), regexp = "snp")
})

test_that("get_ldmat errors when snp column contains non-rs identifiers", {
    region <- data.table::data.table(
        snp     = c("997567", "7289239"),
        dbsnpid = c(997567L,  7289239L),
        allele  = c("T", "C")
    )
    expect_error(get_ldmat(region, refplink_test), regexp = "rs-identifier")
})

test_that("get_ldmat errors when dbsnpid is not integer", {
    region <- data.table::data.table(
        snp     = c("rs997567", "rs7289239"),
        dbsnpid = c(997567,     7289239),    # numeric, not integer
        allele  = c("T", "C")
    )
    expect_error(get_ldmat(region, refplink_test), regexp = "integer")
})

test_that("get_ldmat errors when dbsnpid contains duplicates", {
    region <- data.table::data.table(
        snp     = c("rs997567", "rs997567"),
        dbsnpid = c(997567L,    997567L),
        allele  = c("T", "T")
    )
    expect_error(get_ldmat(region, refplink_test), regexp = "[Dd]uplicat")
})

test_that("get_ldmat errors when reference .bim file does not exist", {
    region <- data.table::data.table(
        snp     = "rs997567",
        dbsnpid = 997567L,
        allele  = "T"
    )
    expect_error(get_ldmat(region, "/nonexistent/path/panel"),
                 regexp = "\\.bim")
})

## ======================================================================
##  Tier 2 — Structural and mathematical properties
## ======================================================================

test_that("get_ldmat returns NULL when no SNPs match the reference panel", {
    skip_if_not(panel_available, "Reference panel not available")

    region <- data.table::data.table(
        snp     = c("rs999999991", "rs999999992"),
        dbsnpid = c(999999991L,   999999992L),
        allele  = c("A", "C")
    )
    expect_null(suppressMessages(get_ldmat(region, refplink_test)))
})

test_that("get_ldmat returns a named list with region and ldmat", {
    skip_if_not(panel_available, "Reference panel not available")

    region <- data.table::data.table(
        snp     = c("rs997567", "rs7289239", "rs8142788"),
        dbsnpid = c(997567L,     7289239L,    8142788L),
        allele  = c("T",         "C",         "A")
    )
    ld <- get_ldmat(region, refplink_test)

    expect_type(ld, "list")
    expect_named(ld, c("region", "ldmat"))
    expect_true(inherits(ld$ldmat, "matrix"))
})

test_that("get_ldmat matrix dimensions and names match returned region", {
    skip_if_not(panel_available, "Reference panel not available")

    region <- data.table::data.table(
        snp     = c("rs997567", "rs7289239", "rs8142788"),
        dbsnpid = c(997567L,     7289239L,    8142788L),
        allele  = c("T",         "C",         "A")
    )
    ld <- get_ldmat(region, refplink_test)
    m  <- ld$ldmat
    n  <- nrow(ld$region)

    expect_equal(dim(m), c(n, n))
    expect_equal(rownames(m), ld$region$snp)
    expect_equal(colnames(m), ld$region$snp)
})

test_that("get_ldmat matrix is symmetric with diagonal == 1", {
    skip_if_not(panel_available, "Reference panel not available")

    region <- data.table::data.table(
        snp     = c("rs997567", "rs7289239", "rs8142788"),
        dbsnpid = c(997567L,     7289239L,    8142788L),
        allele  = c("T",         "C",         "A")
    )
    m <- get_ldmat(region, refplink_test)$ldmat

    expect_equal(m, t(m),                        tolerance = 1e-10)
    expect_equal(unname(diag(m)), rep(1, nrow(m)), tolerance = 1e-10)
})

test_that("get_ldmat off-diagonal values are in [-1, 1]", {
    skip_if_not(panel_available, "Reference panel not available")

    region <- data.table::data.table(
        snp     = c("rs997567", "rs7289239", "rs8142788"),
        dbsnpid = c(997567L,     7289239L,    8142788L),
        allele  = c("T",         "C",         "A")
    )
    m       <- get_ldmat(region, refplink_test)$ldmat
    offdiag <- m[row(m) != col(m)]

    expect_true(all(offdiag >= -1 - 1e-10))
    expect_true(all(offdiag <=  1 + 1e-10))
})

## ======================================================================
##  Tier 3 — Allele-sign correction
##
##  rs7289239 has a1=C, a2=T in the reference panel.
##
##  Run 1: allele=C (matches a1, no flip needed)
##  Run 2: allele=T (matches a2, flip needed)
##
##  Expected:
##    LD(rs997567, rs7289239) in Run 1 == -LD(rs997567, rs7289239) in Run 2
##    Diagonals remain 1 in both runs
##    Matrices remain symmetric in both runs
## ======================================================================

test_that("get_ldmat flips sign correctly when allele matches a2 in reference", {
    skip_if_not(panel_available, "Reference panel not available")

    ## Run 1: rs7289239 allele = C (= a1, no flip)
    region_noflip <- data.table::data.table(
        snp     = c("rs997567", "rs7289239"),
        dbsnpid = c(997567L,    7289239L),
        allele  = c("T",        "C")
    )

    ## Run 2: rs7289239 allele = T (= a2, flip needed)
    region_flip <- data.table::data.table(
        snp     = c("rs997567", "rs7289239"),
        dbsnpid = c(997567L,    7289239L),
        allele  = c("T",        "T")
    )

    ld_noflip <- get_ldmat(region_noflip, refplink_test)
    ld_flip   <- get_ldmat(region_flip,   refplink_test)

    skip_if(is.null(ld_noflip) || is.null(ld_flip),
            "SNPs absent from panel — cannot run sign-flip test")

    r_noflip <- ld_noflip$ldmat["rs997567", "rs7289239"]
    r_flip   <- ld_flip$ldmat[  "rs997567", "rs7289239"]

    expect_equal(r_noflip, -r_flip, tolerance = 1e-10,
                 label = "LD with/without allele flip must have opposite signs")

    expect_equal(unname(diag(ld_noflip$ldmat)), c(1, 1), tolerance = 1e-10)
    expect_equal(unname(diag(ld_flip$ldmat)),   c(1, 1), tolerance = 1e-10)

    expect_equal(ld_noflip$ldmat, t(ld_noflip$ldmat), tolerance = 1e-10)
    expect_equal(ld_flip$ldmat,   t(ld_flip$ldmat),   tolerance = 1e-10)
})
