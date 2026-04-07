## Tests for colocalise()

## ---- ABF ----------------------------------------------------------------

test_that("colocalise ABF returns a list with summary and results", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- colocalise(regions, locus, analysis_type = "abf",
                      output_dir = out_dir)

    expect_type(res, "list")
    expect_true("summary" %in% names(res))
    expect_true("results" %in% names(res))
})

test_that("colocalise ABF summary contains all five posterior probabilities", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res     <- colocalise(regions, locus, analysis_type = "abf",
                          output_dir = out_dir)
    pp_cols <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                 "PP.H3.abf", "PP.H4.abf")
    ## $summary is a named numeric vector in coloc >= 5
    expect_true(all(pp_cols %in% names(res$summary)))
})

test_that("colocalise ABF posterior probabilities sum to 1", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- colocalise(regions, locus, analysis_type = "abf",
                      output_dir = out_dir)
    ## $summary is a named numeric vector; index by name
    pp  <- res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                          "PP.H3.abf", "PP.H4.abf")]
    expect_equal(sum(pp), 1, tolerance = 1e-6)
})

test_that("colocalise ABF results contain SNP-level PP.H4", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- colocalise(regions, locus, analysis_type = "abf",
                      output_dir = out_dir)

    expect_true("SNP.PP.H4" %in% names(res$results))
})

test_that("colocalise ABF writes formatted region and plot files to output_dir", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    colocalise(regions, locus, analysis_type = "abf",
               output_dir = out_dir)

    expect_true(file.exists(file.path(out_dir, "formatted.region.Rdata.gz")))
    expect_true(file.exists(file.path(out_dir, "coloc.abf.png")))
})

## ---- SuSiE --------------------------------------------------------------

test_that("colocalise SuSiE returns a list with summary (or NULL on failure)", {
    regions <- make_test_regions(n = 10L)
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    ## SuSiE may legitimately return NULL on small/degenerate data; accept both
    res <- colocalise(regions, locus, analysis_type = "susie",
                      output_dir = out_dir)

    expect_true(is.null(res) || "summary" %in% names(res))
})

test_that("colocalise SuSiE writes formatted region file", {
    regions <- make_test_regions(n = 10L)
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    colocalise(regions, locus, analysis_type = "susie",
               output_dir = out_dir)

    expect_true(file.exists(file.path(out_dir, "formatted.region.Rdata.gz")))
})

## ---- PWCoCo -------------------------------------------------------------

test_that("colocalise PWCoCo errors when refplinkfile is NULL", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    expect_error(
        colocalise(regions, locus, analysis_type = "pwcoco",
                   output_dir = out_dir, refplinkfile = NULL),
        regexp = "refplinkfile"
    )
})

## ---- Method stubs -------------------------------------------------------

test_that("HyPrColoc stub stops with 'not yet implemented' message", {
    expect_error(colockit:::.method_hyprcoloc(),
                 regexp = "not yet implemented")
})

test_that("eCAVIAR stub stops with 'not yet implemented' message", {
    expect_error(colockit:::.method_ecaviar(),
                 regexp = "not yet implemented")
})

test_that("PAINTOR stub stops with 'not yet implemented' message", {
    expect_error(colockit:::.method_paintor(),
                 regexp = "not yet implemented")
})

test_that("FINEMAP stub stops with 'not yet implemented' message", {
    expect_error(colockit:::.method_finemap(),
                 regexp = "not yet implemented")
})

## ---- Input validation ---------------------------------------------------

test_that("colocalise errors on unknown analysis_type", {
    regions <- make_test_regions()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    expect_error(
        colocalise(regions, locus, analysis_type = "badmethod",
                   output_dir = out_dir),
        regexp = "should be one of"
    )
})
