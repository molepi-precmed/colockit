## Tests for prepare_regions() dispatcher

test_that("prepare_regions errors on unknown backend", {
    expect_error(
        prepare_regions(
            locus        = data.table::data.table(chrom=1L, startpos=1L, endpos=2L),
            gwas1        = data.table::data.table(),
            gwas2        = data.table::data.table(),
            output_dir   = tempdir(),
            refplinkfile = "dummy",
            backend      = "unknown_backend"
        ),
        regexp = "should be one of"
    )
})

test_that("prepare_regions errors when con is NULL for genoscores backend", {
    expect_error(
        prepare_regions(
            locus        = data.table::data.table(chrom=1L, startpos=1L, endpos=2L),
            gwas1        = 1L,
            gwas2        = 2L,
            output_dir   = tempdir(),
            refplinkfile = "dummy",
            backend      = "genoscores",
            con          = NULL
        ),
        regexp = "con"
    )
})

test_that(".check_locus errors when columns missing", {
    expect_error(.check_locus(data.table::data.table(chrom = 1L)),
                 regexp = "startpos")
})

test_that(".valid_dir creates directory if absent", {
    tmp <- file.path(tempdir(), paste0("colockit_test_", Sys.time()))
    expect_false(dir.exists(tmp))
    out <- .valid_dir(tmp)
    expect_true(dir.exists(out))
    unlink(tmp, recursive = TRUE)
})

test_that(".check_refplink errors when .bim file is absent", {
    expect_error(
        colockit:::.check_refplink("/nonexistent/path/panel"),
        regexp = "not found")
})
