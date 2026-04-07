## Tests for shared utility functions

test_that("check_missing_cols passes when all columns present", {
    expect_invisible(check_missing_cols(c("a", "b", "c"), c("a", "b")))
})

test_that("check_missing_cols errors with informative message", {
    expect_error(
        check_missing_cols(c("a", "b"), c("a", "b", "c")),
        regexp = "c"
    )
})

test_that("check_region removes beta == 0", {
    dt <- data.table::data.table(
        dbsnpid = 1:3,
        beta    = c(0.1, 0, -0.2),
        pvalue  = c(0.01, 0.05, 0.001)
    )
    out <- check_region(dt)
    expect_equal(nrow(out), 2L)
    expect_true(all(out$beta != 0))
})

test_that("check_region removes pvalue >= 1", {
    dt <- data.table::data.table(
        dbsnpid = 1:3,
        beta    = c(0.1, 0.2, 0.3),
        pvalue  = c(0.01, 1.0, 0.5)
    )
    out <- check_region(dt)
    expect_equal(nrow(out), 2L)
})

test_that("check_region clamps pvalue == 0", {
    dt <- data.table::data.table(
        dbsnpid = 1:2,
        beta    = c(0.1, 0.2),
        pvalue  = c(0, 0.05)
    )
    out <- check_region(dt)
    expect_equal(nrow(out), 2L)
    expect_equal(out$pvalue[1L], .Machine$double.xmin)
})

test_that("check_region removes bad MAF when maf column present", {
    dt <- data.table::data.table(
        dbsnpid = 1:3,
        beta    = c(0.1, 0.2, 0.3),
        pvalue  = c(0.01, 0.02, 0.03),
        maf     = c(0.1, 0, 1.0)
    )
    out <- check_region(dt)
    expect_equal(nrow(out), 1L)
})

test_that("get_se_varbeta computes varbeta from existing se", {
    dt <- data.table::data.table(
        beta      = c(0.5, -0.3),
        se        = c(0.1,  0.05),
        pvalue    = c(0.01, 0.05),
        n_samples = c(1000L, 1000L)
    )
    out <- get_se_varbeta(dt)
    expect_true("varbeta" %in% names(out))
    expect_equal(out$varbeta, out$se^2)
})

test_that("get_se_varbeta derives se when absent", {
    dt <- data.table::data.table(
        beta      = 0.5,
        pvalue    = 1e-6,
        n_samples = 5000L
    )
    out <- get_se_varbeta(dt)
    expect_true("se"      %in% names(out))
    expect_true("varbeta" %in% names(out))
    expect_true(out$se > 0)
})
