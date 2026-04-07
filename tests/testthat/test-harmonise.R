## Tests for harmonise_gwas()

test_that("harmonise_gwas retains SNPs with matching alleles and leaves beta unchanged", {
    pair    <- make_harmonise_pair()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    expect_true(all(c("rs1", "rs2", "rs3") %in% res$src$snp))
    expect_equal(res$src[snp == "rs1", beta], 0.1)
    expect_equal(res$src[snp == "rs2", beta], 0.2)
    expect_equal(res$src[snp == "rs3", beta], 0.3)
})

test_that("harmonise_gwas flips beta and swaps allele labels when alleles are reversed", {
    pair    <- make_harmonise_pair()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    ## rs4: src allele=A matched trgt oth.allele=A → negate beta, swap labels
    expect_true("rs4" %in% res$src$snp)
    expect_equal(res$src[snp == "rs4", beta],      -0.4)
    expect_equal(res$src[snp == "rs4", allele],     "G")
    expect_equal(res$src[snp == "rs4", oth.allele], "A")
})

test_that("harmonise_gwas drops SNPs with irreconcilable allele mismatches from both GWAS", {
    pair    <- make_harmonise_pair()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    expect_false("rs5" %in% res$src$snp)
    expect_false("rs5" %in% res$trgt$snp)
})

test_that("harmonise_gwas alleles are identical in src and trgt after harmonisation", {
    pair    <- make_harmonise_pair()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    expect_identical(res$src$allele,     res$trgt$allele)
    expect_identical(res$src$oth.allele, res$trgt$oth.allele)
})

test_that("harmonise_gwas drops SNPs with missing oth.allele when no resolver is supplied", {
    pair  <- make_harmonise_pair()
    ## Blank out oth.allele for rs2 in src
    pair$src[snp == "rs2", oth.allele := ""]
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    expect_false("rs2" %in% res$src$snp)
    expect_false("rs2" %in% res$trgt$snp)
})

test_that("harmonise_gwas output is sorted by dbsnpid in both datasets", {
    pair <- make_harmonise_pair()
    ## Shuffle src rows to confirm sorting is enforced
    pair$src  <- pair$src[c(3L, 1L, 5L, 2L, 4L)]
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    res <- harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    expect_equal(res$src$dbsnpid,  sort(res$src$dbsnpid))
    expect_equal(res$trgt$dbsnpid, sort(res$trgt$dbsnpid))
})

test_that("harmonise_gwas writes a locus plot PNG to output_dir", {
    pair    <- make_harmonise_pair()
    locus   <- make_test_locus()
    out_dir <- tempfile(); dir.create(out_dir)

    harmonise_gwas(pair$src, pair$trgt, locus, out_dir)

    expect_true(file.exists(file.path(out_dir, "locus.plot.png")))
})
