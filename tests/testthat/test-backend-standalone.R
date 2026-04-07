## Tests for the standalone backend
##
## Tier 0 — no external dependencies: input validation with synthetic data.
## Tier 1 — requires the 1000 Genomes EUR reference panel at
##   /Users/andrico/Sites/tests-data/kg.2020.hg38.eur
##   (skipped gracefully when the panel is absent).
##
## The 12 SNPs below are confirmed polymorphic (MAF >= 0.05) in the panel,
## located in a 10 kb window on chr22 (29.00–29.01 Mb):
##
##   snp          pos       allele  oth.allele  MAF
##   rs997567     29002010  T       C           0.342
##   rs7289239    29003256  C       T           0.448
##   rs8142788    29004527  A       G           0.172
##   rs469994     29005345  A       G           0.475
##   rs132565     29005510  A       G           0.475
##   rs2347789    29006326  A       G           0.477
##   rs4333023    29006501  T       C           0.477
##   rs5997427    29006551  C       T           0.074
##   rs9613768    29006844  T       G           0.404
##   rs8137499    29007062  C       A           0.476
##   rs9608719    29008663  G       A           0.499
##   rs4580479    29009968  T       G           0.469

panel_path  <- Sys.getenv("COLOCKIT_REFPLINK",
                          unset = "/Users/andrico/Sites/tests-data/kg.2020.hg38.eur")
panel_available <- file.exists(paste0(panel_path, ".bim"))

## ---- shared helpers -------------------------------------------------------

## Build a synthetic GWAS for the 12 chr22 panel SNPs
make_panel_gwas <- function(trait = "trait_A", seed = 1L) {
    set.seed(seed)
    data.table::data.table(
        snp        = c("rs997567","rs7289239","rs8142788","rs469994",
                       "rs132565","rs2347789","rs4333023","rs5997427",
                       "rs9613768","rs8137499","rs9608719","rs4580479"),
        dbsnpid    = as.integer(c(997567L, 7289239L, 8142788L, 469994L,
                                   132565L, 2347789L, 4333023L, 5997427L,
                                   9613768L, 8137499L, 9608719L, 4580479L)),
        chrom      = 22L,
        pos        = as.integer(c(29002010, 29003256, 29004527, 29005345,
                                   29005510, 29006326, 29006501, 29006551,
                                   29006844, 29007062, 29008663, 29009968)),
        allele     = c("T","C","A","A","A","A","T","C","T","C","G","T"),
        oth.allele = c("C","T","G","G","G","G","C","T","G","A","A","G"),
        beta       = rnorm(12L, sd = 0.15),
        pvalue     = runif(12L, 1e-5, 0.5),
        n_samples  = 5000L,
        freq       = c(0.342, 0.448, 0.172, 0.475, 0.475,
                       0.477, 0.477, 0.074, 0.404, 0.476, 0.499, 0.469),
        trait.name = trait
    )
}

## Locus spanning all 12 panel SNPs
make_panel_locus <- function() {
    data.table::data.table(
        chrom    = 22L,
        startpos = 29001000L,
        endpos   = 29011000L
    )
}

## ---- Tier 0: input validation (no panel needed) ---------------------------

test_that("prepare_regions_standalone errors when gwas1 is not a data.table", {
    expect_error(
        colockit:::prepare_regions_standalone(
            locus        = make_panel_locus(),
            gwas1        = list(snp = "rs1"),   # not a data.table
            gwas2        = make_panel_gwas("B"),
            output_dir   = tempdir(),
            refplinkfile = panel_path),
        regexp = "data.table")
})

test_that("prepare_regions_standalone errors when gwas2 is not a data.table", {
    expect_error(
        colockit:::prepare_regions_standalone(
            locus        = make_panel_locus(),
            gwas1        = make_panel_gwas("A"),
            gwas2        = list(snp = "rs1"),
            output_dir   = tempdir(),
            refplinkfile = panel_path),
        regexp = "data.table")
})

test_that("prepare_regions_standalone errors on missing required columns in gwas1", {
    bad <- make_panel_gwas("A")
    bad[, freq := NULL]
    expect_error(
        colockit:::prepare_regions_standalone(
            locus        = make_panel_locus(),
            gwas1        = bad,
            gwas2        = make_panel_gwas("B"),
            output_dir   = tempdir(),
            refplinkfile = panel_path),
        regexp = "freq")
})

test_that("prepare_regions_standalone errors on missing required columns in gwas2", {
    bad <- make_panel_gwas("B")
    bad[, n_samples := NULL]
    expect_error(
        colockit:::prepare_regions_standalone(
            locus        = make_panel_locus(),
            gwas1        = make_panel_gwas("A"),
            gwas2        = bad,
            output_dir   = tempdir(),
            refplinkfile = panel_path),
        regexp = "n_samples")
})

test_that("prepare_regions_standalone errors when gwas1 has no SNPs in locus", {
    gwas1_off <- make_panel_gwas("A")
    gwas1_off[, chrom := 1L]   # wrong chromosome
    expect_error(
        colockit:::prepare_regions_standalone(
            locus        = make_panel_locus(),
            gwas1        = gwas1_off,
            gwas2        = make_panel_gwas("B"),
            output_dir   = tempdir(),
            refplinkfile = panel_path),
        regexp = "no SNPs")
})

test_that("prepare_regions_standalone errors when gwas2 has no SNPs in locus", {
    gwas2_off <- make_panel_gwas("B")
    gwas2_off[, pos := pos + 1e8L]   # shifted far outside locus
    expect_error(
        colockit:::prepare_regions_standalone(
            locus        = make_panel_locus(),
            gwas1        = make_panel_gwas("A"),
            gwas2        = gwas2_off,
            output_dir   = tempdir(),
            refplinkfile = panel_path),
        regexp = "no SNPs")
})

## ---- Tier 1: happy path with real reference panel -------------------------

test_that("prepare_regions_standalone returns named list with real panel", {
    skip_if_not(panel_available, "Reference panel not available")

    out_dir <- tempfile()
    res <- prepare_regions(
        locus        = make_panel_locus(),
        gwas1        = make_panel_gwas("trait_A", seed = 1L),
        gwas2        = make_panel_gwas("trait_B", seed = 2L),
        output_dir   = out_dir,
        refplinkfile = panel_path,
        backend      = "standalone")

    expect_false(is.null(res))
    expect_named(res, c("src", "trgt", "ldmat"))
    expect_true(data.table::is.data.table(res$src))
    expect_true(data.table::is.data.table(res$trgt))
    expect_true(is.matrix(res$ldmat))
})

test_that("prepare_regions_standalone ldmat is square and SNP-consistent", {
    skip_if_not(panel_available, "Reference panel not available")

    out_dir <- tempfile()
    res <- prepare_regions(
        locus        = make_panel_locus(),
        gwas1        = make_panel_gwas("trait_A", seed = 1L),
        gwas2        = make_panel_gwas("trait_B", seed = 2L),
        output_dir   = out_dir,
        refplinkfile = panel_path,
        backend      = "standalone")

    skip_if(is.null(res))
    n <- nrow(res$src)
    expect_equal(nrow(res$ldmat), n)
    expect_equal(ncol(res$ldmat), n)
    expect_equal(rownames(res$ldmat), res$src$snp)
    expect_equal(colnames(res$ldmat), res$trgt$snp)
})

test_that("prepare_regions_standalone writes output files", {
    skip_if_not(panel_available, "Reference panel not available")

    out_dir <- tempfile()
    res <- prepare_regions(
        locus        = make_panel_locus(),
        gwas1        = make_panel_gwas("trait_A", seed = 1L),
        gwas2        = make_panel_gwas("trait_B", seed = 2L),
        output_dir   = out_dir,
        refplinkfile = panel_path,
        backend      = "standalone")

    skip_if(is.null(res))
    expect_true(file.exists(file.path(out_dir, "gwas1.csv")))
    expect_true(file.exists(file.path(out_dir, "gwas2.csv")))
    expect_true(file.exists(file.path(out_dir, "ld.matrix.Rdata.gz")))
})

## ---- End-to-end smoke test: prepare_regions → colocalise(ABF) ------------

test_that("full standalone pipeline runs end-to-end and posteriors sum to 1", {
    skip_if_not(panel_available, "Reference panel not available")

    out_dir <- tempfile()

    ## Stage 1 — prepare regions
    regions <- prepare_regions(
        locus        = make_panel_locus(),
        gwas1        = make_panel_gwas("trait_A", seed = 10L),
        gwas2        = make_panel_gwas("trait_B", seed = 20L),
        output_dir   = out_dir,
        refplinkfile = panel_path,
        backend      = "standalone")

    skip_if(is.null(regions), "No common SNPs after LD filtering")

    ## Stage 2 — colocalise with ABF
    res <- colocalise(
        regions      = regions,
        locus        = make_panel_locus(),
        analysis_type = "abf",
        output_dir   = out_dir,
        refplinkfile = panel_path,
        cutoff       = 0.0)   # cutoff=0 so we always get a result back

    expect_false(is.null(res))

    ## posteriors must sum to 1 (within floating-point tolerance)
    pp_cols <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
    pp_sum  <- sum(res$summary[pp_cols])
    expect_equal(pp_sum, 1.0, tolerance = 1e-6)

    ## plot file produced
    expect_true(file.exists(file.path(out_dir, "coloc.abf.png")))
})
