#' Prepare GWAS regions — standalone backend
#'
#' Accepts user-supplied summary statistics directly with no database
#' dependency. LD is computed from a PLINK reference panel via [get_ldmat()].
#'
#' @param locus A single-row `data.table` with columns `chrom` (integer),
#'   `startpos` (integer), `endpos` (integer).
#' @param gwas1 A `data.table` of summary statistics for trait 1. Required
#'   columns: `snp` (rsid), `dbsnpid` (integer), `chrom`, `pos`, `allele`
#'   (effect allele), `oth.allele`, `beta`, `pvalue`, `n_samples`, `freq`
#'   (effect allele frequency), `trait.name`.
#' @param gwas2 A `data.table` of summary statistics for trait 2 (same
#'   column requirements as `gwas1`).
#' @param output_dir Full path to the output directory.
#' @param refplinkfile Full path to the PLINK reference panel (without
#'   extension).
#'
#' @return A named list with entries `src`, `trgt` (harmonised `data.table`s)
#'   and `ldmat` (LD correlation matrix), or `NULL` if no common SNPs remain.
#' @keywords internal
prepare_regions_standalone <- function(locus, gwas1, gwas2,
                                        output_dir, refplinkfile) {
    ## ------ input validation ------------------------------------------
    .check_locus(locus)
    if (!inherits(gwas1, "data.table"))
        stop("`gwas1` must be a data.table.")
    if (!inherits(gwas2, "data.table"))
        stop("`gwas2` must be a data.table.")

    required <- c("snp", "dbsnpid", "chrom", "pos", "allele", "oth.allele",
                  "beta", "pvalue", "n_samples", "freq", "trait.name")
    check_missing_cols(names(gwas1), required)
    check_missing_cols(names(gwas2), required)

    .check_refplink(refplinkfile)

    ## ------ subset to locus -------------------------------------------
    src  <- gwas1[chrom == locus$chrom & pos > locus$startpos &
                      pos < locus$endpos][order(pos)]
    trgt <- gwas2[chrom == locus$chrom & pos > locus$startpos &
                      pos < locus$endpos][order(pos)]

    if (nrow(src)  == 0) stop("gwas1 contains no SNPs in the specified locus.")
    if (nrow(trgt) == 0) stop("gwas2 contains no SNPs in the specified locus.")

    cat(sprintf("Standalone backend\n Trait 1: %s\n Trait 2: %s\n",
                unique(src$trait.name), unique(trgt$trait.name)))

    ## ------ harmonise -------------------------------------------------
    ## Standalone backend has no DB resolver for missing alleles; SNPs with
    ## missing oth.allele are simply dropped (find_alleles_fn = NULL).
    region <- harmonise_gwas(src, trgt, locus, output_dir,
                              find_alleles_fn = NULL)
    src  <- region$src
    trgt <- region$trgt

    ## ------ se / varbeta ----------------------------------------------
    src  <- get_se_varbeta(src)
    trgt <- get_se_varbeta(trgt)

    ## ------ common SNPs -----------------------------------------------
    common_snps <- intersect(src$snp, trgt$snp)
    if (length(common_snps) == 0) {
        cat("No common SNPs remain between gwas1 and gwas2 in the locus.\n")
        return(NULL)
    }
    cat(length(common_snps), "common SNPs before LD filtering.\n")

    ## ------ LD matrix -------------------------------------------------
    common_region <- src[snp %in% common_snps, .(snp, dbsnpid, allele)]
    ld     <- get_ldmat(common_region, refplinkfile)
    if (is.null(ld)) return(NULL)

    src  <- src[ld$region[, .(snp, dbsnpid)],  on = c("snp", "dbsnpid"),
                nomatch = NULL]
    trgt <- trgt[ld$region[, .(snp, dbsnpid)], on = c("snp", "dbsnpid"),
                 nomatch = NULL]
    ldmat <- ld$ldmat

    ## ------ sanity checks ---------------------------------------------
    src  <- check_region(src)
    trgt <- check_region(trgt)

    common_snps <- intersect(src$snp, trgt$snp)
    src   <- src[snp  %in% common_snps]
    trgt  <- trgt[snp %in% common_snps]
    ldmat <- ldmat[common_snps, common_snps]
    cat(length(common_snps), "common SNPs after all filtering.\n")

    ## ------ save ------------------------------------------------------
    data.table::fwrite(src,  file = file.path(output_dir, "gwas1.csv"))
    data.table::fwrite(trgt, file = file.path(output_dir, "gwas2.csv"))
    save(ldmat, file = file.path(output_dir, "ld.matrix.Rdata.gz"))
    cat("Prepared regions saved to", output_dir, "\n")

    list(src = src, trgt = trgt, ldmat = ldmat)
}
