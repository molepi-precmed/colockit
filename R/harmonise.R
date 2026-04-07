#' Harmonise effect alleles between two GWAS regions
#'
#' Ensures that effect alleles are aligned between `src` and `trgt`.
#' Where the effect allele differs, `beta` and the allele labels in `src`
#' are flipped to match `trgt`. SNPs with irreconcilable allele mismatches
#' are dropped from both datasets.
#'
#' Calls [plot_regions()] as a side-effect to produce a preliminary locus
#' plot highlighting common SNPs.
#'
#' @param src A `data.table` of GWAS summary statistics for trait 1.
#' @param trgt A `data.table` of GWAS summary statistics for trait 2.
#' @param locus A single-row `data.table` with columns `chrom`, `startpos`,
#'   `endpos` defining the region.
#' @param output_dir Full path to the output directory.
#' @param find_alleles_fn A function used to resolve missing `oth.allele`
#'   values. Defaults to `NULL` (missing alleles cause those SNPs to be
#'   dropped). The genoscores backend supplies its own DB-backed resolver.
#'   The function must accept a `data.table` of rows with missing alleles
#'   and return a `data.table` with columns `dbsnpid`, `allele`,
#'   `oth.allele`.
#'
#' @return A named list with entries `src` and `trgt`, each a harmonised
#'   `data.table` containing only SNPs common to both GWAS.
#' @keywords internal
harmonise_gwas <- function(src, trgt, locus, output_dir,
                            find_alleles_fn = NULL) {
    plot_file <- file.path(output_dir, "locus.plot.png")
    region    <- plot_regions(src, trgt, locus, plot_file)
    src   <- region$src
    trgt  <- region$trgt

    ## ------------------------------------------------------------------ ##
    ##  Resolve missing oth.allele                                        ##
    ## ------------------------------------------------------------------ ##
    src  <- .resolve_other_allele(src,  find_alleles_fn)
    trgt <- .resolve_other_allele(trgt, find_alleles_fn)

    ## keep only SNPs present in both after allele resolution
    common_ids <- intersect(src$dbsnpid, trgt$dbsnpid)
    src  <- src[dbsnpid  %in% common_ids]
    trgt <- trgt[dbsnpid %in% common_ids]

    ## ------------------------------------------------------------------ ##
    ##  Align effect alleles: flip beta + labels in src where needed      ##
    ## ------------------------------------------------------------------ ##
    src[trgt, on = c("snp", allele = "allele"),  newbeta := beta]
    src[trgt, on = c("snp", allele = "oth.allele"),
        c("allele", "oth.allele", "newbeta") := list(i.allele, i.oth.allele, -beta)]

    wrong_alleles <- src[is.na(newbeta)]
    if (nrow(wrong_alleles) > 0) {
        cat(nrow(wrong_alleles),
            "variant(s) with incorrectly coded alleles removed.\n")
        bad <- wrong_alleles$dbsnpid
        src  <- src[!dbsnpid  %in% bad]
        trgt <- trgt[!dbsnpid %in% bad]
    }

    src[, beta    := newbeta]
    src[, newbeta := NULL]

    ## enforce integer dbsnpid and sort both by position
    src[,  dbsnpid := as.integer(dbsnpid)]
    trgt[, dbsnpid := as.integer(dbsnpid)]
    src  <- src[order(dbsnpid)]
    trgt <- trgt[order(dbsnpid)]

    stopifnot(identical(src$allele,     trgt$allele))
    stopifnot(identical(src$oth.allele, trgt$oth.allele))

    list(src = src, trgt = trgt)
}

## Internal helper: fill missing oth.allele using the supplied resolver
.resolve_other_allele <- function(gwas, find_alleles_fn) {
    gwas[, oth.allele := as.character(oth.allele)]
    missing <- gwas[oth.allele == "" | is.na(oth.allele)]
    if (nrow(missing) == 0) return(gwas)

    if (is.null(find_alleles_fn)) {
        cat(nrow(missing),
            "SNP(s) with missing oth.allele dropped",
            "(no allele resolver supplied).\n")
        return(gwas[!(oth.allele == "" | is.na(oth.allele))])
    }

    found <- find_alleles_fn(missing)
    gwas[found, on = c(dbsnpid = "dbsnpid", allele = "allele"),
         oth.allele := i.oth.allele]
    gwas <- gwas[!(oth.allele == "" | is.na(oth.allele))]
    gwas
}
