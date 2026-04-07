#' Format a harmonised GWAS region for the chosen colocalisation method
#'
#' Reshapes a harmonised GWAS `data.table` into the exact input structure
#' expected by each method.
#'
#' For `"abf"` and `"susie"` the output is a named list validated by
#' [coloc::check_dataset()]. For `"pwcoco"` the output is a flat
#' `data.table`.
#'
#' @param region A harmonised GWAS `data.table`. Required columns:
#'   `snp`, `pos`, `beta`, `varbeta`, `pvalue`, `freq`, `n_samples`.
#'   Additional columns `allele`, `oth.allele`, `se` are needed for
#'   `"pwcoco"`.
#' @param analysis_type One of `"abf"`, `"susie"`, or `"pwcoco"`.
#' @param ldmat An LD correlation matrix whose row/column names are rsids
#'   matching `region$snp`. May be `NULL` for `"abf"` (LD is optional for
#'   ABF).
#'
#' @return For `"abf"` / `"susie"`: a named list suitable for
#'   [coloc::coloc.abf()] or [coloc::coloc.susie()].
#'   For `"pwcoco"`: a `data.table` with columns `snp`, `eff.allele`,
#'   `oth.allele`, `freq`, `beta`, `se`, `pvalue`, `n_samples`.
#'
#' @importFrom coloc check_dataset
#' @keywords internal
format_region <- function(region, analysis_type, ldmat) {
    check_missing_cols(names(region),
                       c("snp", "beta", "freq", "pvalue", "n_samples"))

    if (analysis_type %in% c("abf", "susie")) {
        out <- list(
            snp      = region$snp,
            position = region$pos,
            beta     = region$beta,
            varbeta  = region$varbeta,
            pvalues  = region$pvalue,
            MAF      = region$freq,
            N        = region[1L, n_samples],
            type     = "quant"
        )
        ## Only include LD when provided; coloc::check_dataset errors on LD=NULL
        if (!is.null(ldmat)) out$LD <- ldmat
        coloc::check_dataset(out)
        return(out)
    }

    if (analysis_type == "pwcoco") {
        return(region[, .(snp,
                           eff.allele = allele,
                           oth.allele,
                           freq,
                           beta,
                           se,
                           pvalue,
                           n_samples)])
    }

    stop("Unknown analysis_type: '", analysis_type,
         "'. Must be one of 'abf', 'susie', 'pwcoco'.")
}
