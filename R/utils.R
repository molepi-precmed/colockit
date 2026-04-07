#' Check that required columns are present in a data table
#'
#' Stops with an informative error message listing any missing columns.
#'
#' @param cols Character vector of column names present in the data.
#' @param cols_expect Character vector of column names that must be present.
#'
#' @return Invisibly `NULL`; called for its side-effect.
#' @keywords internal
check_missing_cols <- function(cols, cols_expect) {
    cols_missing <- cols_expect[!cols_expect %in% cols]
    if (length(cols_missing) != 0) {
        stop("Required columns missing: ",
             paste(cols_missing, collapse = ", "), ".")
    }
    invisible(NULL)
}

#' Sanity-check betas, p-values, and MAF in a GWAS region
#'
#' Removes SNPs with `beta == 0` or p-values outside `(0, 1)`.
#' P-values that are exactly 0 are clamped to `.Machine$double.xmin`.
#' If a `maf` column is present, SNPs with MAF outside `(0, 1)` are removed.
#'
#' @param region A `data.table` containing the GWAS region. Required columns:
#'   `dbsnpid`, `beta`, `pvalue`. Optional column: `maf`.
#'
#' @return A filtered `data.table`.
#' @keywords internal
check_region <- function(region) {
    if (length(region[beta == 0, dbsnpid]) != 0) {
        cat("The following SNPs have beta = 0 and will be removed:\n")
        print(region[beta == 0, dbsnpid])
        region <- region[beta != 0]
    }
    if (length(region[pvalue >= 1, dbsnpid]) != 0) {
        cat("The following SNPs have pvalue >= 1 and will be removed:\n")
        print(region[pvalue <= 0 | pvalue >= 1, dbsnpid])
        region <- region[pvalue > 0 & pvalue < 1]
    }
    if (length(region[pvalue == 0, dbsnpid]) != 0) {
        cat("The following SNPs have pvalue == 0; clamping to .Machine$double.xmin:\n")
        print(region[pvalue == 0, dbsnpid])
        region[pvalue == 0, pvalue := .Machine$double.xmin]
    }
    if ("maf" %in% names(region)) {
        wrong_maf <- region[maf <= 0 | maf >= 1]
        if (nrow(wrong_maf) > 0) {
            cat("SNPs with MAF outside (0, 1) detected and removed:\n")
            print(wrong_maf)
            region <- region[maf > 0 & maf < 1]
        }
    }
    region
}

#' Add standard error and variance of beta to a GWAS region
#'
#' If `se` is absent it is derived from `beta` and `pvalue` via the
#' t-distribution with `n_samples - 2` degrees of freedom.
#' `varbeta` is always (re-)computed as `se^2`.
#'
#' @param region A `data.table` with columns `beta`, `pvalue`, `n_samples`.
#'   An `se` column may be present; if so it is used directly.
#'
#' @return The same `data.table` with `se` and `varbeta` columns added or
#'   updated.
#' @keywords internal
get_se_varbeta <- function(region) {
    if (!"se" %in% names(region))
        region[, se := abs(beta / qt(pvalue, df = n_samples - 2))]
    region[, varbeta := se^2]
    region
}
