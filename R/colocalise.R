#' Run colocalisation analysis on prepared GWAS regions
#'
#' Takes the output of [prepare_regions()] and runs the chosen
#' colocalisation method. Saves formatted input data and a colocalisation
#' plot to `output_dir`.
#'
#' @param regions A named list with entries `src`, `trgt` (harmonised GWAS
#'   `data.table`s) and `ldmat` (LD correlation matrix), as returned by
#'   [prepare_regions()].
#' @param locus A single-row `data.table` with columns `chrom`, `startpos`,
#'   `endpos`.
#' @param analysis_type Colocalisation method. One of:
#'   * `"abf"` — coloc Approximate Bayes Factor
#'   * `"susie"` — coloc SuSiE (handles multiple causal variants)
#'   * `"pwcoco"` — PWCoCo (external command-line tool)
#' @param output_dir Full path to the directory where results will be saved.
#' @param refplinkfile Full path to the PLINK reference panel. Required for
#'   `"pwcoco"`. Default `NULL`.
#' @param cutoff P-value cut-off for the stepwise SNP selection in PWCoCo.
#'   Default `0.05`.
#'
#' @return For `"abf"` / `"susie"`: the coloc result object (list with
#'   `$summary` and `$results`). For `"pwcoco"`: a `data.table` of PWCoCo
#'   output.  Returns `NULL` if SuSiE fails.
#'
#' @examples
#' \dontrun{
#' ## After running prepare_regions():
#' result <- colocalise(
#'   regions       = regions,
#'   locus         = locus,
#'   analysis_type = "abf",
#'   output_dir    = "path/to/output"
#' )
#'
#' ## Inspect posterior probabilities:
#' result$summary
#'
#' ## SuSiE (handles multiple causal variants):
#' result_susie <- colocalise(
#'   regions       = regions,
#'   locus         = locus,
#'   analysis_type = "susie",
#'   output_dir    = "path/to/output"
#' )
#' }
#'
#' @export
#' @importFrom coloc coloc.abf coloc.susie
#' @importFrom data.table fwrite fread setDT setnames
colocalise <- function(regions,
                        locus,
                        analysis_type = c("abf", "susie", "pwcoco"),
                        output_dir,
                        refplinkfile = NULL,
                        cutoff = 0.05) {
    analysis_type <- match.arg(analysis_type)
    output_dir    <- .valid_dir(output_dir)

    src   <- regions$src
    trgt  <- regions$trgt
    ldmat <- regions$ldmat

    src_fmt  <- format_region(src,  analysis_type, ldmat)
    trgt_fmt <- format_region(trgt, analysis_type, ldmat)

    formatted <- list(src = src_fmt, trgt = trgt_fmt)
    save(formatted,
         file = file.path(output_dir, "formatted.region.Rdata.gz"))
    cat("Formatted regions saved to", output_dir, "\n")

    ## ------------------------------------------------------------------ ##
    ##  ABF                                                               ##
    ## ------------------------------------------------------------------ ##
    if (analysis_type == "abf") {
        result <- coloc::coloc.abf(src_fmt, trgt_fmt)
        plot_file <- file.path(output_dir, "coloc.abf.png")
        plot_regions(src, trgt, locus, plot_file,
                     coloc = result$results, ldmat = ldmat)
        print(result$summary)
        return(invisible(result))
    }

    ## ------------------------------------------------------------------ ##
    ##  SuSiE                                                             ##
    ## ------------------------------------------------------------------ ##
    if (analysis_type == "susie") {
        result <- coloc::coloc.susie(src_fmt, trgt_fmt)
        if (!"summary" %in% names(result)) {
            cat("SuSiE failed to complete colocalisation analysis.\n")
            return(invisible(NULL))
        }
        ## select the credible set with the highest PP.H4
        best_idx  <- which.max(result$summary$PP.H4.abf)
        cols      <- c(1L, 1L + best_idx)
        best_res  <- result$results[, ..cols]
        data.table::setnames(best_res, 2L, "SNP.PP.H4")

        plot_file <- file.path(output_dir, "coloc.susie.png")
        plot_regions(src, trgt, locus, plot_file,
                     coloc = best_res, ldmat = ldmat)
        print(result$summary)
        return(invisible(result))
    }

    ## ------------------------------------------------------------------ ##
    ##  PWCoCo                                                            ##
    ## ------------------------------------------------------------------ ##
    if (analysis_type == "pwcoco") {
        if (is.null(refplinkfile))
            stop("`refplinkfile` must be provided for PWCoCo analysis.")

        wd <- getwd()
        on.exit(setwd(wd), add = TRUE)
        setwd(output_dir)

        if (file.exists("pwcoco_out.coloc"))
            stop("PWCoCo output already exists in '", output_dir, "'. ",
                 "Remove it or choose a different output directory.")

        gwas1_file <- file.path(output_dir, "pwcoco.gwas1.csv")
        gwas2_file <- file.path(output_dir, "pwcoco.gwas2.csv")
        data.table::fwrite(src_fmt,  file = gwas1_file)
        data.table::fwrite(trgt_fmt, file = gwas2_file)

        cmd <- sprintf(
            "pwcoco --bfile %s --sum_stats1 %s --sum_stats2 %s \\\n                --maf 0.05 --p_cutoff %s --chr %s --verbose",
            refplinkfile, gwas1_file, gwas2_file, cutoff, locus$chrom)
        system(cmd)

        result <- data.table::fread("pwcoco_out.coloc", header = TRUE)
        return(invisible(result))
    }
}

## ---- placeholder stubs for future methods ----------------------------

#' @keywords internal
.method_hyprcoloc <- function(...) {
    stop("HyPrColoc support is not yet implemented in this version of colockit.")
}

#' @keywords internal
.method_ecaviar <- function(...) {
    stop("eCAVIAR support is not yet implemented in this version of colockit.")
}

#' @keywords internal
.method_paintor <- function(...) {
    stop("PAINTOR support is not yet implemented in this version of colockit.")
}

#' @keywords internal
.method_finemap <- function(...) {
    stop("FINEMAP support is not yet implemented in this version of colockit.")
}
