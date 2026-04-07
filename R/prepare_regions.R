#' Prepare GWAS regions for colocalisation analysis
#'
#' Stage 1 of the colockit pipeline. Harmonises two GWAS datasets within a
#' specified genomic locus, validates allele frequencies, computes an LD
#' matrix, and writes intermediate files to `output_dir`.
#'
#' Dispatches to one of two backends depending on the `backend` argument:
#' * `"standalone"` — no database dependency; all inputs supplied directly.
#' * `"genoscores"` — integrates with the Genoscores database to retrieve
#'   published GWAS or generate private GWAS via SNPTEST.
#'
#' @param locus A single-row `data.table` with integer columns `chrom`,
#'   `startpos`, and `endpos` defining the genomic region of interest.
#' @param gwas1 For `backend = "standalone"`: a `data.table` of summary
#'   statistics for trait 1 (see Details).
#'   For `backend = "genoscores"`: an integer gwasid, or a `data.table` of
#'   private summary statistics when `gwas_type = "private"`.
#' @param gwas2 Same as `gwas1` but for trait 2.
#' @param output_dir Full path to the directory where intermediate files will
#'   be saved (`gwas1.csv`, `gwas2.csv`, `ld.matrix.Rdata.gz`,
#'   `locus.plot.png`).
#' @param refplinkfile Full path to the PLINK reference panel used to compute
#'   the LD matrix (without file extension; `.bed`, `.bim`, `.fam` must
#'   exist).
#' @param backend Which data backend to use. One of `"standalone"` (default)
#'   or `"genoscores"`.
#' @param gwas_type *(genoscores backend only)* One of `"published"`,
#'   `"mixed"`, or `"private"`.
#' @param con *(genoscores backend only)* A DBI database connection to the
#'   Genoscores database.
#'
#' @details
#' **Standalone backend — required columns for `gwas1` / `gwas2`:**
#' \describe{
#'   \item{`snp`}{rsid string, e.g. `"rs123456"`}
#'   \item{`dbsnpid`}{integer dbSNP identifier}
#'   \item{`chrom`}{integer chromosome number}
#'   \item{`pos`}{integer genomic position (hg38)}
#'   \item{`allele`}{effect allele (A/C/G/T)}
#'   \item{`oth.allele`}{other allele}
#'   \item{`beta`}{effect size}
#'   \item{`pvalue`}{association p-value}
#'   \item{`n_samples`}{sample size}
#'   \item{`freq`}{effect allele frequency}
#'   \item{`trait.name`}{character label for the trait}
#' }
#'
#' @return A named list with entries:
#'   * `src` — harmonised `data.table` for trait 1
#'   * `trgt` — harmonised `data.table` for trait 2
#'   * `ldmat` — LD correlation matrix (rsid-named)
#'
#'   Returns `NULL` if no common SNPs remain after all filtering steps.
#'
#' @seealso [colocalise()] for Stage 2.
#'
#' @examples
#' \dontrun{
#' ## Standalone backend ------------------------------------------------
#' locus <- data.table::data.table(
#'   chrom    = 22L,
#'   startpos = 29001000L,
#'   endpos   = 29011000L
#' )
#'
#' ## gwas1 and gwas2 must be data.tables with columns:
#' ## snp, dbsnpid, chrom, pos, allele, oth.allele,
#' ## beta, pvalue, n_samples, freq, trait.name
#' regions <- prepare_regions(
#'   locus        = locus,
#'   gwas1        = gwas1,
#'   gwas2        = gwas2,
#'   output_dir   = "path/to/output",
#'   refplinkfile = "path/to/kg.2020.hg38.eur",
#'   backend      = "standalone"
#' )
#'
#' ## Genoscores backend ------------------------------------------------
#' con <- DBI::dbConnect(...)
#' regions <- prepare_regions(
#'   locus        = locus,
#'   gwas1        = 42L,   # integer gwasid
#'   gwas2        = 57L,
#'   output_dir   = "path/to/output",
#'   refplinkfile = "path/to/kg.2020.hg38.eur",
#'   backend      = "genoscores",
#'   gwas_type    = "published",
#'   con          = con
#' )
#' }
#'
#' @export
prepare_regions <- function(locus,
                             gwas1,
                             gwas2,
                             output_dir,
                             refplinkfile,
                             backend   = c("standalone", "genoscores"),
                             gwas_type = c("published", "mixed", "private"),
                             con       = NULL) {
    backend   <- match.arg(backend)
    gwas_type <- match.arg(gwas_type)

    if (backend == "standalone") {
        prepare_regions_standalone(
            locus        = locus,
            gwas1        = gwas1,
            gwas2        = gwas2,
            output_dir   = output_dir,
            refplinkfile = refplinkfile
        )
    } else {
        if (is.null(con))
            stop("A DBI database connection (`con`) must be supplied for ",
                 "the genoscores backend.")
        prepare_regions_genoscores(
            locus        = locus,
            gwas1        = gwas1,
            gwas2        = gwas2,
            gwas_type    = gwas_type,
            output_dir   = output_dir,
            refplinkfile = refplinkfile,
            con          = con
        )
    }
}

## ======================================================================
##  Shared internal validators
## ======================================================================

#' @keywords internal
.check_locus <- function(locus) {
    if (!data.table::is.data.table(locus))
        stop("`locus` must be a data.table.")
    check_missing_cols(names(locus), c("chrom", "startpos", "endpos"))
}

#' @keywords internal
.check_refplink <- function(refplinkfile) {
    if (is.null(refplinkfile))
        stop("`refplinkfile` must be provided.")
    if (!file.exists(paste0(refplinkfile, ".bed")))
        stop("PLINK .bed file not found: ", paste0(refplinkfile, ".bed"))
}

#' @keywords internal
.valid_dir <- function(path) {
    if (!dir.exists(path))
        dir.create(path, recursive = TRUE)
    normalizePath(path)
}
