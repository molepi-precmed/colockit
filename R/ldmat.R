#' Compute an LD correlation matrix from a PLINK reference panel
#'
#' Reads a PLINK binary reference panel using `bigsnpr`, filters to SNPs
#' present in `region`, removes monomorphic variants, and returns a
#' signed LD correlation matrix with allele orientation matched to the
#' effect alleles in `region`.
#'
#' This function is shared by both the standalone and genoscores backends.
#'
#' @param region A `data.table` with columns `snp` (rsid, e.g. `"rs123"`),
#'   `dbsnpid` (integer), and `allele` (effect allele, one of A/C/G/T).
#' @param refplinkfile Full path to the PLINK reference panel (without
#'   extension; `.bed`, `.bim`, `.fam` files must exist at that path).
#'
#' @return A named list with two entries:
#'   * `region` — the input `data.table` filtered to SNPs that could be
#'     matched in the reference panel and are not monomorphic.
#'   * `ldmat` — a symmetric numeric matrix of LD correlations (r, not r²),
#'     with row and column names set to rsids.
#'   Returns `NULL` if no SNPs from `region` are found in the reference.
#'
#' @importFrom bigsnpr snp_readBed snp_attach snp_cor snp_MAF
#' @importFrom data.table data.table setDT setnames fread fwrite
#' @keywords internal
get_ldmat <- function(region, refplinkfile) {
    check_missing_cols(names(region), c("snp", "dbsnpid", "allele"))

    if (!all(grepl("^rs", region$snp)))
        stop("All entries in the `snp` column must be rs-identifiers.")

    if (!is.integer(region$dbsnpid))
        stop("`dbsnpid` must be of class integer.")

    if (anyDuplicated(region$dbsnpid))
        stop("Duplicated SNPs detected. Please remove them before proceeding.")

    ## ------------------------------------------------------------------ ##
    ##  1. Read reference BIM (no genotypes) to find overlapping SNPs     ##
    ## ------------------------------------------------------------------ ##
    bim_file <- paste0(refplinkfile, ".bim")
    if (!file.exists(bim_file))
        stop("Reference .bim file not found: ", bim_file)

    ref_bim <- data.table::fread(
        bim_file,
        col.names  = c("chrom", "snp", "cm", "pos", "a1", "a2"),
        colClasses = list(character = c(2L, 5L, 6L))
    )

    ## keep only SNPs present in both region and reference
    region <- region[snp %in% ref_bim$snp]
    if (nrow(region) == 0) {
        cat("No SNPs from the region were found in the reference panel.\n")
        return(NULL)
    }
    ref_bim <- ref_bim[snp %in% region$snp]

    ## ------------------------------------------------------------------ ##
    ##  2. Read genotypes into a bigSNP object                            ##
    ## ------------------------------------------------------------------ ##
    plink_bin <- Sys.which("plink")
    if (nchar(plink_bin) == 0)
        stop("plink executable not found. ",
             "Ensure plink is installed and on your PATH.")

    tmp_bed <- tempfile(fileext = ".bed")
    rds_file <- sub("\\.bed$", ".rds", tmp_bed)

    snp_list_file <- tempfile()
    data.table::fwrite(region[, .(snp)], file = snp_list_file,
                       col.names = FALSE)
    cmd <- sprintf(
        "%s --silent --bfile %s --extract %s --make-bed --out %s",
        plink_bin, refplinkfile, snp_list_file,
        sub("\\.bed$", "", tmp_bed)
    )
    system(cmd)

    if (!file.exists(tmp_bed))
        stop("PLINK failed to extract the region from the reference panel.")

    snp_obj  <- bigsnpr::snp_readBed(tmp_bed, backingfile = sub("\\.bed$", "", rds_file))
    snp_data <- bigsnpr::snp_attach(snp_obj)

    geno  <- snp_data$genotypes   # FBM object (individuals x SNPs)
    bim   <- data.table::setDT(snp_data$map)
    data.table::setnames(bim,
                         c("chr.autosome", "marker.ID", "genetic.dist",
                           "physical.pos", "allele1", "allele2"),
                         c("chrom", "snp", "cm", "pos", "a1", "a2"),
                         skip_absent = TRUE)

    ## ------------------------------------------------------------------ ##
    ##  3. Remove monomorphic SNPs                                        ##
    ## ------------------------------------------------------------------ ##
    maf <- bigsnpr::snp_MAF(geno)
    mono_idx <- which(maf == 0)
    if (length(mono_idx) > 0) {
        mono_snps <- bim$snp[mono_idx]
        cat("Removed", length(mono_idx),
            "SNP(s) monomorphic in the reference panel.\n")
        region  <- region[!snp %in% mono_snps]
        keep_idx <- which(!bim$snp %in% mono_snps)
        bim   <- bim[keep_idx]
        geno_cols <- keep_idx          # column indices to pass to snp_cor
    } else {
        geno_cols <- seq_len(ncol(geno))
    }

    if (nrow(region) == 0) {
        cat("No polymorphic SNPs remain after monomorphic filtering.\n")
        return(NULL)
    }

    ## ------------------------------------------------------------------ ##
    ##  4. Compute LD correlation matrix                                  ##
    ## ------------------------------------------------------------------ ##
    cat("Computing LD matrix for", nrow(region), "SNPs ...\n")
    ldmat <- as.matrix(
        bigsnpr::snp_cor(geno, ind.col = geno_cols)
    )

    ## re-order rows/cols to match region order
    rownames(ldmat) <- bim$snp
    colnames(ldmat) <- bim$snp
    ldmat <- ldmat[region$snp, region$snp]

    ## ------------------------------------------------------------------ ##
    ##  5. Allele-sign correction                                         ##
    ##                                                                    ##
    ##  bigsnpr codes alleles as a1 (ALT/effect-coded) vs a2 (REF).       ##
    ##  If the effect allele in `region` matches a2 (the REF allele in    ##
    ##  the reference), the correlation column/row needs to be flipped    ##
    ##  (multiplied by -1) so that all correlations are expressed         ##
    ##  relative to the same effect allele strand.                        ##
    ## ------------------------------------------------------------------ ##
    flip <- region$allele != bim[match(region$snp, bim$snp), a1]
    if (any(flip)) {
        flip_idx <- which(flip)
        ldmat[flip_idx, ] <- -ldmat[flip_idx, ]
        ldmat[, flip_idx] <- -ldmat[, flip_idx]
    }

    ## set dimnames to rsids
    rownames(ldmat) <- region$snp
    colnames(ldmat) <- region$snp

    cat("Computed", paste0(nrow(ldmat), "x", ncol(ldmat)),
        "LD matrix.\n")

    ## clean up temp files
    unlink(c(tmp_bed,
             sub("\\.bed$", ".bim", tmp_bed),
             sub("\\.bed$", ".fam", tmp_bed),
             sub("\\.bed$", ".log", tmp_bed),
             rds_file,
             snp_list_file),
           force = TRUE)

    list(region = region, ldmat = ldmat)
}
