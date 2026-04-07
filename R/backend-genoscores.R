#' Prepare GWAS regions — genoscores backend
#'
#' Integrates with the Genoscores database and infrastructure to retrieve
#' published GWAS summary statistics and/or generate private GWAS via SNPTEST.
#' All existing genoscores-coupled logic is preserved here unchanged.
#'
#' @param locus A single-row `data.table` with columns `chrom`, `startpos`,
#'   `endpos`.
#' @param gwas1 Integer gwasid if `gwas_type` is `"published"` or `"mixed"`;
#'   a `data.table` of private summary stats if `gwas_type` is `"private"`.
#' @param gwas2 Integer gwasid if `gwas_type` is `"published"`; a `data.table`
#'   of private summary stats if `gwas_type` is `"mixed"` or `"private"`.
#' @param gwas_type One of `"published"`, `"mixed"`, or `"private"`.
#' @param output_dir Full path to the output directory.
#' @param refplinkfile Full path to the PLINK reference panel.
#' @param con A DBI database connection to the Genoscores database.
#'
#' @return A named list with entries `src`, `trgt`, and `ldmat`, or `NULL`
#'   if no common SNPs remain.
#' @keywords internal
prepare_regions_genoscores <- function(locus, gwas1, gwas2,
                                        gwas_type = c("published",
                                                       "mixed",
                                                       "private"),
                                        output_dir, refplinkfile, con) {
    gwas_type  <- match.arg(gwas_type)
    output_dir <- .valid_dir(output_dir)
    .check_locus(locus)
    .check_refplink(refplinkfile)

    ## ------ load GWAS depending on type --------------------------------
    all_gwas <- .validate_input_genoscores(locus, gwas1, gwas2,
                                            gwas_type, refplinkfile, con)
    src  <- all_gwas$src
    trgt <- all_gwas$trgt
    cat(sprintf("%s\n Trait 1: %s\n Trait 2: %s\n",
                all_gwas$msg,
                gsub(" ", "_", unique(src$trait.name)),
                gsub(" ", "_", unique(trgt$trait.name))))

    ## ------ subset to locus -------------------------------------------
    src  <- src[chrom  == locus$chrom & pos > locus$startpos &
                    pos < locus$endpos][order(pos)]
    trgt <- trgt[chrom == locus$chrom & pos > locus$startpos &
                     pos < locus$endpos][order(pos)]

    if (nrow(src)  == 0 || nrow(trgt) == 0)
        stop("One or both GWAS contain no data in the specified locus.")

    ## ------ harmonise --------------------------------------------------
    ## Genoscores backend resolves missing alleles via dbSNP table
    resolver <- function(missing) find_missing_alleles_gs(missing, con)
    region <- harmonise_gwas(src, trgt, locus, output_dir,
                              find_alleles_fn = resolver)
    src  <- region$src
    trgt <- region$trgt

    ## ------ allele frequencies from DB --------------------------------
    freqs <- get_freqs_gs(src, con)
    src   <- freqs[src,  on = c("dbsnpid", "allele"), nomatch = NULL]
    trgt  <- freqs[trgt, on = c("dbsnpid", "allele"), nomatch = NULL]

    ## for private GWAS: cross-check SNPTEST frequencies against dbSNP
    if ("maf" %in% names(src))  src  <- compare_freqs(src)
    if ("maf" %in% names(trgt)) trgt <- compare_freqs(trgt)

    ## ------ se / varbeta ----------------------------------------------
    src  <- get_se_varbeta(src)
    trgt <- get_se_varbeta(trgt)

    ## ------ common SNPs -----------------------------------------------
    common_snps <- intersect(src$snp, trgt$snp)
    if (length(common_snps) == 0) {
        cat("No common SNPs remain between GWAS 1 and GWAS 2 in the locus.\n")
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

    ## ------ save -------------------------------------------------------
    data.table::fwrite(src,  file = file.path(output_dir, "gwas1.csv"))
    data.table::fwrite(trgt, file = file.path(output_dir, "gwas2.csv"))
    save(ldmat, file = file.path(output_dir, "ld.matrix.Rdata.gz"))
    cat("Prepared regions saved to", output_dir, "\n")

    list(src = src, trgt = trgt, ldmat = ldmat)
}

## ======================================================================
##  Genoscores-specific helpers (DB-coupled)
## ======================================================================

#' Validate user input and load GWAS — genoscores backend
#' @keywords internal
.validate_input_genoscores <- function(locus, src, trgt, gwas_type,
                                        refplinkfile, con) {
    if (gwas_type == "published") {
        if (!is.numeric(src) || !is.numeric(trgt))
            stop("`gwas1` and `gwas2` must be integer gwasids for ",
                 "gwas_type = 'published'.")
        msg      <- "Using published GWAS from Genoscores database."
        all_gwas <- get_published_gwas(c(src, trgt), con)
        src  <- all_gwas[gwasid == src]
        trgt <- all_gwas[gwasid == trgt]
    }
    if (gwas_type == "mixed") {
        if (!is.numeric(src))
            stop("`gwas1` must be an integer gwasid for gwas_type = 'mixed'.")
        if (!inherits(trgt, "data.table"))
            stop("`gwas2` must be a data.table for gwas_type = 'mixed'.")
        msg <- "Using published GWAS (gwas1) + private GWAS (gwas2)."
        src <- get_published_gwas(src, con)
    }
    if (gwas_type == "private") {
        if (!inherits(src, "data.table") || !inherits(trgt, "data.table"))
            stop("`gwas1` and `gwas2` must be data.tables for ",
                 "gwas_type = 'private'.")
        msg <- "Using private GWAS."
    }
    list(msg = msg, src = src, trgt = trgt)
}

#' Read published GWAS flatfiles from the Genoscores datastore
#'
#' @param gwasids Integer vector of gwasids to retrieve.
#' @param con DBI connection to the Genoscores database.
#' @keywords internal
get_published_gwas <- function(gwasids, con) {
    query <- sprintf(
        "select g.gwasid, t.mapped_value, t.trait_name,
         gm.intermediate_file, gm.trait_variable, gm.n_samples
         from gwas g
         inner join gwas_meta gm on g.gwasid = gm.gwasid
         inner join traits t on g.traitid = t.traitid
         where g.gwasid in (%s)",
        paste(sprintf("'%s'", gwasids), collapse = ", "))
    meta <- data.table::setDT(DBI::dbGetQuery(con, query))

    if (nrow(meta) == 0)
        stop("No traits found for the provided gwasids. ",
             "Check that they are correct.")

    intermediatedir <- "/opt/datastore/genoscores/flatfiles"
    gwas_list <- vector("list", nrow(meta))
    for (idx in seq_len(nrow(meta))) {
        this <- meta[idx]
        nm   <- paste0("X_", this$gwasid)
        dt   <- data.table::fread(
            file.path(intermediatedir, this$intermediate_file))
        dt[, snp := paste0("rs", dbsnpid)]
        if ("mapped.value"    %in% names(dt))
            dt <- dt[mapped.value    == this$mapped_value]
        if ("trait.variable"  %in% names(dt))
            dt <- dt[trait.variable  == this$trait_variable]
        dt <- dt[, .(snp, dbsnpid, chrom, pos, allele, oth.allele,
                      beta, pvalue)]
        ## check uniqueness
        dt[, g := .GRP, by = .(dbsnpid, chrom, pos)]
        if (nrow(dt[duplicated(g)]) > 0) {
            cat("Duplicated entries detected in gwasid", this$gwasid, "\n")
            print(dt[duplicated(g)])
            stop("Duplicates indicate possible GWAS integrity issues.")
        }
        dt[, g := NULL]
        dt[, `:=`(gwasid     = this$gwasid,
                   flatfile   = this$intermediate_file,
                   trait.name = this$trait_name,
                   n_samples  = this$n_samples)]
        gwas_list[[idx]] <- dt
    }
    data.table::rbindlist(gwas_list)
}

#' Retrieve effect allele frequencies from the Genoscores r_freqs table
#'
#' @param gwas A `data.table` with columns `dbsnpid` and `allele`.
#' @param con DBI connection to the Genoscores database.
#' @keywords internal
get_freqs_gs <- function(gwas, con) {
    check_missing_cols(names(gwas), c("dbsnpid", "allele"))
    cat("Retrieving allele frequencies from Genoscores database ...\n")
    genoscores:::dbTemporaryTable(con, "mygwas",
                                   gwas[, .(dbsnpid, allele)])
    DBI::dbExecute(con, "create index dbsnpid on mygwas (dbsnpid)")
    freqs <- data.table::setDT(DBI::dbGetQuery(con,
        "select g.dbsnpid, g.allele, f.freq
         from mygwas g
         inner join r_freqs as f on g.dbsnpid = f.dbsnpid
                                 and g.allele = f.allele"))
    DBI::dbExecute(con, "drop temporary table if exists mygwas")
    cat("Retrieved", nrow(freqs), "frequencies.\n")

    wrong <- freqs[freq <= 0 | freq >= 1]
    if (nrow(wrong) > 0) {
        cat("Removing", nrow(wrong), "SNP(s) with freq outside (0, 1).\n")
        freqs <- freqs[freq > 0 & freq < 1]
    }
    freqs
}

#' Resolve missing other alleles by querying the Genoscores dbSNP table
#'
#' @param missing A `data.table` of rows lacking `oth.allele`.
#' @param con DBI connection to the Genoscores database.
#' @keywords internal
find_missing_alleles_gs <- function(missing, con) {
    cat("Querying dbSNP for", nrow(missing), "missing allele(s) ...\n")
    genoscores:::dbTemporaryTable(con, "mymissing",
                                   missing[, .(dbsnpid, allele)])
    DBI::dbExecute(con, "create index dbsnpid on mymissing (dbsnpid)")
    found <- data.table::setDT(DBI::dbGetQuery(con,
        "select d.dbsnpid, d.ref, d.alt, m.allele
         from dbsnp d
         inner join mymissing m on d.dbsnpid = m.dbsnpid
         and (d.ref = m.allele or d.alt = m.allele)"))
    DBI::dbExecute(con, "drop temporary table if exists mymissing")

    found[allele == ref, oth.allele := alt]
    found[allele == alt, oth.allele := ref]
    found[, `:=`(ref = NULL, alt = NULL)]

    unresolved <- missing[!dbsnpid %in% found$dbsnpid]
    if (nrow(unresolved) > 0) {
        cat(nrow(unresolved),
            "SNP(s) not matched in dbSNP and will be removed:\n")
        print(unresolved[, .(snp, chrom, pos, allele)])
    }
    found
}

#' Compare SNPTEST-derived allele frequencies against dbSNP frequencies
#'
#' SNPs where the absolute difference exceeds `tolerance` are removed.
#'
#' @param region A `data.table` with columns `freq` (dbSNP) and `eff.freq`
#'   (SNPTEST).
#' @param tolerance Maximum allowed absolute difference. Default `0.1`.
#' @keywords internal
compare_freqs <- function(region, tolerance = 0.1) {
    wrong <- region[abs(freq - eff.freq) > tolerance]
    if (nrow(wrong) > 0) {
        cat(nrow(wrong), "SNP(s) removed: allele frequency differs from",
            "dbSNP by more than", tolerance * 100, "%\n")
        print(wrong)
        region <- region[abs(freq - eff.freq) <= tolerance]
    }
    region[, freq    := eff.freq]
    region[, `:=`(eff.freq = NULL, oth.freq = NULL)]
    region
}

#' Create a private regional GWAS using SNPTEST
#'
#' Runs PLINK to extract the locus, then SNPTEST to generate regional
#' association summary statistics.
#'
#' @param outcome Name of the phenotype column in `samplefile`.
#' @param samplefile Full path to the SNPTEST `.sample` file.
#' @param plinkfile Full path to the PLINK genotype file (without extension).
#' @param locus A single-row `data.table` with `chrom`, `startpos`, `endpos`.
#' @param output_dir Full path to the output directory.
#'
#' @return A `data.table` of regional GWAS summary statistics.
#' @export
create_private_gwas <- function(outcome, samplefile, plinkfile,
                                 locus, output_dir) {
    output_dir <- .valid_dir(output_dir)

    bim   <- genoscores:::readbim(plinkfile)
    build <- guess.genome.build(bim)
    if (build != "hg38")
        stop("Genotypes are not in build hg38. ",
             "Lift positions before proceeding.")

    ## extract locus from PLINK file
    range_file   <- tempfile()
    locus_plink  <- tempfile()
    locus[, range := 0]
    data.table::fwrite(locus, file = range_file, col.names = FALSE, sep = "\t")
    system(sprintf(
        "plink --silent --bfile %s --extract range %s --make-bed --out %s",
        plinkfile, range_file, locus_plink))

    ## keep only individuals present in the sample file
    sample     <- data.table::fread(samplefile, select = c("ID_1", "ID_2"),
                                    header = TRUE)
    ids_file   <- tempfile()
    data.table::fwrite(sample[2:.N], file = ids_file, sep = "\t",
                       col.names = FALSE)
    system(sprintf(
        "plink --silent --bfile %s --keep %s --make-bed --out %s",
        locus_plink, ids_file, locus_plink))

    ## run SNPTEST
    out_file <- tempfile()
    log_file <- paste0(out_file, ".log")
    system(sprintf(
        paste("snptest -data %s %s -o %s -log %s",
              "-method score -frequentist 1 -pheno %s",
              "-sex_column Gender -debug -cov_all",
              "-quantile_normalise_phenotypes -renorm"),
        paste0(locus_plink, ".bed"), samplefile,
        out_file, log_file, outcome),
        ignore.stdout = TRUE)

    gwas <- read_snptest(out_file)
    gwas[, dbsnpid := genoscores:::striprs(snp)]
    gwas <- process_private_gwas(gwas, outcome)

    out_csv <- file.path(output_dir, paste0(outcome, ".private.gwas.csv"))
    data.table::fwrite(gwas, file = out_csv)
    cat("Private GWAS saved to:", out_csv, "\n")
    gwas
}

## ---- SNPTEST file reader ------------------------------------------------

#' Read a SNPTEST output file
#' @keywords internal
read_snptest <- function(gwas_file) {
    data.table::fread(
        gwas_file, header = TRUE,
        select = c("rsid", "chromosome", "position",
                   "alleleA", "alleleB",
                   "frequentist_add_beta_1", "frequentist_add_se_1",
                   "frequentist_add_pvalue",
                   "all_AA", "all_AB", "all_BB", "all_maf", "all_total"),
        col.names = c("snp", "chrom", "pos", "oth.allele", "allele",
                      "beta", "se", "pvalue",
                      "n_AA", "n_AB", "n_BB", "maf", "n_samples"),
        skip = "alternate_ids")
}

#' Post-process a raw SNPTEST GWAS
#' @keywords internal
process_private_gwas <- function(gwas, outcome) {
    if (any(!grepl("^rs", gwas$snp)))
        stop("Variants without rsids detected. ",
             "Run find.rsids() to resolve them or remove them first.")

    check_missing_cols(names(gwas),
                       c("snp", "dbsnpid", "chrom", "pos", "oth.allele",
                         "allele", "beta", "se", "pvalue", "n_samples", "maf"))

    if (!"eaf" %in% names(gwas)) {
        check_missing_cols(names(gwas), c("n_AA", "n_AB", "n_BB"))
        gwas[, eff.freq := (2 * n_BB + n_AB) / (2 * n_samples)]
        gwas[, oth.freq := (2 * n_AA + n_AB) / (2 * n_samples)]
        n_eff <- gwas[eff.freq < 0.5, .N]
        n_oth <- gwas[oth.freq < 0.5, .N]
        gwas[, minor := if (n_eff > n_oth) eff.freq else oth.freq]
        diff_check <- gwas[abs(minor - maf) > 0.01]
        if (nrow(diff_check) > 0) {
            cat("SNPs with MAF mismatch vs SNPTEST:\n")
            print(diff_check)
        }
        gwas[, `:=`(minor = NULL, n_AA = NULL, n_AB = NULL,
                     n_BB = NULL, oth.freq = NULL)]
        data.table::setnames(gwas, "eff.freq", "eaf")
    }

    gwas <- gwas[!is.na(pvalue)]
    gwas[, gwasid     := 0L]
    gwas[, trait.name := outcome]
    check_region(gwas)
}

#' Lift genomic positions to hg38 using liftOver
#'
#' @param bim A `data.table` with columns `snp`, `chrom`, `pos`.
#' @param liftover_file Full path to the liftOver chain file.
#' @keywords internal
liftoverbim <- function(bim, liftover_file) {
    check_missing_cols(names(bim), c("snp", "chrom", "pos"))
    bim[, pos := as.integer(pos)]
    old_file    <- file.path(tempdir(), "oldbim.bed")
    lifted_file <- file.path(tempdir(), "bim_hg38.bed")
    data.table::fwrite(
        bim[, .(paste0("chr", chrom), startpos = pos - 1L, endpos = pos, snp)],
        file = old_file, sep = "\t", col.names = FALSE, scipen = 999)
    system(paste("liftOver", old_file, liftover_file, lifted_file, "/dev/null"))
    if (!file.exists(lifted_file) || file.info(lifted_file)$size == 0)
        return(NULL)
    lifted <- data.table::fread(lifted_file, select = c(4, 2),
                                 col.names = c("snp", "pos.lifted"))
    lifted[, snp       := as.character(snp)]
    lifted[, pos.lifted := as.integer(pos.lifted) + 1L]
    unlink(c(old_file, lifted_file))
    cat(nrow(lifted), "positions lifted to hg38.\n")
    lifted
}
