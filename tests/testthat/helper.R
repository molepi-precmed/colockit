## Shared test helpers for colockit

#' Build a minimal valid standalone GWAS data.table for testing
#'
#' @param n       Number of SNPs (default 20).
#' @param trait   Name to assign to the trait (default "trait_A").
#' @param chrom   Chromosome (default 1L).
#' @param start   Start position for SNP placement (default 1e6).
#' @param seed    Random seed for reproducibility (default 42L).
make_test_gwas <- function(n     = 20L,
                            trait = "trait_A",
                            chrom = 1L,
                            start = 1e6,
                            seed  = 42L) {
    set.seed(seed)
    data.table::data.table(
        snp        = paste0("rs", seq_len(n)),
        dbsnpid    = as.integer(seq_len(n)),
        chrom      = as.integer(chrom),
        pos        = as.integer(seq(start, by = 1000L, length.out = n)),
        allele     = sample(c("A", "C"), n, replace = TRUE),
        oth.allele = sample(c("G", "T"), n, replace = TRUE),
        beta       = rnorm(n, sd = 0.1),
        pvalue     = runif(n, 1e-5, 0.5),
        n_samples  = 1000L,
        freq       = runif(n, 0.05, 0.45),   # keep < 0.5 for MAF compliance
        trait.name = trait
    )
}

#' Build a minimal locus data.table for testing
#'
#' @param chrom    Chromosome (default 1L).
#' @param startpos Start position (default 999000L).
#' @param endpos   End position (default 1025000L).
make_test_locus <- function(chrom    = 1L,
                             startpos = 999000L,
                             endpos   = 1025000L) {
    data.table::data.table(
        chrom    = as.integer(chrom),
        startpos = as.integer(startpos),
        endpos   = as.integer(endpos)
    )
}

#' Build a pair of GWAS with controlled allele situations for harmonise tests
#'
#' Five SNPs:
#'   rs1-rs3 — identical alleles in src and trgt  (no flip)
#'   rs4     — src allele=A matches trgt oth.allele=A  (flip needed)
#'   rs5     — irreconcilable mismatch  (dropped)
make_harmonise_pair <- function() {
    snps <- paste0("rs", 1:5)
    base <- list(
        snp       = snps,
        dbsnpid   = 1:5L,
        chrom     = 1L,
        pos       = as.integer(seq(1e6, by = 1000L, length.out = 5L)),
        pvalue    = rep(0.01, 5L),
        n_samples = 1000L,
        freq      = rep(0.3, 5L)
    )

    src <- data.table::as.data.table(c(base, list(
        allele     = c("A", "A", "A", "A", "A"),
        oth.allele = c("G", "G", "G", "G", "G"),
        beta       = c(0.1, 0.2, 0.3, 0.4, 0.5),
        trait.name = "src_trait"
    )))

    trgt <- data.table::as.data.table(c(base, list(
        allele     = c("A", "A", "A", "G", "C"),   # rs4 flipped, rs5 mismatch
        oth.allele = c("G", "G", "G", "A", "T"),
        beta       = c(0.15, 0.25, 0.35, -0.4, 0.55),
        trait.name = "trgt_trait"
    )))

    list(src = src, trgt = trgt)
}

#' Build a complete regions object suitable for passing to colocalise()
#'
#' Returns list(src, trgt, ldmat) with n SNPs sharing the same rsids,
#' a diagonal LD matrix, and all columns required by format_region().
#'
#' @param n Number of SNPs (default 10L).
make_test_regions <- function(n = 10L) {
    set.seed(1L)
    snps    <- paste0("rs", seq_len(n))
    gwas_base <- data.table::data.table(
        snp        = snps,
        dbsnpid    = as.integer(seq_len(n)),
        chrom      = 1L,
        pos        = as.integer(seq(1e6, by = 1000L, length.out = n)),
        allele     = rep("A", n),
        oth.allele = rep("G", n),
        beta       = rnorm(n, sd = 0.2),
        pvalue     = runif(n, 1e-6, 0.1),
        n_samples  = 5000L,
        freq       = runif(n, 0.05, 0.45),
        trait.name = "trait_X"
    )
    src  <- data.table::copy(gwas_base)
    trgt <- data.table::copy(gwas_base)
    trgt[, `:=`(beta       = rnorm(n, sd = 0.2),
                 pvalue     = runif(n, 1e-6, 0.1),
                 trait.name = "trait_Y")]

    src  <- get_se_varbeta(src)
    trgt <- get_se_varbeta(trgt)

    ldmat <- diag(n)
    rownames(ldmat) <- snps
    colnames(ldmat) <- snps

    list(src = src, trgt = trgt, ldmat = ldmat)
}
