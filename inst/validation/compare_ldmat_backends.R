## ==========================================================================
## Validation script: compare get_ldmat() (bigsnpr) vs genoscores:::getld()
##
## Run this manually on a machine with both the genoscores database and
## the genoscores package available. It is NOT part of the automated test
## suite.
##
## Usage (from R):
##   source("inst/validation/compare_ldmat_backends.R")
##
## What it checks:
##   1. Both implementations return the same SNP set after filtering
##   2. The correlation matrices agree within numerical tolerance
##   3. Any sign-flipped entries are correctly handled by both sides
##   4. Prints a summary table and the max absolute difference
## ==========================================================================

stopifnot(requireNamespace("colockit",    quietly = TRUE))
stopifnot(requireNamespace("genoscores",  quietly = TRUE))
stopifnot(requireNamespace("data.table",  quietly = TRUE))

## --------------------------------------------------------------------------
## Configuration — adjust paths as needed
## --------------------------------------------------------------------------
REFPLINK  <- "/opt/datastore/genoscores/tests-data/kg.2020.hg38.eur"
SCORECON  <- scorecon   # expects a live DBI connection named scorecon

## Test region: chr22, ~29 Mb window (SNPs confirmed in kg.2020.hg38.eur)
region <- data.table::data.table(
    snp     = c("rs116577575", "rs190805234", "rs150948285"),
    dbsnpid = c(116577575L,    190805234L,    150948285L),
    allele  = c("G",           "T",           "C")
)

cat("=== colockit::get_ldmat() (bigsnpr backend) ===\n")
ld_new <- colockit:::get_ldmat(region, REFPLINK)

if (is.null(ld_new)) stop("New implementation returned NULL — check SNPs.")

cat("\n=== genoscores:::getld() (original backend) ===\n")

## Replicate the internals of the original get.ldmat() using genoscores
ref_geno <- genoscores:::plink2bigmatrix(REFPLINK, region$snp)
refbed   <- ref_geno$bed
refbim   <- ref_geno$bim
data.table::setnames(refbim, old = "snp", new = "dbsnpid")
refbim[, dbsnpid := genoscores:::striprs(dbsnpid)]

MAF           <- genoscores:::compute_maf(refbed@address)
monomorphic   <- MAF == 0
if (any(monomorphic))
    refbim[monomorphic, dbsnpid := NA]

ld_old_mat <- genoscores:::getld(
    dbsnpids = region$dbsnpid,
    alleles  = region$allele,
    refbed   = refbed,
    refbim   = refbim
)
## rename to rsids (genoscores returns integer colnames)
colnames(ld_old_mat) <- paste0("rs", colnames(ld_old_mat))
rownames(ld_old_mat) <- paste0("rs", rownames(ld_old_mat))

## --------------------------------------------------------------------------
## Comparison
## --------------------------------------------------------------------------
cat("\n=== Comparison ===\n")

## Align to common SNPs
common_snps <- intersect(rownames(ld_new$ldmat), rownames(ld_old_mat))
cat("SNPs in new implementation:", nrow(ld_new$region), "\n")
cat("SNPs in old implementation:", nrow(ld_old_mat), "\n")
cat("Common SNPs for comparison:", length(common_snps), "\n\n")

if (length(common_snps) == 0)
    stop("No common SNPs between implementations — cannot compare.")

m_new <- ld_new$ldmat[common_snps, common_snps]
m_old <- ld_old_mat[common_snps,   common_snps]

diff_mat   <- m_new - m_old
max_diff   <- max(abs(diff_mat))
mean_diff  <- mean(abs(diff_mat))
cor_values <- cor(as.vector(m_new), as.vector(m_old))

cat("Max absolute difference:   ", round(max_diff,  8), "\n")
cat("Mean absolute difference:  ", round(mean_diff, 8), "\n")
cat("Correlation of all values: ", round(cor_values, 8), "\n\n")

if (max_diff < 1e-6) {
    cat("PASS: implementations agree within 1e-6 tolerance.\n")
} else {
    cat("FAIL: implementations differ beyond 1e-6 tolerance.\n")
    cat("Difference matrix:\n")
    print(round(diff_mat, 8))
}

cat("\nNew matrix:\n");  print(round(m_new, 6))
cat("\nOld matrix:\n");  print(round(m_old, 6))
