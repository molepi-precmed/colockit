#' colockit: A Toolkit for Principled Colocalisation Analysis
#'
#' Provides a principled two-stage pipeline for genomic colocalisation
#' analysis with two interchangeable backends (standalone and genoscores)
#' and a growing menu of colocalisation methods.
#'
#' @section Workflow:
#' 1. **Prepare regions** — harmonise two GWAS, validate allele frequencies,
#'    and compute an LD matrix via [prepare_regions()].
#' 2. **Colocalise** — run the chosen method on the prepared regions
#'    via [colocalise()].
#'
#' @section Backends:
#' * `"standalone"` — accepts user-supplied summary statistics and a reference
#'   PLINK file; no external database required.
#' * `"genoscores"` — integrates with the Genoscores database for published
#'   and private GWAS.
#'
#' @section Supported methods:
#' * `"abf"` — coloc Approximate Bayes Factor ([coloc::coloc.abf()])
#' * `"susie"` — coloc SuSiE ([coloc::coloc.susie()])
#' * `"pwcoco"` — PWCoCo (external command-line tool)
#'
#' @references
#' Robinson JW, Hemani G, Babaei MS, Huang Y, Baird DA, Tsai EA, Chen C-Y,
#' Gaunt TR, Zheng J (2022). An efficient and robust tool for colocalisation:
#' Pair-wise Conditional and Colocalisation (PWCoCo). *bioRxiv*.
#' \doi{10.1101/2022.08.08.503158}
#'
#' @aliases colockit-package
#' @importFrom stats qt
"_PACKAGE"

## Suppress R CMD check notes for data.table's non-standard evaluation
## and other internal symbols.
utils::globalVariables(c(
    ## data.table specials
    ".", ".N", ".GRP", ".SD", ":=", "..cols",
    ## GWAS column names
    "snp", "dbsnpid", "chrom", "pos", "allele", "oth.allele",
    "beta", "pvalue", "se", "varbeta", "freq", "eaf", "maf",
    "n_samples", "gwasid", "trait.name", "flatfile",
    "n_AA", "n_AB", "n_BB", "eff.freq", "oth.freq", "minor", "diff",
    "mapped.value", "trait.variable", "g",
    "newbeta", "SNP.PP.H4", "PP.H4.abf",
    "hg38pos", "pos.lifted", "startpos", "endpos",
    "i.allele", "i.oth.allele", "ref", "alt",
    ## bigsnpr / reference BIM
    "a1", "a2", "cm",
    ## genoscores (only used in backend-genoscores.R, loaded conditionally)
    "guess.genome.build"
))
