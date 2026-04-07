# colockit

`colockit` is an R package for genomic colocalisation analysis with a clear two-stage workflow:

1. `prepare_regions()` — harmonises GWAS summary statistics within a locus, aligns effect alleles, and computes an LD matrix from a PLINK reference panel.
2. `colocalise()` — runs the chosen colocalisation method on the prepared data.

**Supported methods:**

| Method | Description |
|---|---|
| `"abf"` | coloc Approximate Bayes Factor |
| `"susie"` | coloc SuSiE (multiple causal variants) |
| `"pwcoco"` | PWCoCo — external command-line tool ([Robinson et al. 2022](https://doi.org/10.1101/2022.08.08.503158)) |

**Supported backends:**

| Backend | Description |
|---|---|
| `"standalone"` | Supply GWAS summary statistics directly — no database required |
| `"genoscores"` | Integrate with the Genoscores database for published/private GWAS |

## Installation

```r
remotes::install_github("molepi-precmed/colockit")
```

Or from a local checkout:

```r
remotes::install_local("path/to/colockit")
```

## Quick start

### Standalone backend

```r
library(colockit)
library(data.table)

locus <- data.table(
  chrom    = 22L,
  startpos = 29001000L,
  endpos   = 29011000L
)

## gwas1 and gwas2 are data.tables with columns:
## snp, dbsnpid, chrom, pos, allele, oth.allele,
## beta, pvalue, n_samples, freq, trait.name
gwas1 <- fread("trait1_summary_stats.csv")
gwas2 <- fread("trait2_summary_stats.csv")

## Stage 1 — harmonise and compute LD
regions <- prepare_regions(
  locus        = locus,
  gwas1        = gwas1,
  gwas2        = gwas2,
  output_dir   = "results/locus_22_29Mb",
  refplinkfile = "path/to/kg.2020.hg38.eur",
  backend      = "standalone"
)

## Stage 2 — colocalise
result <- colocalise(
  regions       = regions,
  locus         = locus,
  analysis_type = "abf",
  output_dir    = "results/locus_22_29Mb"
)

result$summary
```

### Genoscores backend

```r
con <- DBI::dbConnect(...)   # Genoscores database connection

regions <- prepare_regions(
  locus        = locus,
  gwas1        = 42L,         # integer gwasid
  gwas2        = 57L,
  output_dir   = "results/locus_22_29Mb",
  refplinkfile = "path/to/kg.2020.hg38.eur",
  backend      = "genoscores",
  gwas_type    = "published",
  con          = con
)
```

Private GWAS can be created with `create_private_gwas()` and used with
`gwas_type = "mixed"` or `gwas_type = "private"`.

## Output files

`prepare_regions()` writes to `output_dir`:

- `gwas1.csv`, `gwas2.csv` — harmonised summary statistics
- `ld.matrix.Rdata.gz` — LD correlation matrix
- `locus.plot.png` — preview locus plot

`colocalise()` additionally writes:

- `formatted.region.Rdata.gz` — formatted coloc input
- `coloc.abf.png` / `coloc.susie.png` — colocalisation plot coloured by PP.H4

## Development

Tests and `R CMD check` are run manually on a server with the 1000 Genomes
reference panel and PWCoCo installed:

```bash
./run-tests.sh            # pull, restore deps, run tests + R CMD check
./run-tests.sh --no-pull  # skip git pull (use current working tree)
```

## Learn more

- `vignette("getting-started")` — full walkthrough of both flows
- `vignette("simulated-colocalisation")` — reproducible example using `coloc::coloc_test_data` with known causal structure (PP.H4 ≈ 1)

## References

Robinson JW, Hemani G, Babaei MS, Huang Y, Baird DA, Tsai EA, Chen C-Y,
Gaunt TR, Zheng J (2022). An efficient and robust tool for colocalisation:
Pair-wise Conditional and Colocalisation (PWCoCo). *bioRxiv*.
<https://doi.org/10.1101/2022.08.08.503158>

## License

MIT © 2026 Andrii Iakovliev
