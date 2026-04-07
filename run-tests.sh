#!/usr/bin/env bash
## run-tests.sh — run the full colockit test suite and R CMD check locally.
##
## Usage:
##   ./run-tests.sh                  # pull latest, restore deps, test + check
##   ./run-tests.sh --no-pull        # skip git pull (use current working tree)
##   ./run-tests.sh --no-check       # skip R CMD check, run tests only
##
## Environment variables:
##   COLOCKIT_REFPLINK  path to PLINK reference panel (without extension)
##                      e.g. /opt/datastore/genome/1000G/release_2020/kg.2020.hg38.eur

set -euo pipefail

PULL=true
CHECK=true

for arg in "$@"; do
  case "$arg" in
    --no-pull)  PULL=false ;;
    --no-check) CHECK=false ;;
    *) echo "Unknown argument: $arg"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

## ── 1. Pull latest code ──────────────────────────────────────────────────────
if [ "$PULL" = true ]; then
  echo ">>> git pull"
  git pull
fi

## ── 2. Sync renv library ─────────────────────────────────────────────────────
echo ">>> renv::restore()"
Rscript -e '
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/focal/latest"))
  options(pkgbuild.has_compiler = TRUE)
  renv::restore(prompt = FALSE)
'

## ── 3. Run tests ─────────────────────────────────────────────────────────────
echo ">>> testthat"
Rscript -e '
  options(pkgbuild.has_compiler = TRUE)
  pkgload::load_all(".", quiet = TRUE)
  results <- testthat::test_local(".", reporter = "progress")
  if (any(as.data.frame(results)$failed > 0))
    stop("One or more tests failed.")
'

## ── 4. R CMD check ───────────────────────────────────────────────────────────
if [ "$CHECK" = true ]; then
  echo ">>> R CMD check"
  Rscript -e '
    options(pkgbuild.has_compiler = TRUE)
    rcmdcheck::rcmdcheck(
      args      = c("--no-manual", "--no-build-vignettes"),
      error_on  = "warning",
      check_dir = "check",
      env       = c("_R_CHECK_FORCE_SUGGESTS_" = "false")
    )
  '
fi

echo ""
echo "All done."
