#!/usr/bin/env bash
# check_suggests_absent.sh -- run the test suite as CI runs it, with the
# heavy Suggests absent.
#
# Why this exists. The maintainer's machine has brms and INLA installed,
# so `devtools::test()` here exercises every engine path and reports
# green. CI installs neither reliably: brms needs a C++ toolchain and
# INLA is off-CRAN. A suite that is green locally can therefore be red on
# the first push, and was -- 20 failures on 2026-08-28, every one a test
# that reached an engine without guarding on it.
#
# The trick is a symlink farm: a library holding every installed package
# EXCEPT the ones named, so nothing is uninstalled and nothing is built.
#
#   tools/check_suggests_absent.sh              # mask brms + INLA
#   tools/check_suggests_absent.sh brms         # mask brms only
#
# Exit code is the suite's: non-zero on any failure or warning.
set -euo pipefail

MASK=("${@:-brms INLA}")
[ $# -gt 0 ] && MASK=("$@") || MASK=(brms INLA)

PKG_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LIB_SRC="$(Rscript -e 'cat(.libPaths()[1])')"
SCR="$(mktemp -d)"
LIB="$SCR/lib"
mkdir -p "$LIB"

for p in "$LIB_SRC"/*/; do
  b="$(basename "$p")"
  skip=0
  for m in "${MASK[@]}"; do [ "$b" = "$m" ] && skip=1; done
  [ "$skip" -eq 0 ] && ln -sfn "$p" "$LIB/$b"
done

echo "masked: ${MASK[*]}"
echo "library: $LIB"
for m in "${MASK[@]}"; do
  printf "  %-8s visible? %s\n" "$m" \
    "$(R_LIBS_USER="$LIB" Rscript -e "cat(requireNamespace('$m', quietly=TRUE))" 2>/dev/null)"
done

R_LIBS_USER="$LIB" NOT_CRAN=true Rscript -e '
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  r <- as.data.frame(testthat::test_dir("tests/testthat",
                                        reporter = "silent",
                                        stop_on_failure = FALSE))
  cat(sprintf("SUGGESTS-ABSENT SUITE failed=%d warned=%d passed=%d skipped=%d\n",
              sum(r$failed), sum(r$warning), sum(r$passed), sum(r$skipped)))
  bad <- r[r$failed > 0 | r$error, c("file", "test")]
  if (nrow(bad)) {
    cat("failing:\n")
    print(utils::head(bad, 30), row.names = FALSE)
  }
  quit(status = as.integer(sum(r$failed) > 0 || sum(r$warning) > 0))
'
