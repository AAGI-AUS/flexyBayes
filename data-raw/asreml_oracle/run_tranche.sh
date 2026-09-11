#!/bin/zsh
# run_tranche.sh -- run paired ASReml / flexyBayes oracle scenarios.
#
# Durable companion to run_one.R. The environment this needs is not the
# default one: ASReml and INLA are installed under R 4.5.2 on this machine
# while R's Current link points at 4.6, so the run goes through the
# user-space 4.5 home that r45_home.sh builds. That shim is ephemeral and is
# rebuilt here rather than assumed.
#
# Usage, from the package root:
#
#   data-raw/asreml_oracle/run_tranche.sh <work-dir> [scenario ...]
#
# With no scenario named, every scenario in scenarios.R is run. Results land
# as one .rds per scenario in <work-dir>/oracle_out, and the manifest fields
# are listed in SPEC.md.
set -e

WORK_IN="${1:?usage: run_tranche.sh <work-dir> [scenario ...]}"
shift
# Absolute. The R 4.5 home the shim builds resolves R_HOME against this path
# and every child process inherits it, so a relative work directory produces
# an R that cannot find its own library. Exercised: passing ./oracle_work
# fails at install with "cannot open .../library/utils/R/sysdata.rdb".
mkdir -p "$WORK_IN"
WORK="$(cd "$WORK_IN" && pwd)"
PKG="$(cd "$(dirname "$0")/../.." && pwd)"
SHIM_SRC="$PKG/../applications/00_preliminaries/tools/r45_home.sh"
ASREML_LIB="${FLEXYBAYES_ASREML_LIB:-$HOME/Library/R/arm64/4.5/library}"

OUT="$WORK/oracle_out"; mkdir -p "$OUT"

R45HOME="$WORK/r45home_$(date +%Y%m%d_%H%M%S)"
sh "$SHIM_SRC" "$R45HOME"

LIB="$WORK/oraclelib"; mkdir -p "$LIB"
cd "$PKG"
R CMD build --no-build-vignettes . >/dev/null
TARBALL="$(ls -t flexyBayes_*.tar.gz | head -1)"
TSHA="$(shasum -a 256 "$TARBALL" | cut -d' ' -f1)"
SHA="$(git rev-parse HEAD)"
R_PROFILE_USER=/dev/null R_LIBS_USER="$LIB:$ASREML_LIB" \
  "$R45HOME/bin/R" CMD INSTALL --library="$LIB" --no-test-load "$TARBALL" >/dev/null

SCENARIOS=("$@")
if [ ${#SCENARIOS[@]} -eq 0 ]; then
  SCENARIOS=($("$R45HOME/bin/R" --no-echo --no-restore -e \
    'source("data-raw/asreml_oracle/scenarios.R"); cat(names(scenario_registry()))'))
fi

echo "source SHA   : $SHA"
echo "tarball SHA  : $TSHA"
for sc in "${SCENARIOS[@]}"; do
  FB_LIB="$LIB" FB_ORACLE_DIR="$PKG/data-raw/asreml_oracle" FB_OUT="$OUT" \
  FB_SCENARIO="$sc" FB_SOURCE_SHA="$SHA" FB_TARBALL_SHA="$TSHA" \
  R_PROFILE_USER=/dev/null "$R45HOME/bin/R" --no-echo --no-restore \
    --file="$PKG/data-raw/asreml_oracle/run_one.R" 2>&1 | tail -1
done
echo "TRANCHE COMPLETE -- results in $OUT"
