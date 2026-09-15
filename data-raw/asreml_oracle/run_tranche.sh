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
#
# PROVENANCE (review finding A21). The first tranche built the tarball from
# the working directory and stamped it with `git rev-parse HEAD`. Those are
# two different things whenever the tree is dirty, and the tree was dirty:
# five banked fixtures carry source_sha 5e4e0ef while the code they exercised
# was what later became 9b1d93c. The content may well have been identical --
# the stored evidence gives a future reader no way to establish that, which
# is the whole point of storing it. So:
#
#   * the tree must be clean, or the run is refused by name;
#   * the tarball is built from `git archive <sha>`, never from the working
#     directory, so the tarball hash is a function of the commit;
#   * the comparison code (scenarios.R, compare.R, run_one.R) is taken from
#     that same archive and hashed into the record.
#
# Set FB_ORACLE_ALLOW_DIRTY=1 to run anyway. The record then carries
# dirty = TRUE and the tarball is still built from the commit, so the run is
# honest about being an exploration rather than evidence.
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

cd "$PKG"
SHA="$(git rev-parse HEAD)"
DIRTY_FILES="$(git status --porcelain)"
DIRTY="FALSE"
if [ -n "$DIRTY_FILES" ]; then
  if [ "${FB_ORACLE_ALLOW_DIRTY:-0}" != "1" ]; then
    echo "REFUSED: the working tree is not clean, so a tarball built from" >&2
    echo "$SHA would not be the code this tree holds. Commit, or re-run" >&2
    echo "with FB_ORACLE_ALLOW_DIRTY=1 to bank the result as exploration." >&2
    echo "$DIRTY_FILES" >&2
    exit 1
  fi
  DIRTY="TRUE"
fi

OUT="$WORK/oracle_out"; mkdir -p "$OUT"

STAMP="$(date +%Y%m%d_%H%M%S)"
R45HOME="$WORK/r45home_$STAMP"
sh "$SHIM_SRC" "$R45HOME"

# The source under test, extracted from the commit rather than copied from
# the tree. Everything downstream -- the tarball, the scenario code, the
# comparator -- comes from here. The directory is stamped rather than reused,
# so a previous extraction is never silently built on and nothing is removed.
SRC="$WORK/src_${SHA:0:12}_$STAMP"
mkdir -p "$SRC"
git archive --format=tar "$SHA" | tar -x -C "$SRC"

LIB="$WORK/oraclelib"; mkdir -p "$LIB"
cd "$SRC"
R CMD build --no-build-vignettes . >/dev/null
TARBALL="$SRC/$(ls -t flexyBayes_*.tar.gz | head -1)"
TSHA="$(shasum -a 256 "$TARBALL" | cut -d' ' -f1)"
SCEN_SHA="$(shasum -a 256 "$SRC/data-raw/asreml_oracle/scenarios.R" | cut -d' ' -f1)"
CMP_SHA="$(shasum -a 256 "$SRC/data-raw/asreml_oracle/compare.R" | cut -d' ' -f1)"
RUN_SHA="$(shasum -a 256 "$SRC/data-raw/asreml_oracle/run_one.R" | cut -d' ' -f1)"

R_PROFILE_USER=/dev/null R_LIBS_USER="$LIB:$ASREML_LIB" \
  "$R45HOME/bin/R" CMD INSTALL --library="$LIB" --no-test-load "$TARBALL" >/dev/null

SCENARIOS=("$@")
if [ ${#SCENARIOS[@]} -eq 0 ]; then
  # R_PROFILE_USER=/dev/null, exactly as the per-scenario calls below set it.
  # Without it this line reads ~/.Rprofile, which sources a config that calls
  # requireNamespace("data.table"); under the 4.5 shim that dyn.load
  # segfaults, and the run dies before a single scenario starts:
  #
  #   *** caught segfault *** address 0x0, cause 'invalid permissions'
  #   1: dyn.load(file, DLLpath = DLLpath, ...)
  #   8: requireNamespace("data.table", quietly = TRUE)
  #
  # The no-argument form is the documented default and had never been run --
  # every previous tranche named its scenarios, which takes this branch out of
  # the path entirely.
  SCENARIOS=($(R_PROFILE_USER=/dev/null "$R45HOME/bin/R" --no-echo --no-restore -e \
    "source(\"$SRC/data-raw/asreml_oracle/scenarios.R\"); cat(names(scenario_registry()))"))
fi

RUN_CMD="data-raw/asreml_oracle/run_tranche.sh <work-dir> ${SCENARIOS[*]}"

echo "source SHA    : $SHA (dirty: $DIRTY)"
echo "tarball SHA   : $TSHA"
echo "scenarios SHA : $SCEN_SHA"
for sc in "${SCENARIOS[@]}"; do
  FB_LIB="$LIB" FB_ORACLE_DIR="$SRC/data-raw/asreml_oracle" FB_OUT="$OUT" \
  FB_SCENARIO="$sc" FB_SOURCE_SHA="$SHA" FB_TARBALL_SHA="$TSHA" \
  FB_SCENARIOS_SHA="$SCEN_SHA" FB_COMPARE_SHA="$CMP_SHA" FB_RUNNER_SHA="$RUN_SHA" \
  FB_DIRTY="$DIRTY" FB_RUN_CMD="$RUN_CMD" \
  R_PROFILE_USER=/dev/null "$R45HOME/bin/R" --no-echo --no-restore \
    --file="$SRC/data-raw/asreml_oracle/run_one.R" 2>&1 | tail -1
done
echo "TRANCHE COMPLETE -- results in $OUT"
