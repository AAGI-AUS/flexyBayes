# data-raw/asreml_classify_random.R
#
# ASReml 4.2 predict(classify = ) on models where the classify factor, or a
# term built from it, is random. The stored predicted values and standard
# errors are the oracle for tests/testthat/test-predict-classify-random.R.
#
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     Rscript --vanilla data-raw/asreml_classify_random.R
#
# Unset, the default library is used. The data are stored with the
# predictions, so the tests fit flexyBayes to the same rows ASReml saw.
asreml_lib <- Sys.getenv("FLEXYBAYES_ASREML_LIB", unset = "")
if (nzchar(asreml_lib)) {
  if (!dir.exists(asreml_lib)) {
    stop(
      "FLEXYBAYES_ASREML_LIB is set to a directory that does not exist: ",
      asreml_lib,
      call. = FALSE
    )
  }
  .libPaths(c(asreml_lib, .libPaths()))
}
suppressMessages(library(asreml))

stamp <- list(
  date = as.character(Sys.Date()),
  asreml = as.character(utils::packageVersion("asreml")),
  agridat = as.character(utils::packageVersion("agridat")),
  r = R.version.string
)

pvals <- function(m, classify) {
  p <- predict(m, classify = classify)$pvals
  out <- as.data.frame(p)
  for (v in names(out)) {
    if (is.factor(out[[v]])) out[[v]] <- as.character(out[[v]])
  }
  out
}

varcomp <- function(m) {
  s <- summary(m)$varcomp
  stats::setNames(s$component, rownames(s))
}

# --- 1. RCB: john.alpha with replicates as complete blocks -------------- #

data(john.alpha, package = "agridat")
rcb <- droplevels(john.alpha[, c("rep", "block", "gen", "yield")])

m_rcb <- asreml(yield ~ rep, random = ~gen, data = rcb, trace = FALSE)
# The incomplete blocks enter only through a random term that names a
# factor outside classify, so ASReml ignores them in the prediction.
m_rcb_blk <- asreml(
  yield ~ rep,
  random = ~ gen + rep:block,
  data = rcb,
  trace = FALSE
)

# --- 2. MET: besag.met, 16 genotypes in six counties --------------------- #

data(besag.met, package = "agridat")
met <- besag.met[besag.met$gen %in% levels(besag.met$gen)[1:16], ]
met <- met[!is.na(met$yield), c("county", "rep", "block", "gen", "yield")]
met <- droplevels(met)

m_met <- asreml(
  yield ~ county + county:rep,
  random = ~ gen + gen:county,
  data = met,
  trace = FALSE
)

# The same model with every plot of one genotype in one county removed, so
# one gen:county cell has no data.
met_gap <- droplevels(met[!(met$gen == "G01" & met$county == "C1"), ])
m_met_gap <- asreml(
  yield ~ county + county:rep,
  random = ~ gen + gen:county,
  data = met_gap,
  trace = FALSE
)

# Genotype fixed, genotype-by-county random: ASReml's default rules
# include gen:county in predict(classify = "gen"), averaged over county.
m_met_fixed <- asreml(
  yield ~ gen + county + county:rep,
  random = ~gen:county,
  data = met,
  trace = FALSE
)

ref <- list(
  stamp = stamp,
  data = list(rcb = rcb, met = met, met_gap = met_gap),
  rcb = list(
    formula = "yield ~ rep, random = ~gen",
    varcomp = varcomp(m_rcb),
    gen = pvals(m_rcb, "gen")
  ),
  rcb_blk = list(
    formula = "yield ~ rep, random = ~gen + rep:block",
    varcomp = varcomp(m_rcb_blk),
    gen = pvals(m_rcb_blk, "gen")
  ),
  met = list(
    formula = "yield ~ county + county:rep, random = ~gen + gen:county",
    varcomp = varcomp(m_met),
    gen = pvals(m_met, "gen"),
    gen_county = pvals(m_met, "gen:county"),
    county = pvals(m_met, "county")
  ),
  met_gap = list(
    formula = "yield ~ county + county:rep, random = ~gen + gen:county",
    varcomp = varcomp(m_met_gap),
    gen = pvals(m_met_gap, "gen"),
    gen_county = pvals(m_met_gap, "gen:county")
  ),
  met_fixed = list(
    formula = "yield ~ gen + county + county:rep, random = ~gen:county",
    varcomp = varcomp(m_met_fixed),
    gen = pvals(m_met_fixed, "gen")
  )
)

out <- file.path("tests", "testthat", "fixtures", "asreml_classify_random.rds")
saveRDS(ref, out, version = 3)
message("wrote ", out)
