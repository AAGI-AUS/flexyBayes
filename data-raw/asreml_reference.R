# data-raw/asreml_reference.R
#
# The ASReml side of vignette("flexyBayes-04-illustration"), computed once
# and stored, so the vignette builds for a reader with no ASReml licence.
#
# Run under R 4.5.2 with a licensed ASReml:
#   R_HOME=/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources \
#     $R_HOME/bin/exec/R --vanilla -f data-raw/asreml_reference.R
#
# Every number the vignette attributes to ASReml comes from this file and
# from no other source. Re-run it and the vignette's comparison moves with
# it; edit the stored numbers by hand and the provenance is gone.

# ASReml is licensed and installed outside the default library on the
# maintainer's machine. The path is read from the environment rather than
# written here: a hard-coded home directory makes this recipe unrunnable by
# anyone else and leaks a workstation path into the repository.
#
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     R --vanilla -f data-raw/asreml_reference.R
#
# Unset, the default library is used, which is correct where asreml is
# installed normally.
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
  r = R.version.string
)

wald <- function(x) {
  c(lower = x$Estimate - 1.96 * x$SE, upper = x$Estimate + 1.96 * x$SE)
}

# --- 1. john.alpha: an alpha design, whole and degraded ----------------- #

data(john.alpha, package = "agridat")

alpha_row <- function(label, d) {
  m <- asreml(yield ~ rep, random = ~ rep:block + gen, data = d,
              trace = FALSE)
  h <- vpredict(m, h2 ~ V2 / (V2 + V3))
  w <- wald(h)
  data.frame(
    design = label, n = nrow(d),
    h2 = h$Estimate, se = h$SE,
    lower = unname(w[["lower"]]), upper = unname(w[["upper"]]),
    stringsAsFactors = FALSE
  )
}

alpha <- rbind(
  alpha_row("24 genotypes, 3 reps (whole)", john.alpha),
  do.call(rbind, lapply(c(12L, 8L, 6L), function(g) {
    keep <- levels(john.alpha$gen)[seq_len(g)]
    alpha_row(paste0(g, " genotypes, 3 reps"),
              droplevels(subset(john.alpha, gen %in% keep)))
  })),
  do.call(rbind, lapply(c(0.25, 0.50, 0.65), function(frac) {
    set.seed(2026)
    d <- droplevels(john.alpha[-sample(nrow(john.alpha),
                                       round(frac * nrow(john.alpha))), ])
    alpha_row(paste0(100 * frac, "% of plots dropped at random"), d)
  }))
)

# --- 2. dasilva.maize: the factor-analytic multi-environment trial ------ #

data(dasilva.maize, package = "agridat")
d <- dasilva.maize[order(dasilva.maize$env, dasilva.maize$gen), ]

met_fit <- asreml(yield ~ env, random = ~ fa(env, 1):gen,
                  residual = ~ dsum(~ units | env), data = d,
                  trace = FALSE, maxit = 40)
stopifnot(isTRUE(met_fit$converge))

vc <- summary(met_fit)$varcomp
met_varcomp <- data.frame(
  term = rownames(vc), component = vc$component, se = vc$std.error,
  z = vc$z.ratio, stringsAsFactors = FALSE
)

# Components 1-9 are the environment-specific variances, 10-18 the
# loadings, 19-27 the per-environment residuals, in environment order.
met_functional <- function(label, form) {
  x <- vpredict(met_fit, form)
  w <- wald(x)
  data.frame(functional = label, estimate = x$Estimate, se = x$SE,
             lower = unname(w[["lower"]]), upper = unname(w[["upper"]]),
             stringsAsFactors = FALSE)
}

met_derived <- rbind(
  met_functional("r_E5_E9", r ~ V14 * V18 / sqrt((V14^2 + V5) * (V18^2 + V9))),
  met_functional("r_E3_E8", r ~ V12 * V17 / sqrt((V12^2 + V3) * (V17^2 + V8))),
  met_functional("h2_E5", h2 ~ (V14^2 + V5) / (V14^2 + V5 + V23))
)

asreml_reference <- list(
  stamp = stamp, alpha = alpha,
  met_varcomp = met_varcomp, met_derived = met_derived
)

saveRDS(asreml_reference, "inst/extdata/asreml_reference.rds")
cat("written inst/extdata/asreml_reference.rds\n")
print(stamp)
