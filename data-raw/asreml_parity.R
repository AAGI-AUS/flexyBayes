# data-raw/asreml_parity.R
#
# Stored ASReml 4.2 fits for the eight model shapes of the ASReml oracle
# suite (data-raw/asreml_oracle/SPEC.md, rows 01 to 07 and 15), read by
# tests/testthat/test-asreml-parity.R and by
# vignette("flexyBayes-07-from-asreml").
#
# Run from the package root with a licensed ASReml:
#
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     Rscript --vanilla data-raw/asreml_parity.R
#
# FLEXYBAYES_ASREML_LIB may be left unset where asreml is in the default
# library. The data come from tests/testthat/helper-asreml-parity-data.R,
# so the tests fit the same rows; each frame's MD5 is stored with its fit.

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

helper <- "tests/testthat/helper-asreml-parity-data.R"
source(helper)

# Variance components as ASReml reports them, on the variance scale.
varcomp_table <- function(m) {
  vc <- summary(m)$varcomp
  data.frame(
    term = rownames(vc),
    component = vc$component,
    std.error = vc$std.error,
    bound = as.character(vc$bound),
    stringsAsFactors = FALSE
  )
}

# Fixed effects, with the reference levels ASReml sets to zero left out.
fixed_table <- function(m) {
  cf <- summary(m, coef = TRUE)$coef.fixed
  out <- data.frame(
    term = rownames(cf),
    estimate = cf[, "solution"],
    std.error = cf[, "std error"],
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$std.error) & out$std.error > 0, ]
  out <- out[rev(seq_len(nrow(out))), ]
  rownames(out) <- NULL
  out
}

vpredict_table <- function(m, form) {
  x <- vpredict(m, form)
  data.frame(
    functional = rownames(x),
    estimate = x$Estimate,
    std.error = x$SE,
    stringsAsFactors = FALSE
  )
}

record <- function(shape, call, m, extra = list()) {
  d <- .ap_data(shape)
  c(
    list(
      shape = shape,
      call = call,
      n = nrow(d),
      data_md5 = .ap_data_md5(d),
      converged = isTRUE(m$converge),
      loglik = m$loglik,
      varcomp = varcomp_table(m),
      fixed = fixed_table(m)
    ),
    extra
  )
}

fits <- list()

# --- john.alpha: iid, nested and weighted -------------------------------- #

d <- .ap_data("01_iid_intercept")
m <- asreml(yield ~ gen, random = ~ rep, data = d, trace = FALSE)
fits[["01_iid_intercept"]] <- record(
  "01_iid_intercept", "asreml(yield ~ gen, random = ~ rep)", m
)

d <- .ap_data("02_nested_stratum")
m <- asreml(yield ~ gen, random = ~ rep + rep:block, data = d, trace = FALSE)
p <- predict(m, classify = "gen", trace = FALSE)$pvals
fits[["02_nested_stratum"]] <- record(
  "02_nested_stratum", "asreml(yield ~ gen, random = ~ rep + rep:block)", m,
  list(
    predict_gen = data.frame(
      gen = as.character(p$gen),
      predicted.value = p$predicted.value,
      std.error = p$std.error,
      stringsAsFactors = FALSE
    )
  )
)

d <- .ap_data("15_weights")
m <- asreml(yield ~ gen, random = ~ rep, weights = wt, data = d,
            trace = FALSE)
fits[["15_weights"]] <- record(
  "15_weights", "asreml(yield ~ gen, random = ~ rep, weights = wt)", m
)

# --- besag.met: diagonal genotype variances, alone and with dsum() ------- #

d <- .ap_data("03_diag_gxe")
m <- asreml(yield ~ env, random = ~ diag(env):gen, data = d, trace = FALSE)
m_at <- asreml(yield ~ env, random = ~ at(env):gen, data = d, trace = FALSE)
fits[["03_diag_gxe"]] <- record(
  "03_diag_gxe", "asreml(yield ~ env, random = ~ diag(env):gen)", m,
  # at(env):gen is the same model in ASReml: one genotype variance per
  # environment, no covariance.
  list(at = list(
    call = "asreml(yield ~ env, random = ~ at(env):gen)",
    converged = isTRUE(m_at$converge),
    loglik = m_at$loglik,
    varcomp = varcomp_table(m_at)
  ))
)

d <- .ap_data("05_combined_met")
m <- asreml(yield ~ env, random = ~ diag(env):gen,
            residual = ~ dsum(~ units | env), data = d, trace = FALSE)
fits[["05_combined_met"]] <- record(
  "05_combined_met",
  "asreml(yield ~ env, random = ~ diag(env):gen, residual = ~ dsum(~ units | env))",
  m
)

# --- yan.winterwheat: one residual variance per environment -------------- #

d <- .ap_data("04_dsum_residual")
m <- asreml(yield ~ env, random = ~ gen, residual = ~ dsum(~ units | env),
            data = d, trace = FALSE)
fits[["04_dsum_residual"]] <- record(
  "04_dsum_residual",
  "asreml(yield ~ env, random = ~ gen, residual = ~ dsum(~ units | env))",
  m,
  # V1 is the genotype variance, V2 the residual variance in BH93.
  list(vpredict = vpredict_table(m, gen_share_BH93 ~ V1 / (V1 + V2)))
)

# --- besag.met, two and three counties: unstructured G ------------------- #

for (shape in c("06_us_dim2", "07_us_dim3")) {
  d <- .ap_data(shape)
  # ASReml warns each time it adjusts the US matrix to keep it positive
  # definite; the count is stored with the fit.
  n_pd_adjust <- 0L
  m <- withCallingHandlers(
    {
      m <- asreml(yield ~ env, random = ~ us(env):gen, data = d,
                  trace = FALSE, maxit = 50)
      # A second call from the first one's estimates; ASReml reports
      # convergence on the call that finishes.
      if (!isTRUE(m$converge)) {
        m <- update(m)
      }
      m
    },
    warning = function(w) {
      if (grepl("remain positive definite", conditionMessage(w))) {
        n_pd_adjust <<- n_pd_adjust + 1L
        invokeRestart("muffleWarning")
      }
    }
  )
  # Genetic correlations, from the US terms in summary() order:
  # C1:C1, C2:C1, C2:C2 and, over three counties, C3:C1, C3:C2, C3:C3.
  cors <- if (identical(shape, "06_us_dim2")) {
    list(cor_gen__envC1__envC2 ~ V2 / sqrt(V1 * V3))
  } else {
    list(
      cor_gen__envC1__envC2 ~ V2 / sqrt(V1 * V3),
      cor_gen__envC1__envC3 ~ V4 / sqrt(V1 * V6),
      cor_gen__envC2__envC3 ~ V5 / sqrt(V3 * V6)
    )
  }
  fits[[shape]] <- record(
    shape, "asreml(yield ~ env, random = ~ us(env):gen)", m,
    list(
      pd_adjustments = n_pd_adjust,
      correlations = do.call(rbind, lapply(cors, vpredict_table, m = m))
    )
  )
}

asreml_parity <- list(
  stamp = list(
    date = as.character(Sys.Date()),
    asreml = as.character(utils::packageVersion("asreml")),
    agridat = as.character(utils::packageVersion("agridat")),
    r = R.version.string,
    platform = R.version$platform,
    script_md5 = unname(tools::md5sum("data-raw/asreml_parity.R")),
    helper_md5 = unname(tools::md5sum(helper))
  ),
  fits = fits
)

saveRDS(asreml_parity, "inst/extdata/asreml_parity.rds")
cat("written inst/extdata/asreml_parity.rds\n")
for (f in fits) {
  cat(sprintf("%-18s converged %-5s loglik %.3f\n", f$shape, f$converged,
              f$loglik))
}
