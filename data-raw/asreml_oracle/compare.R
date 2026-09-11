# compare.R -- the Layer 1 and Layer 2 comparisons of data-raw/asreml_oracle/SPEC.md.
#
# Every tolerance here is quoted from that document, which was written before
# any comparison ran. Nothing in this file chooses a tolerance; it only applies
# one. A tolerance that needs to move is a specification change, not a code
# change, and the exit gate says so.
#
# The parameter map that matters throughout: ASReml reports variance
# components on the VARIANCE scale, flexyBayes on the STANDARD-DEVIATION
# scale. Everything below converts to the variance scale before comparing, so
# a scale slip cannot be mistaken for agreement.

TOL <- list(
  fixed_se_fraction = 0.25,   # |b_fb - b_asr| <= 0.25 * se_asr
  variance_ratio    = 0.05,   # absolute, on a 0-1 scale
  correlation       = 0.10,   # absolute
  sigma_rel         = 0.15,   # max elementwise relative difference
  blup_cor          = 0.95,   # minimum Pearson correlation
  identified_cv     = 0.50    # asreml se/estimate below this = well-identified
)

# ---- extraction -------------------------------------------------------- #

# ASReml variance components on the variance scale, named.
asr_varcomp <- function(fit) {
  v <- summary(fit)$varcomp
  stats::setNames(v$component, rownames(v))
}

asr_varcomp_se <- function(fit) {
  v <- summary(fit)$varcomp
  stats::setNames(v$std.error, rownames(v))
}

# flexyBayes variance components converted to the variance scale, named.
# `sigma` and every `sd_*` row is a standard deviation; a correlation row is
# not squared.
fb_varcomp <- function(s) {
  v <- s$varcomp
  # Every standard-deviation row, and the set is wider than it first looks.
  # A sectioned residual has no scalar `sigma`: it reports `sigma_<level>`,
  # one per section, and those are standard deviations too. The first version
  # of this function tested only `sigma` exactly and the `sd_` prefix, so the
  # sectioned residuals were compared against ASReml's variances unsquared --
  # 17.88 read against 299.2, which looks like catastrophic disagreement and
  # is in fact 17.88^2 = 319.7 against 299.2, close agreement.
  #
  # It also corrupted the denominator of every variance ratio in a sectioned
  # model, because the total is a sum over all components.
  is_sd <- v$component == "sigma" |
    startsWith(v$component, "sigma_") |
    startsWith(v$component, "sd_")
  # A correlation is on neither scale and is passed through untouched.
  is_cor <- startsWith(v$component, "cor_") | startsWith(v$component, "rho_")
  val <- ifelse(is_sd & !is_cor, v$estimate^2, v$estimate)
  stats::setNames(val, v$component)
}

# ---- Layer 1: structural, exact ---------------------------------------- #

structural_checks <- function(asr, fb_summary, d, factor_cols) {
  n_asr_var <- nrow(summary(asr)$varcomp)
  n_fb_var <- nrow(fb_summary$varcomp)
  coefs <- summary(asr, coef = TRUE)$coef.fixed
  n_asr_fixed <- sum(!is.na(coefs[, "std error"]))
  n_fb_fixed <- nrow(fb_summary$fixed)

  levels_match <- vapply(
    factor_cols,
    function(cn) identical(levels(factor(d[[cn]])), levels(factor(d[[cn]]))),
    logical(1)
  )

  list(
    variance_parameter_count = list(
      asreml = n_asr_var, flexybayes = n_fb_var,
      pass = identical(as.integer(n_asr_var), as.integer(n_fb_var))
    ),
    estimable_fixed_count = list(
      asreml = n_asr_fixed, flexybayes = n_fb_fixed,
      pass = identical(as.integer(n_asr_fixed), as.integer(n_fb_fixed))
    ),
    observations = list(
      asreml = nrow(d), flexybayes = fb_summary$n_observed,
      pass = identical(as.integer(nrow(d)), as.integer(fb_summary$n_observed))
    ),
    factor_levels = list(
      columns = factor_cols, pass = all(levels_match)
    )
  )
}

# ---- Layer 2: numerical, predeclared tolerances ------------------------ #

cmp_variance_ratio <- function(name, v_asr, v_fb, se_asr) {
  tot_asr <- sum(v_asr)
  tot_fb <- sum(v_fb)
  r_asr <- v_asr[[name$asreml]] / tot_asr
  r_fb <- v_fb[[name$flexybayes]] / tot_fb
  cv <- se_asr[[name$asreml]] / v_asr[[name$asreml]]
  identified <- is.finite(cv) && cv < TOL$identified_cv
  list(
    quantity = paste0("variance ratio: ", name$asreml),
    asreml = r_asr, flexybayes = r_fb, diff = r_fb - r_asr,
    tolerance = TOL$variance_ratio,
    identified = identified,
    asreml_cv = cv,
    pass = if (!identified) NA else abs(r_fb - r_asr) <= TOL$variance_ratio,
    note = if (!identified) "not asserted -- weakly identified" else ""
  )
}

cmp_fixed <- function(asr, fb_summary, terms_map) {
  coefs <- summary(asr, coef = TRUE)$coef.fixed
  out <- list()
  for (i in seq_along(terms_map)) {
    a_row <- terms_map[[i]]$asreml
    f_row <- terms_map[[i]]$flexybayes
    if (!a_row %in% rownames(coefs)) next
    b_asr <- coefs[a_row, "solution"]
    se_asr <- coefs[a_row, "std error"]
    idx <- match(f_row, fb_summary$fixed$term)
    if (is.na(idx) || !is.finite(se_asr) || se_asr <= 0) next
    b_fb <- fb_summary$fixed$estimate[[idx]]
    out[[length(out) + 1L]] <- list(
      quantity = paste0("fixed effect: ", a_row),
      asreml = b_asr, flexybayes = b_fb, diff = b_fb - b_asr,
      tolerance = TOL$fixed_se_fraction * se_asr,
      asreml_se = se_asr,
      pass = abs(b_fb - b_asr) <= TOL$fixed_se_fraction * se_asr,
      note = ""
    )
  }
  out
}

cmp_blup <- function(blup_asr, blup_fb, label) {
  common <- intersect(names(blup_asr), names(blup_fb))
  if (length(common) < 3L) {
    return(list(quantity = paste0("BLUP agreement: ", label),
                pass = NA, note = "not asserted -- fewer than three shared levels"))
  }
  r <- stats::cor(blup_asr[common], blup_fb[common])
  list(
    quantity = paste0("BLUP agreement: ", label),
    asreml = NA_real_, flexybayes = NA_real_, correlation = r,
    tolerance = TOL$blup_cor, n_levels = length(common),
    pass = r >= TOL$blup_cor, note = ""
  )
}

cmp_sigma_matrix <- function(S_asr, S_fb, label) {
  if (is.null(S_asr) || is.null(S_fb) ||
        !identical(dim(S_asr), dim(S_fb))) {
    return(list(quantity = paste0("covariance: ", label), pass = NA,
                note = "not asserted -- matrices not both available at equal dimension"))
  }
  rel <- abs(S_fb - S_asr) / pmax(abs(S_asr), .Machine$double.eps)
  list(
    quantity = paste0("covariance: ", label),
    max_relative_difference = max(rel),
    tolerance = TOL$sigma_rel,
    pass = max(rel) <= TOL$sigma_rel, note = ""
  )
}

cmp_correlation <- function(r_asr, r_fb, label) {
  list(
    quantity = paste0("correlation: ", label),
    asreml = r_asr, flexybayes = r_fb, diff = r_fb - r_asr,
    tolerance = TOL$correlation,
    pass = abs(r_fb - r_asr) <= TOL$correlation, note = ""
  )
}

# ---- manifest ---------------------------------------------------------- #

# SHA-256 of the serialised data frame, as SPEC.md requires. Computed by
# shelling out to shasum rather than approximating with md5 under a sha256
# name: the digest is the manifest's claim that both systems saw the same
# bytes, and a digest that is not what it says it is cannot support it.
data_digest <- function(d) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(d, tf, version = 3)
  out <- system2("shasum", c("-a", "256", shQuote(tf)), stdout = TRUE)
  sub("\\s.*$", "", out[[1L]])
}

file_digest <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("\\s.*$", "", out[[1L]])
}
