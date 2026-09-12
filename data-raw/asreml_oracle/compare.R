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


# ---- Layer 3: the selection-validation comparisons ---------------------- #
#
# SPEC.md Amendment 5. A third layer rather than a relaxation of the second:
# Layer 2 asks whether the translation is faithful and answers it on a quarter
# of a standard error; Layer 3 asks whether a reader of the two systems'
# reports would draw the same scientific conclusion, and allows a posterior
# summary and a REML estimate to differ inside the width of the interval
# either one reports.

TOL_V1 <- list(
  fixed_posterior_sd = 2,      # |b_fb - b_asr| <= 2 * posterior sd
  fixed_relative     = 0.05,   # relative, where the effect is resolved
  fixed_resolved_z   = 2,      # |z_asr| at or above which "resolved" holds
  varcomp_relative   = 0.15,   # posterior mean AND median, variance scale
  blup_cor           = 0.99    # name-matched levels only
)


# .fb_inla_variance_moments() --- the posterior mean and median of every
# variance component, on the VARIANCE scale, taken from INLA's own marginals.
#
# Not read off `summary(fit)$varcomp`. That table reports a standard
# deviation, and the mean of a standard deviation is not the square root of
# the mean of a variance, so squaring it would compare a quantity neither
# system reports. The transform is applied to the marginal density and the
# moments are taken afterwards, which is the only order that gives the mean
# and the median of the variance itself.
#
# A precision marginal transforms as 1/x; a correlation hyperparameter is on
# neither scale and passes through untouched.
.fb_inla_variance_moments <- function(fit) {
  marg <- fit$inla$marginals.hyperpar
  if (!is.list(marg) || length(marg) == 0L) {
    return(list())
  }
  out <- lapply(names(marg), function(nm) {
    m <- marg[[nm]]
    is_precision <- grepl("^Precision", nm)
    tm <- if (is_precision) {
      INLA::inla.tmarginal(function(x) 1 / x, m)
    } else {
      m
    }
    list(
      hyperparameter = nm,
      scale = if (is_precision) "variance" else "as reported",
      mean = INLA::inla.emarginal(function(x) x, tm),
      median = INLA::inla.qmarginal(0.5, tm),
      q025 = INLA::inla.qmarginal(0.025, tm),
      q975 = INLA::inla.qmarginal(0.975, tm)
    )
  })
  stats::setNames(out, names(marg))
}


# .fb_asr_fixed_table() --- ASReml's estimable fixed effects, named.
.fb_asr_fixed_table <- function(asr) {
  coefs <- summary(asr, coef = TRUE)$coef.fixed
  keep <- !is.na(coefs[, "std error"]) & coefs[, "std error"] > 0
  data.frame(
    term = rownames(coefs)[keep],
    solution = unname(coefs[keep, "solution"]),
    std_error = unname(coefs[keep, "std error"]),
    stringsAsFactors = FALSE
  )
}


# cmp_fixed_v1() --- V1-F1 and V1-F2 of Amendment 5.
cmp_fixed_v1 <- function(asr, fb_summary, model_label) {
  tab <- .fb_asr_fixed_table(asr)
  out <- list()
  for (i in seq_len(nrow(tab))) {
    term <- tab$term[[i]]
    idx <- match(term, fb_summary$fixed$term)
    if (is.na(idx)) next
    b_asr <- tab$solution[[i]]
    se_asr <- tab$std_error[[i]]
    b_fb <- fb_summary$fixed$estimate[[idx]]
    sd_fb <- fb_summary$fixed$std.error[[idx]]
    diff <- b_fb - b_asr
    z_asr <- b_asr / se_asr
    resolved <- is.finite(z_asr) && abs(z_asr) >= TOL_V1$fixed_resolved_z
    out[[length(out) + 1L]] <- list(
      quantity = paste0("V1-F1 fixed effect (", model_label, "): ", term),
      asreml = b_asr, flexybayes = b_fb, diff = diff,
      posterior_sd = sd_fb, asreml_se = se_asr,
      tolerance = TOL_V1$fixed_posterior_sd * sd_fb,
      pass = abs(diff) <= TOL_V1$fixed_posterior_sd * sd_fb,
      note = ""
    )
    out[[length(out) + 1L]] <- list(
      quantity = paste0("V1-F2 fixed effect relative (", model_label, "): ",
                        term),
      asreml = b_asr, flexybayes = b_fb,
      relative = abs(diff) / abs(b_asr),
      asreml_z = z_asr,
      tolerance = TOL_V1$fixed_relative,
      pass = if (!resolved) NA else
        abs(diff) / abs(b_asr) <= TOL_V1$fixed_relative,
      note = if (!resolved)
        "not asserted -- the ASReml effect is not resolved from zero" else ""
    )
  }
  out
}


# cmp_varcomp_v1() --- V1-V1 and V1-V2 of Amendment 5.
#
# `map` is a named character vector: names are ASReml varcomp rownames, values
# are the INLA hyperparameter the same quantity arrives under. The map is
# written per scenario rather than inferred, because a wrong automatic match
# between two naming systems is the exact defect this suite exists to catch.
cmp_varcomp_v1 <- function(asr, moments, map, model_label) {
  v <- summary(asr)$varcomp
  out <- list()
  for (nm in names(map)) {
    hyper <- map[[nm]]
    if (!nm %in% rownames(v)) {
      out[[length(out) + 1L]] <- list(
        quantity = paste0("V1-V variance component (", model_label, "): ", nm),
        pass = NA,
        note = "not asserted -- absent from the ASReml varcomp table"
      )
      next
    }
    if (!hyper %in% names(moments)) {
      out[[length(out) + 1L]] <- list(
        quantity = paste0("V1-V variance component (", model_label, "): ", nm),
        pass = NA,
        note = paste0("not asserted -- no INLA hyperparameter named `", hyper,
                      "` on this fit")
      )
      next
    }
    v_asr <- v[nm, "component"]
    bound <- as.character(v[nm, "bound"])
    at_boundary <- !identical(bound, "P") || !is.finite(v_asr) || v_asr <= 0
    m <- moments[[hyper]]
    for (stat in c("mean", "median")) {
      value <- m[[stat]]
      rel <- abs(value - v_asr) / abs(v_asr)
      out[[length(out) + 1L]] <- list(
        quantity = paste0("V1-V", if (identical(stat, "mean")) "1" else "2",
                          " variance component (", model_label, "), posterior ",
                          stat, ": ", nm),
        asreml = v_asr, flexybayes = value, relative = rel,
        asreml_bound = bound, inla_hyperparameter = hyper,
        tolerance = TOL_V1$varcomp_relative,
        pass = if (at_boundary) NA else rel <= TOL_V1$varcomp_relative,
        note = if (at_boundary)
          paste0("not asserted -- the REML estimate is at a boundary (bound = ",
                 bound, ")") else ""
      )
    }
  }
  out
}


# cmp_blup_v1() --- V1-B1 of Amendment 5.
#
# Levels are matched by name and never by position. `strip` removes the
# decoration one system puts on a level label so the two vocabularies can be
# compared as the same words; where one side returns integer indices there is
# no name to match and the comparison is recorded as blocked rather than
# guessed.
cmp_blup_v1 <- function(blup_asr, blup_fb, label) {
  if (length(blup_fb) == 0L || length(blup_asr) == 0L) {
    return(list(quantity = paste0("V1-B1 BLUP agreement: ", label), pass = NA,
                note = "not asserted -- one system returned no predictions"))
  }
  common <- intersect(names(blup_asr), names(blup_fb))
  if (length(common) < 3L) {
    return(list(
      quantity = paste0("V1-B1 BLUP agreement: ", label), pass = NA,
      n_asreml = length(blup_asr), n_flexybayes = length(blup_fb),
      note = paste0("not asserted -- levels cannot be matched by name (",
                    length(common), " of ", length(blup_asr),
                    " ASReml levels found among the flexyBayes labels)")
    ))
  }
  r <- stats::cor(blup_asr[common], blup_fb[common])
  list(
    quantity = paste0("V1-B1 BLUP agreement: ", label),
    correlation = r,
    max_abs_difference = max(abs(blup_asr[common] - blup_fb[common])),
    sd_asreml_blup = stats::sd(blup_asr[common]),
    n_levels = length(common),
    tolerance = TOL_V1$blup_cor,
    pass = r >= TOL_V1$blup_cor, note = ""
  )
}


# .fb_blup_asreml() / .fb_blup_flexybayes() --- one named numeric vector each,
# with the naming decoration removed so the two can be matched by name.
.fb_blup_asreml <- function(asr, prefix) {
  ra <- asreml::summary.asreml(asr, coef = TRUE)$coef.random
  keep <- grep(paste0("^", prefix, "_"), rownames(ra))
  if (!length(keep)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  stats::setNames(ra[keep, "solution"],
                  sub(paste0("^", prefix, "_"), "", rownames(ra)[keep]))
}

.fb_blup_flexybayes <- function(fit, group) {
  rf_all <- ranef(fit)
  rf <- if (is.data.frame(rf_all)) {
    rf_all[rf_all$group == group, , drop = FALSE]
  } else {
    rf_all[[group]]
  }
  if (is.null(rf) || !nrow(rf)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  stats::setNames(rf$estimate, as.character(rf$level))
}
