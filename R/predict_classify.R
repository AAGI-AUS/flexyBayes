# predict_classify.R -- the marginal-means table an ASReml user takes to a
# meeting.
#
# `predict(asreml_fit, classify = "Variety")` is the object a breeder
# reads, and the flexyBayes equivalent is `predict(fit, classify =
# "Variety")`. Two routes build it, chosen by ASReml's default rules.
#
# When no random term is built only from the classify factors and the
# factors averaged over, the table is the fixed part alone and comes from
# the emmeans seam this package already registers (recover_data /
# emm_basis in R/emmeans_support.R). On an INLA fit that interval is the
# Gaussian approximation of the joint fixed-effect posterior, and the
# banner says so.
#
# When such a random term exists (`gen` random and classify = "gen"; or
# `gen:env` random, averaged over `env`), the table is the posterior of
# fixed part plus those random effects, averaged over the hyper-table
# cells. That route is in the second half of this file.


# ---------------------------------------------------------------- #
# Reading the classify argument                                     #
# ---------------------------------------------------------------- #

# .fb_classify_vars() --- the factors a classify request names.
#
# Accepts ASReml's own spelling ("Variety", "Variety:env") and the
# one-sided formula an R user reaches for (~ Variety, ~ Variety * env).
# The operator carries no meaning here: what a classify names is a set of
# factors whose combinations the table has one row for.
#
# @noRd
# @keywords internal
.fb_classify_vars <- function(classify) {
  if (inherits(classify, "formula")) {
    if (length(classify) != 2L) {
      stop(
        "predict(): `classify` given as a two-sided formula. It names the ",
        "factors the means table is broken down by, so it takes the ",
        "one-sided form -- classify = ~ Variety, or ~ Variety * env.",
        call. = FALSE
      )
    }
    vars <- all.vars(classify)
  } else if (is.character(classify)) {
    vars <- unlist(strsplit(classify, "[:*+]"))
    vars <- trimws(vars)
    vars <- vars[nzchar(vars)]
  } else {
    stop(
      "predict(): `classify` must be a character value (\"Variety\", ",
      "\"Variety:env\") or a one-sided formula (~ Variety), not an object ",
      "of class <", paste(class(classify), collapse = ", "), ">.",
      call. = FALSE
    )
  }
  if (length(vars) == 0L) {
    stop(
      "predict(): `classify` names no variables, so there is no table to ",
      "build.",
      call. = FALSE
    )
  }
  unique(vars)
}

# .fb_classify_column() --- one emmeans output column under our name.
#
# emmeans names the interval bounds by the degrees of freedom it used:
# `asymp.LCL` / `asymp.UCL` on the z-intervals this package's emm_basis
# asks for (df = Inf), `lower.CL` / `upper.CL` when a finite df is in
# play. Both spellings are read so the returned column names do not
# depend on which one emmeans chose.
#
# @noRd
# @keywords internal
.fb_classify_column <- function(tab, candidates) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit) == 0L) {
    return(rep(NA_real_, nrow(tab)))
  }
  as.numeric(tab[[hit[[1L]]]])
}


# ---------------------------------------------------------------- #
# The table                                                         #
# ---------------------------------------------------------------- #

# .fb_predict_classify() --- build the marginal-means table.
#
# @noRd
# @keywords internal
.fb_predict_classify <- function(object, classify, level = 0.95, ...) {
  .fb_check_level(level)
  .fb_check_classify_arguments(names(list(...)))
  .fb_check_classify_emmeans()

  vars <- .fb_classify_vars(classify)
  plan <- .fb_classify_plan(object, vars)
  if (length(plan$random) > 0L) {
    return(.fb_predict_classify_random(object, plan, level))
  }
  .fb_check_classify_random_factor(object, vars)
  spec <- stats::as.formula(paste("~", paste(vars, collapse = " * ")))
  emm <- emmeans::emmeans(object, specs = spec, level = level)
  tab <- as.data.frame(summary(emm, level = level))

  present <- vars[vars %in% names(tab)]
  out <- tab[, present, drop = FALSE]
  # C4/FS-26: on an INLA fit whose model touched a non-syntactic factor
  # level, the reference grid emmeans built (via R/emmeans_support.R's
  # recover_data(), which reads object$data) carries the legalised
  # (make.names()) labels; restore the user's own here so the table a
  # breeder reads shows what they wrote. object$level_labels is NULL
  # for every fit this legalisation never touched (every non-INLA
  # engine, and an INLA fit with no non-syntactic level), so this is a
  # no-op there.
  for (v in present) {
    if (is.factor(out[[v]])) {
      levels(out[[v]]) <- .inla_restore_level_labels(
        object$level_labels,
        v,
        levels(out[[v]])
      )
    } else if (is.character(out[[v]])) {
      out[[v]] <- .inla_restore_level_labels(object$level_labels, v, out[[v]])
    }
  }
  out$estimate <- .fb_classify_column(tab, c("emmean", "response", "rate"))
  out$std.error <- .fb_classify_column(tab, c("SE", "std.error"))
  out$conf.low <- .fb_classify_column(tab, c("asymp.LCL", "lower.CL"))
  out$conf.high <- .fb_classify_column(tab, c("asymp.UCL", "upper.CL"))
  # emmeans returns NA for a mean the fitting rows do not estimate.
  if (!is.null(object$extras$fixed_alias)) {
    out$status <- ifelse(is.na(out$estimate), "Aliased", "Estimable")
  }
  rownames(out) <- NULL

  structure(
    out,
    class = c("fb_predict_classify", "data.frame"),
    classify = vars,
    level = level,
    engine = .fb_fit_engine(object)
  )
}

# .fb_emmeans_available() --- is the estimation seam's package present?
#
# A named predicate rather than an inline requireNamespace(), so the
# refusal below is reachable from a test without uninstalling emmeans.
#
# @noRd
# @keywords internal
.fb_emmeans_available <- function() {
  requireNamespace("emmeans", quietly = TRUE)
}

# .fb_check_classify_emmeans() --- refuse by name, not by namespace.
#
# emmeans is a suggested package because most of flexyBayes does not need
# it. A bare requireNamespace() error carries no reason code and cannot
# be caught by class, so the guard raises the registered refusal and
# names what the fit already answers without it.
#
# @noRd
# @keywords internal
.fb_check_classify_emmeans <- function() {
  if (.fb_emmeans_available()) {
    return(invisible(TRUE))
  }
  stop(.fb_refusal_condition(
    reason_code = "classify_requires_emmeans",
    message = paste0(
      "predict(classify = ) builds the marginal-means table through ",
      "emmeans, which is not installed. emmeans is a suggested package ",
      "rather than a hard dependency because most of flexyBayes does not ",
      "need it. Install it with install.packages(\"emmeans\"), or read the ",
      "fixed effects directly from summary(fit)$fixed and the ",
      "random-effect predictions from coef(fit, what = \"random\")."
    ),
    family_class = "flexybayes_classify_requires_emmeans"
  ))
}

# .fb_check_classify_random_factor() --- refuse a random classify factor
# whose terms the default rules all leave out.
#
# Reached only when no random term is built from classify and the fixed
# variables alone. A factor that is in a random term all the same (block in
# `rep:block` with rep random) has no term to predict from, and emmeans
# would die on it with "No variable named <f> in the reference grid". A
# name the model does not carry at all is left to emmeans, whose message
# for that case is already about the right thing.
#
# @noRd
# @keywords internal
.fb_check_classify_random_factor <- function(object, vars) {
  grid_vars <- tryCatch(
    all.vars(.fb_fixef_terms(object)),
    error = function(e) character(0)
  )
  missing_vars <- setdiff(vars, grid_vars)
  if (length(missing_vars) == 0L) {
    return(invisible(TRUE))
  }
  terms <- .fb_fit_ir(object)$random_terms %||% list()
  term_vars <- lapply(terms, .fb_classify_term_vars)
  holds <- vapply(
    term_vars,
    function(v) any(missing_vars %in% v),
    logical(1L)
  )
  if (!any(holds)) {
    return(invisible(TRUE))
  }
  f <- intersect(missing_vars, unlist(term_vars[holds]))[[1L]]
  other <- setdiff(term_vars[holds][[1L]], c(vars, grid_vars))
  reason <- paste0(
    "every random term containing ", f, " also contains ",
    paste(other, collapse = ", "),
    ", which is not in classify or the fixed part, so ASReml's default ",
    "rules leave the term out. Use classify = \"",
    paste(term_vars[holds][[1L]], collapse = ":"),
    "\" instead"
  )
  entry <- .lookup_refusal("classify_random_factor_not_supported")
  stop(.fb_refusal_condition(
    reason_code = "classify_random_factor_not_supported",
    message = sprintf(entry$message_template, f, reason),
    factor = f
  ))
}

# .fb_classify_newdata_note() --- the two prediction paths do not mix.
#
# `newdata` predicts the rows a caller hands over; `classify` averages
# the fitted model over the levels it saw. Supplying both is a request
# for two different things, and only one of them is answered, so the
# other is not passed over in silence.
#
# @noRd
# @keywords internal
.fb_classify_newdata_note <- function(newdata) {
  if (is.null(newdata)) {
    return(invisible(NULL))
  }
  warning(
    "predict(): `newdata` is ignored on the classify path. The means ",
    "table averages the fitted model over the levels it was fitted to. ",
    "To predict rows you supply -- untested lines, a new site -- call ",
    "predict(fit, newdata = ...) without `classify`; on a brms fit, rows at ",
    "levels the fit did not see also need allow_new_levels = TRUE.",
    call. = FALSE
  )
  invisible(NULL)
}

# .fb_check_level() --- the credible level is a proportion.
#
# @noRd
# @keywords internal
.fb_check_level <- function(level) {
  ok <- is.numeric(level) && length(level) == 1L && !is.na(level) &&
    level > 0 && level < 1
  if (!ok) {
    stop(
      "predict(): `level` is the credible level as a proportion, strictly ",
      "between 0 and 1 -- 0.95 for a 95% interval.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------- #
# ASReml's default rules                                            #
# ---------------------------------------------------------------- #
#
# ASReml-R 4.2 reference manual, section 6.2 ("The prediction process"):
# the averaging set is the variables of the fixed terms not in classify;
# variables that define only random terms are ignored; every model term
# built entirely from classify and averaging variables is used; cells of
# the hyper-table are averaged with equal weights; a simple covariate is
# held at its mean. Those five rules are matched here. The arguments that
# change them are refused by name (.FB_CLASSIFY_ASREML_ARGS).

# predict.asreml() arguments that change the default rules or ask for
# output this table does not carry.
.FB_CLASSIFY_ASREML_ARGS <- c(
  "levels", "present", "ignore", "use", "except", "only", "associate",
  "average", "vcov", "sed", "parallel", "aliased", "design.points"
)

# Random-term types whose effects are one value per level combination,
# indexed by the factors the term names.
.FB_CLASSIFY_RANDOM_TYPES <- c("simple", "ide", "nested", "combo")

# .fb_check_classify_arguments() --- refuse a predict.asreml() argument.
#
# @noRd
# @keywords internal
.fb_check_classify_arguments <- function(supplied) {
  hit <- intersect(supplied, .FB_CLASSIFY_ASREML_ARGS)
  if (length(hit) == 0L) {
    return(invisible(TRUE))
  }
  entry <- .lookup_refusal("classify_argument_not_supported")
  stop(.fb_refusal_condition(
    reason_code = "classify_argument_not_supported",
    message = sprintf(entry$message_template, hit[[1L]]),
    argument = hit[[1L]]
  ))
}

# .fb_classify_term_vars() --- the variables one IR random term names.
#
# @noRd
# @keywords internal
.fb_classify_term_vars <- function(term) {
  fields <- c("var", "outer", "inner", "vars", "slope_var", "row_var",
              "col_var")
  vals <- unlist(term[intersect(fields, names(term))], use.names = FALSE)
  unique(as.character(vals[nzchar(vals)]))
}

# .fb_classify_plan() --- the averaging set and the random terms used.
#
# @noRd
# @keywords internal
.fb_classify_plan <- function(object, vars) {
  fixed_vars <- tryCatch(
    all.vars(.fb_fixef_terms(object)),
    error = function(e) character(0)
  )
  average <- setdiff(fixed_vars, vars)
  terms <- .fb_fit_ir(object)$random_terms %||% list()
  term_vars <- lapply(terms, .fb_classify_term_vars)
  used <- vapply(
    term_vars,
    function(v) length(v) > 0L && all(v %in% c(vars, average)),
    logical(1L)
  )
  list(
    vars = vars,
    average = average,
    random = terms[used],
    random_vars = term_vars[used],
    ignored = setdiff(unlist(term_vars[!used]), c(vars, average))
  )
}


# ---------------------------------------------------------------- #
# The table with random effects                                     #
# ---------------------------------------------------------------- #

# .fb_predict_classify_random() --- fixed part plus the random effects the
# rules include, averaged over the hyper-table.
#
# @noRd
# @keywords internal
.fb_predict_classify_random <- function(object, plan, level) {
  .fb_check_classify_random_model(object, plan)
  data <- .fb_fit_data(object)
  unknown <- setdiff(plan$vars, names(data))
  if (length(unknown) > 0L) {
    stop(
      "predict(): `classify` names ", unknown[[1L]], ", which is not a ",
      "variable of the fitted model.",
      call. = FALSE
    )
  }
  values <- lapply(c(plan$vars, plan$average), function(v) {
    x <- data[[v]]
    if (is.factor(x) || is.character(x)) {
      factor(levels(as.factor(x)), levels = levels(as.factor(x)))
    } else {
      mean(x, na.rm = TRUE)
    }
  })
  names(values) <- c(plan$vars, plan$average)
  cells <- expand.grid(values, KEEP.OUT.ATTRS = FALSE)
  table_rows <- expand.grid(values[plan$vars], KEEP.OUT.ATTRS = FALSE)
  # Each hyper-table cell averages into the classify row it shares
  # classify values with; expand.grid() varies the first column fastest,
  # so the classify columns repeat in blocks of nrow(table_rows).
  row_of_cell <- rep(seq_len(nrow(table_rows)), length.out = nrow(cells))
  avg <- Matrix::sparseMatrix(
    i = row_of_cell,
    j = seq_len(nrow(cells)),
    x = 1,
    dims = c(nrow(table_rows), nrow(cells))
  )
  avg <- Matrix::Diagonal(x = 1 / Matrix::rowSums(avg)) %*% avg

  bhat <- stats::coef(object)
  x_cells <- .fb_fixef_model_matrix(
    .fb_fixef_terms(object),
    cells,
    names(bhat),
    data
  )
  x_rows <- as.matrix(avg %*% x_cells)
  status <- .fb_fixed_alias_status(object, x_rows)

  summ <- switch(
    .fb_fit_engine(object),
    inla = .fb_classify_inla(object, plan, cells, x_cells, avg, level),
    brms = .fb_classify_brms(object, plan, cells, avg, level)
  )

  out <- table_rows
  for (v in plan$vars) {
    if (is.factor(out[[v]])) {
      out[[v]] <- as.character(out[[v]])
      out[[v]] <- factor(out[[v]], levels = unique(out[[v]]))
    }
  }
  out$estimate <- summ$estimate
  out$std.error <- summ$std.error
  out$conf.low <- summ$conf.low
  out$conf.high <- summ$conf.high
  if (!is.null(status)) {
    aliased <- status == "Aliased"
    out[aliased, c("estimate", "std.error", "conf.low", "conf.high")] <-
      NA_real_
    out$status <- status
  }
  rownames(out) <- NULL

  structure(
    out,
    class = c("fb_predict_classify", "data.frame"),
    classify = plan$vars,
    level = level,
    engine = .fb_fit_engine(object),
    random = vapply(plan$random_vars, paste, character(1L), collapse = ":"),
    average = plan$average,
    ignored = as.character(plan$ignored)
  )
}

# .fb_check_classify_random_model() --- refuse what the route cannot build.
#
# @noRd
# @keywords internal
.fb_check_classify_random_model <- function(object, plan) {
  ir <- .fb_fit_ir(object)
  types <- vapply(plan$random, function(t) t$type %||% "", character(1L))
  bad <- which(!types %in% .FB_CLASSIFY_RANDOM_TYPES)
  reason <- NULL
  if (length(bad) > 0L) {
    reason <- paste0(
      "the random term ", .fb_term_human_string(plan$random[[bad[[1L]]]]),
      " does not have one effect per level combination"
    )
  } else if (!identical(ir$family %||% "gaussian", "gaussian") ||
    !identical(ir$link %||% "identity", "identity")) {
    reason <- paste0(
      "classify with random terms supports only the identity link, and ",
      "this ", ir$family, " model uses a ", ir$link, " link"
    )
  } else if (length(ir$offsets %||% character(0)) > 0L) {
    reason <- "classify with random terms does not support an offset"
  } else if (!.fb_fit_engine(object) %in% c("inla", "brms")) {
    reason <- "classify with random terms needs an INLA or brms fit"
  }
  if (is.null(reason)) {
    return(invisible(TRUE))
  }
  factor_name <- intersect(plan$vars, unlist(plan$random_vars))[1L]
  entry <- .lookup_refusal("classify_random_factor_not_supported")
  stop(.fb_refusal_condition(
    reason_code = "classify_random_factor_not_supported",
    message = sprintf(entry$message_template, factor_name, reason),
    factor = factor_name
  ))
}

# .fb_classify_cell_keys() --- each cell's level of one random term, in the
# term's own variable order, as the fit labels interaction levels.
#
# @noRd
# @keywords internal
.fb_classify_cell_keys <- function(cells, term_vars) {
  do.call(
    paste,
    c(lapply(term_vars, function(v) as.character(cells[[v]])), sep = ":")
  )
}

# .fb_classify_inla() --- posterior summaries from INLA's configurations.
#
# Each predicted value is a linear function of the latent field, so its
# posterior is the mixture, over INLA's hyperparameter configurations, of
# the Gaussian approximations INLA stores. Means, SDs and quantiles are
# read off that mixture exactly, with no random draws.
#
# @noRd
# @keywords internal
.fb_classify_inla <- function(object, plan, cells, x_cells, avg, level) {
  cs <- object$inla$misc$configs
  if (is.null(cs)) {
    stop(
      "predict(classify = ) needs an INLA fit built with ",
      "control.compute = list(config = TRUE).",
      call. = FALSE
    )
  }
  ct <- cs$contents
  offset <- sum(ct$length) - cs$n
  field <- ct$start > offset
  start <- stats::setNames(ct$start[field] - offset, ct$tag[field])
  rows <- unlist(Map(
    function(tag, len) paste0(tag, ":", seq_len(len)),
    ct$tag[field],
    ct$length[field]
  ), use.names = FALSE)
  map <- .fb_inla_draw_map(object, rows)

  # Fixed effects: the latent field names them by the legalised label, and
  # coef() restores the user's label position by position.
  fixed_tags <- rownames(object$inla$summary.fixed)
  lat <- Matrix::sparseMatrix(
    i = rep(seq_len(nrow(cells)), ncol(x_cells)),
    j = rep(unname(start[fixed_tags]), each = nrow(cells)),
    x = as.vector(x_cells),
    dims = c(nrow(cells), cs$n)
  )

  fb <- object$extras$fb_terms %||% object$fb
  known <- object$extras$call_info$known_matrices
  group_labels <- .fb_inla_group_labels(fb, known)
  group_levels <- .fb_inla_group_levels(
    fb,
    object$data,
    object$level_labels,
    known
  )
  new_levels <- list()
  for (k in seq_along(plan$random)) {
    term_vars <- plan$random_vars[[k]]
    label <- paste(term_vars, collapse = ":")
    tag <- names(group_labels)[group_labels == label]
    tag <- if (length(tag) == 1L) tag else label
    tab <- object$inla$summary.random[[tag]]
    if (is.null(tab) || is.na(start[tag])) {
      stop(
        "flexyBayes internal error: no INLA field for the random term `",
        label, "`.",
        call. = FALSE
      )
    }
    spec <- group_levels[[tag]]
    if (is.null(spec)) {
      pos <- seq_len(nrow(tab))
      labs <- .inla_restore_level_labels(
        object$level_labels,
        tag,
        as.character(tab$ID)
      )
    } else {
      pos <- spec$rows
      labs <- spec$levels
    }
    keys <- .fb_classify_cell_keys(cells, term_vars)
    hit <- match(keys, labs)
    seen <- !is.na(hit)
    lat <- lat + Matrix::sparseMatrix(
      i = which(seen),
      j = start[[tag]] - 1L + pos[hit[seen]],
      x = 1,
      dims = c(nrow(cells), cs$n)
    )
    # A level combination the data never held has no node in the field;
    # its effect is a fresh draw from the term's normal distribution.
    if (!all(seen)) {
      unseen <- factor(keys[!seen])
      new_levels[[tag]] <- avg[, !seen, drop = FALSE] %*%
        Matrix::sparseMatrix(
          i = seq_along(unseen),
          j = as.integer(unseen),
          x = 1,
          dims = c(length(unseen), nlevels(unseen))
        )
    }
  }

  # Original scale: each latent value maps to a + b * x.
  lin <- avg %*% lat
  shift <- as.numeric(lin %*% map$a)
  lin_b <- lin %*% Matrix::Diagonal(x = map$b)
  tlin <- Matrix::t(lin_b)

  w <- .fb_inla_config_weights(cs)
  n_rows <- nrow(lin)
  internal <- rownames(object$inla$internal.summary.hyperpar)
  hyper_names <- rownames(object$inla$summary.hyperpar)[
    match(names(cs$config[[1L]]$theta), internal)
  ]
  hyper_h <- .fb_inla_draw_map(object, character(0), hyper_names)$h
  new_weight <- lapply(new_levels, function(m) Matrix::rowSums(m^2))
  new_hyper <- match(
    sprintf("Precision for %s", names(new_levels) %||% character(0)),
    hyper_names
  )
  if (anyNA(new_hyper)) {
    stop(
      "flexyBayes internal error: no INLA precision for the random term `",
      names(new_levels)[is.na(new_hyper)][[1L]], "`.",
      call. = FALSE
    )
  }
  means <- matrix(0, n_rows, length(w))
  vars <- matrix(0, n_rows, length(w))
  for (k in seq_along(cs$config)) {
    g <- .fb_inla_config_gaussian(cs, k)
    mu <- g$mu
    z <- Matrix::solve(
      g$chol_q,
      Matrix::solve(g$chol_q, tlin, system = "P"),
      system = "L"
    )
    v <- Matrix::colSums(z^2)
    if (!is.null(g$gain)) {
      mu <- mu - drop(g$gain %*% (g$a_mat %*% mu - cs$constr$e))
      u <- as.matrix(Matrix::crossprod(tlin, g$sigma_at))
      v <- v - rowSums((u %*% solve(g$a_mat %*% g$sigma_at)) * u)
    }
    for (m in seq_along(new_levels)) {
      j <- new_hyper[[m]]
      prec <- do.call(
        object$inla$misc$from.theta[[j]],
        list(cs$config[[k]]$theta[[j]])
      )
      v <- v + new_weight[[m]] / (prec * hyper_h[[j]])
    }
    means[, k] <- shift + as.numeric(lin_b %*% mu)
    vars[, k] <- pmax(as.numeric(v), 0)
  }
  estimate <- drop(means %*% w)
  second <- drop((vars + means^2) %*% w)
  sds <- sqrt(vars)
  quant <- function(p, i) {
    if (all(sds[i, ] == 0)) {
      return(estimate[[i]])
    }
    f <- function(q) sum(w * stats::pnorm(q, means[i, ], sds[i, ])) - p
    span <- range(means[i, ]) + c(-10, 10) * max(sds[i, ])
    stats::uniroot(f, span, tol = 1e-10)$root
  }
  alpha <- (1 - level) / 2
  list(
    estimate = estimate,
    std.error = sqrt(pmax(second - estimate^2, 0)),
    conf.low = vapply(seq_len(n_rows), function(i) quant(alpha, i), 0),
    conf.high = vapply(seq_len(n_rows), function(i) quant(1 - alpha, i), 0)
  )
}

# .fb_classify_brms() --- posterior summaries from the brms draws.
#
# brms evaluates the fixed part and the named random terms on the cells.
# A level combination the data never held (a genotype not grown in one
# environment) is drawn from its fitted normal distribution, as its
# prediction-error variance in ASReml is the full variance component.
#
# @noRd
# @keywords internal
.fb_classify_brms <- function(object, plan, cells, avg, level) {
  re_terms <- vapply(
    plan$random_vars,
    function(v) paste0("(1 | ", paste(v, collapse = ":"), ")"),
    character(1L)
  )
  re_formula <- stats::as.formula(
    paste("~", paste(re_terms, collapse = " + "))
  )
  draws <- brms::posterior_linpred(
    object$brms,
    newdata = cells,
    re_formula = re_formula,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )
  draws <- as.matrix(draws %*% Matrix::t(avg))
  alpha <- (1 - level) / 2
  list(
    estimate = colMeans(draws),
    std.error = apply(draws, 2L, stats::sd),
    conf.low = apply(draws, 2L, stats::quantile, probs = alpha, names = FALSE),
    conf.high = apply(
      draws, 2L, stats::quantile, probs = 1 - alpha, names = FALSE
    )
  )
}

# .populate_refusal_registry_classify_random() --- the 0.16.0 classify code.
#
# @noRd
# @keywords internal
.populate_refusal_registry_classify_random <- function() {
  .register_refusal(
    reason_code = "classify_argument_not_supported",
    description = paste0(
      "predict(classify = ) was given an argument of ASReml's predict() ",
      "that changes the default prediction rules or asks for output the ",
      "table does not carry."
    ),
    message_template = paste0(
      "predict(classify = ) does not support ASReml's `%1$s` argument. ",
      "Predictions follow ASReml's defaults: equal-weight averaging over ",
      "the fixed-term variables not in classify, variables that appear ",
      "only in random terms ignored, and every random term whose variables ",
      "are all in classify or averaged over included. For other cells or ",
      "weights, use predict(fit, newdata = ) and average the rows yourself."
    ),
    registered_in_adr = "PLAN 0.16.0 C3",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
}


# ---------------------------------------------------------------- #
# Printing                                                          #
# ---------------------------------------------------------------- #

#' Print a classify means table
#'
#' Prints the table with a banner naming what the numbers are. The
#' estimate is a posterior mean, and the table is not a
#' standard-error-of-difference table -- ASReml's `predict()` prints one
#' beside the means and this does not, so the banner says so rather than
#' leaving the absence to be noticed.
#'
#' When the table includes random effects, the banner names the random
#' terms included, the variables averaged over and the variables left
#' out, as ASReml's own notes do. Otherwise, on an INLA fit, the banner
#' says the interval comes from the Gaussian approximation of the joint
#' fixed-effect posterior, which is what the estimation basis supplies to
#' emmeans, and not from INLA's own marginal densities.
#'
#' @param x An `fb_predict_classify` table, as returned by
#'   `predict(fit, classify = )`.
#' @param digits Number of significant digits for the printed table.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the table it prints.
#' @export
print.fb_predict_classify <- function(x, digits = 4L, ...) {
  level <- attr(x, "level") %||% 0.95
  random <- attr(x, "random")
  cat(sprintf(
    "Predicted %s: %s  (%g%% credible interval)\n",
    if (is.null(random)) "means" else "values",
    paste(attr(x, "classify") %||% character(0), collapse = " x "),
    100 * level
  ))
  if (is.null(random)) {
    cat(
      "Estimate = posterior mean of the marginal mean. Not an SED table.\n"
    )
    if (identical(attr(x, "engine"), "inla")) {
      cat(
        "INLA: intervals from the Gaussian approximation of the fixed ",
        "effects.\n",
        sep = ""
      )
    }
  } else {
    cat(
      "Estimate = posterior mean of the fixed part plus ",
      paste(random, collapse = ", "),
      ". Not an SED table.\n",
      sep = ""
    )
    note <- function(label, v) {
      if (length(v) > 0L) {
        cat(label, paste(v, collapse = ", "), "\n", sep = "")
      }
    }
    note("Averaged over: ", attr(x, "average"))
    note("Ignored: ", attr(x, "ignored"))
  }
  body <- as.data.frame(x)
  num <- vapply(body, is.numeric, logical(1L))
  body[num] <- lapply(body[num], round, digits = digits)
  print(body, row.names = FALSE)
  invisible(x)
}
