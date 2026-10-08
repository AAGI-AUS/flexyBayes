# predict_new_levels.R -- the posterior predictive at levels a fit did not
# see: a new trial, a new locality.
#
# For a set P of rows of `newdata` (a trial, a cell), the mean response over
# P is predicted from the posterior. Conditional on one posterior draw it is
# normal:
#
#   mean   x_P' beta + offset_P + the mean over P of every effect at a level
#          the fit saw (that level's own draw)
#   var    sum over terms of sigma2_k * sum_l (n_l / |P|)^2 over the levels l
#          the fit did not see, plus sigma2_field * the mean of the AR1 x AR1
#          correlations over P x P within each trial the fit did not see,
#          plus the mean over P of the residual variance / |P|
#
# The prediction is the mixture of these normals over the draws: its mean and
# standard deviation are exact, its quantiles are solved from the mixture
# CDF, and joint draws are simulated from the levels up.
#
# On INLA the fixed effects and the effects at seen levels are drawn jointly
# (R/inla_standardise.R, .fb_inla_config_sample()), and the variances from
# their continuous joint posterior (fb_variance_draws()), independently of the
# first. On brms the draws are brms's own,
# posterior_predict(allow_new_levels = TRUE, sample_new_levels = "gaussian").

# ---------------------------------------------------------------- #
# Arguments                                                         #
# ---------------------------------------------------------------- #

# .fb_predict_new_args() --- the arguments the new-level path takes.
#
# @noRd
# @keywords internal
.fb_predict_new_args <- function(newdata, n_draws, level) {
  if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
    stop(
      "predict(): `allow_new_levels = TRUE` predicts the rows of `newdata`, ",
      "so `newdata` must be a data frame with at least one row.",
      call. = FALSE
    )
  }
  ok <- is.numeric(n_draws) &&
    length(n_draws) == 1L &&
    is.finite(n_draws) &&
    n_draws >= 2 &&
    n_draws == round(n_draws)
  if (!ok) {
    stop(
      "predict(): `n_draws` must be a whole number of posterior draws, at ",
      "least 2.",
      call. = FALSE
    )
  }
  .fb_check_level(level)
  invisible(TRUE)
}

# .fb_predict_new_sets() --- the sets of rows whose means are predicted: one
# set per distinct combination of the `mean_by` columns, in order of first
# appearance, or one set per row when `mean_by` is NULL.
#
# @noRd
# @keywords internal
.fb_predict_new_sets <- function(newdata, mean_by) {
  n <- nrow(newdata)
  if (is.null(mean_by)) {
    return(list(
      id = seq_len(n),
      size = rep(1L, n),
      labels = data.frame(row = seq_len(n)),
      vars = character(0)
    ))
  }
  vars <- if (inherits(mean_by, "formula")) {
    if (length(mean_by) != 2L) {
      stop(
        "predict(): `mean_by` takes a one-sided formula, such as ",
        "mean_by = ~ TRIAL.",
        call. = FALSE
      )
    }
    all.vars(mean_by)
  } else if (is.character(mean_by)) {
    mean_by
  } else {
    stop(
      "predict(): `mean_by` must be NULL, a character vector of column ",
      "names, or a one-sided formula (~ TRIAL).",
      call. = FALSE
    )
  }
  missing_vars <- setdiff(vars, names(newdata))
  if (length(vars) == 0L || length(missing_vars) > 0L) {
    stop(
      "predict(): `mean_by` names column(s) not in `newdata`: ",
      paste(if (length(vars)) missing_vars else "(none)", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  key <- do.call(
    paste,
    c(lapply(vars, function(v) as.character(newdata[[v]])), sep = "\r")
  )
  id <- match(key, unique(key))
  first <- match(seq_len(max(id)), id)
  labels <- newdata[first, vars, drop = FALSE]
  rownames(labels) <- NULL
  list(id = id, size = tabulate(id), labels = labels, vars = vars)
}


# ---------------------------------------------------------------- #
# The plan: which effects a row draws, and from which variance      #
# ---------------------------------------------------------------- #

# .fb_predict_new_refuse() --- a refusal of the new-level path.
#
# @noRd
# @keywords internal
.fb_predict_new_refuse <- function(code, message, terms) {
  stop(.fb_refusal_condition(
    reason_code = code,
    message = paste0("predict(allow_new_levels = TRUE): ", message),
    terms = terms
  ))
}

# .fb_predict_new_key() --- one key per row from the named columns.
#
# @noRd
# @keywords internal
.fb_predict_new_key <- function(data, vars) {
  do.call(
    paste,
    c(lapply(vars, function(v) as.character(data[[v]])), sep = "\r")
  )
}

# .fb_predict_new_position() --- for each newdata key, the latent position a
# training row with the same key loads, or NA for a key the fit did not see.
#
# @noRd
# @keywords internal
.fb_predict_new_position <- function(train_key, train_position, new_key) {
  first <- !duplicated(train_key)
  train_position[first][match(new_key, train_key[first])]
}

# .fb_predict_new_plan() --- the random and residual structure of a fit read
# for the rows of `newdata`. Refuses a structure whose new-level draw is not
# built, and a new level with no fitted variance to draw it from. With
# `positions`, also the INLA latent position of every seen level.
#
# @noRd
# @keywords internal
.fb_predict_new_plan <- function(object, newdata, positions = FALSE) {
  ir <- .fb_fit_ir(object)
  train <- .fb_fit_data(object)
  legal <- object$data %||% train
  n <- nrow(newdata)

  weights <- .fb_ir_weights(ir) %||% object$extras$call_info$weights
  blocked <- c(
    if (.fb_weights_nonconstant(weights)) "observation weights",
    if (!is.null(.fb_known_variance_term(ir))) "a known variance per row"
  )

  random_terms <- ir$random_terms %||% list()
  for (t in random_terms) {
    type <- t$type %||% "unknown"
    built <- type %in%
      c("simple", "id", "ide", "nested", "combo") ||
      (identical(type, "at_simple") && is.null(t$level)) ||
      (identical(type, "ar1_spatial") && isTRUE(t$col_ar1))
    if (!built) {
      named <- unlist(t[c(
        "var",
        "outer",
        "inner",
        "row_var",
        "col_var",
        "vars",
        "at_var",
        "level",
        "slope_var"
      )])
      blocked <- c(
        blocked,
        paste0(type, "(", paste(named, collapse = ", "), ")")
      )
    }
  }
  for (t in ir$residual_terms %||% list()) {
    type <- t$type %||% "unknown"
    if (
      !(identical(type, "units") ||
        (identical(type, "at_units") && is.null(t$level)))
    ) {
      blocked <- c(blocked, paste0("residual ", type))
    }
  }
  if (length(blocked) > 0L) {
    .fb_predict_new_refuse(
      "predict_new_levels_structure_not_supported",
      paste0(
        "a level the fit did not see is drawn for random intercepts, ",
        "interactions such as TRIAL:REP, diag(f):g and ",
        "id(trial):ar1(row):ar1(col), over an independent residual or one ",
        "variance per section (dsum(~ units | f)). This model also carries ",
        paste(unique(blocked), collapse = ", "),
        ", for which the draw is not built. predict() without ",
        "allow_new_levels gives the population-level prediction."
      ),
      unique(blocked)
    )
  }

  needed <- unique(c(
    unlist(
      lapply(random_terms, function(t) {
        c(.fb_predict_new_term_vars(t), t$at_var)
      }),
      use.names = FALSE
    ),
    unlist(lapply(ir$residual_terms %||% list(), function(t) t$var))
  ))
  absent <- setdiff(needed, names(newdata))
  if (length(absent) > 0L) {
    stop(
      "predict(): `newdata` lacks column(s) the random or residual part of ",
      "the model reads: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  na_cols <- needed[vapply(needed, function(v) anyNA(newdata[[v]]), TRUE)]
  if (length(na_cols) > 0L) {
    stop(
      "predict(): `newdata` has missing values in column(s) ",
      paste(na_cols, collapse = ", "),
      ", which the random or residual part of the model reads.",
      call. = FALSE
    )
  }

  plan <- list(
    n = n,
    seen = list(),
    iid = list(),
    field = NULL,
    residual = NULL
  )
  random_ids <- names(object$inla$summary.random %||% list())
  seen_of <- function(tag, train_position, train_key, new_key) {
    list(
      tag = tag,
      position = .fb_predict_new_position(train_key, train_position, new_key)
    )
  }
  inla_id <- function(tag) {
    if (!tag %in% random_ids) {
      stop(
        "flexyBayes internal error: the INLA fit carries no latent field `",
        tag,
        "`.",
        call. = FALSE
      )
    }
    object$inla$summary.random[[tag]]$ID
  }
  add_iid <- function(component, new_key) {
    lev <- rep(NA_integer_, n)
    lev[!is.na(new_key)] <- match(
      new_key[!is.na(new_key)],
      unique(new_key[!is.na(new_key)])
    )
    plan$iid[[length(plan$iid) + 1L]] <<- list(
      component = component,
      level = lev
    )
  }

  for (t in random_terms) {
    type <- t$type
    if (type %in% c("simple", "id", "ide")) {
      var <- t$var
      train_key <- as.character(train[[var]])
      new_key <- as.character(newdata[[var]])
      is_new <- !new_key %in% train_key
      if (positions) {
        pos <- match(as.character(legal[[var]]), as.character(inla_id(var)))
        plan$seen[[length(plan$seen) + 1L]] <- seen_of(
          var,
          pos,
          train_key,
          ifelse(is_new, NA_character_, new_key)
        )
      }
      add_iid(paste0("sd_", var), ifelse(is_new, new_key, NA_character_))
      next
    }
    if (type %in% c("nested", "combo")) {
      vars <- .fb_nested_vars(t)
      train_key <- .fb_predict_new_key(train, vars)
      new_key <- .fb_predict_new_key(newdata, vars)
      is_new <- !new_key %in% train_key
      if (positions) {
        tag <- .fb_inla_interaction_index(vars)
        idx <- as.integer(
          interaction(legal[vars], drop = TRUE, lex.order = TRUE)
        )
        plan$seen[[length(plan$seen) + 1L]] <- seen_of(
          tag,
          match(idx, inla_id(tag)),
          train_key,
          ifelse(is_new, NA_character_, new_key)
        )
      }
      add_iid(
        paste0("sd_", paste(vars, collapse = ":")),
        ifelse(is_new, new_key, NA_character_)
      )
      next
    }
    if (identical(type, "at_simple")) {
      outer <- t$outer
      inner <- t$inner
      at_new <- as.character(newdata[[outer]])
      unseen <- setdiff(unique(at_new), as.character(train[[outer]]))
      if (length(unseen) > 0L) {
        .fb_predict_new_refuse(
          "predict_new_levels_no_variance",
          paste0(
            "diag(",
            outer,
            "):",
            inner,
            " carries one variance per level of `",
            outer,
            "`, and `newdata` has level(s) the fit did not see: ",
            paste(unseen, collapse = ", "),
            ". There is no fitted variance to draw their effects from."
          ),
          paste0("diag(", outer, "):", inner)
        )
      }
      train_key <- .fb_predict_new_key(train, c(outer, inner))
      new_key <- .fb_predict_new_key(newdata, c(outer, inner))
      is_new <- !new_key %in% train_key
      inner_pos <- as.integer(factor(legal[[inner]]))
      for (l in unique(as.character(train[[outer]]))) {
        here <- at_new == l
        if (!any(here)) {
          next
        }
        if (positions) {
          legal_l <- unique(as.character(legal[[outer]])[
            as.character(train[[outer]]) == l
          ])
          tag <- .fb_inla_at_index(outer, inner, legal_l)
          pos <- match(inner_pos, inla_id(tag))
          pos[as.character(train[[outer]]) != l] <- NA_integer_
          plan$seen[[length(plan$seen) + 1L]] <- seen_of(
            tag,
            pos,
            train_key,
            ifelse(here & !is_new, new_key, NA_character_)
          )
        }
        add_iid(
          paste0("sd_", inner, "__", outer, l),
          ifelse(here & is_new, new_key, NA_character_)
        )
      }
      next
    }
    # ar1_spatial with an AR1 column: the separable field, one per trial.
    plan$field <- .fb_predict_new_field_plan(
      t,
      newdata,
      train,
      legal,
      positions,
      inla_id
    )
    if (positions) {
      plan$seen[[length(plan$seen) + 1L]] <- plan$field$seen
    }
  }

  plan$residual <- .fb_predict_new_residual_plan(ir, newdata, train)
  plan$seen <- Filter(Negate(is.null), plan$seen)
  plan
}

# .fb_predict_new_term_vars() --- the data columns a random term reads.
#
# @noRd
# @keywords internal
.fb_predict_new_term_vars <- function(t) {
  type <- t$type %||% ""
  if (type %in% c("nested", "combo")) {
    return(.fb_nested_vars(t))
  }
  if (identical(type, "at_simple")) {
    return(c(t$outer, t$inner))
  }
  if (identical(type, "ar1_spatial")) {
    return(c(t$row_var, t$col_var))
  }
  t$var
}

# .fb_predict_new_field_plan() --- the AR1 x AR1 field read for the rows of
# `newdata`: each row's lattice row and column on the fitted lattice, the
# new trial it falls in (NA at a trial the fit saw), and, with `positions`,
# the latent node of a row at a seen trial.
#
# @noRd
# @keywords internal
.fb_predict_new_field_plan <- function(
  t,
  newdata,
  train,
  legal,
  positions,
  inla_id
) {
  label <- paste0(
    if (!is.null(t$at_var)) paste0("id(", t$at_var, "):") else "",
    "ar1(",
    t$row_var,
    "):ar1(",
    t$col_var,
    ")"
  )
  lattice <- function(v) {
    id <- as.integer(factor(legal[[v]]))
    new_id <- .fb_predict_new_position(
      as.character(train[[v]]),
      id,
      as.character(newdata[[v]])
    )
    if (anyNA(new_id)) {
      .fb_predict_new_refuse(
        "predict_new_levels_no_variance",
        paste0(
          label,
          " is fitted on the lattice of `",
          v,
          "` levels in the ",
          "data, and `newdata` has level(s) off it: ",
          paste(
            unique(as.character(newdata[[v]])[is.na(new_id)]),
            collapse = ", "
          ),
          ". The field's correlation is not defined there."
        ),
        label
      )
    }
    list(id = id, new_id = new_id, n = max(id))
  }
  row <- lattice(t$row_var)
  col <- lattice(t$col_var)

  block <- rep(NA_integer_, nrow(newdata))
  trial_pos <- rep(1L, nrow(newdata))
  train_trial <- rep(1L, nrow(train))
  if (!is.null(t$at_var)) {
    train_key <- as.character(train[[t$at_var]])
    new_key <- as.character(newdata[[t$at_var]])
    train_trial <- as.integer(factor(legal[[t$at_var]]))
    trial_pos <- .fb_predict_new_position(train_key, train_trial, new_key)
    is_new <- is.na(trial_pos)
    block[is_new] <- match(new_key[is_new], unique(new_key[is_new]))
  }

  seen <- NULL
  if (positions) {
    # The latent node of a grouped, replicated INLA field: row first, then
    # group (column), then replicate (trial).
    tag <- paste0(t$row_var, "_id")
    node <- row$new_id +
      row$n * (col$new_id - 1L) +
      row$n * col$n * (trial_pos - 1L)
    node[!is.na(block)] <- NA_integer_
    if (length(inla_id(tag)) != row$n * col$n * max(train_trial)) {
      stop(
        "flexyBayes internal error: the INLA field `",
        tag,
        "` does not ",
        "have one node per lattice cell and trial.",
        call. = FALSE
      )
    }
    seen <- list(tag = tag, position = node)
  }
  list(
    component = "sd_spatial",
    rho = c(row = "rho_row", col = "rho_col"),
    n_row = row$n,
    n_col = col$n,
    row = row$new_id,
    col = col$new_id,
    block = block,
    seen = seen
  )
}

# .fb_predict_new_residual_plan() --- each row's residual-variance component.
#
# @noRd
# @keywords internal
.fb_predict_new_residual_plan <- function(ir, newdata, train) {
  section <- Filter(
    function(t) identical(t$type, "at_units"),
    ir$residual_terms %||% list()
  )
  if (length(section) == 0L) {
    return(list(component = rep("sigma", nrow(newdata))))
  }
  var <- section[[1L]]$var
  lev <- as.character(newdata[[var]])
  unseen <- setdiff(unique(lev), as.character(train[[var]]))
  if (length(unseen) > 0L) {
    .fb_predict_new_refuse(
      "predict_new_levels_no_variance",
      paste0(
        "the residual carries one variance per level of `",
        var,
        "` (dsum(~ units | ",
        var,
        ")), and `newdata` has level(s) the fit ",
        "did not see: ",
        paste(unseen, collapse = ", "),
        ". There is no fitted residual variance for them."
      ),
      paste0("dsum(~ units | ", var, ")")
    )
  }
  list(component = paste0("sigma_", lev))
}

# .fb_predict_new_family() --- the new-level path is Gaussian on the
# identity link.
#
# @noRd
# @keywords internal
.fb_predict_new_family <- function(object) {
  fam <- family(object)
  # A sectioned residual on INLA carries one family entry per section.
  name <- unique(as.character(if (is.list(fam)) fam$family else fam))
  link <- unique(as.character(if (is.list(fam)) fam$link else "identity"))
  if (!identical(name, "gaussian") || !identical(link, "identity")) {
    .fb_predict_new_refuse(
      "predict_new_levels_family_not_supported",
      paste0(
        "the posterior predictive at new levels is built for a Gaussian ",
        "response on the identity link, where a set mean is normal given ",
        "the parameters. This fit is family '",
        name,
        "' with link '",
        link,
        "'."
      ),
      paste0(name, "(", link, ")")
    )
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------- #
# Moments of each set mean, draw by draw                            #
# ---------------------------------------------------------------- #

# .fb_predict_new_share2() --- for each set, the sum over groups of
# (rows of the set in the group / rows in the set)^2; rows with an NA group
# load no group.
#
# @noRd
# @keywords internal
.fb_predict_new_share2 <- function(sets, group) {
  out <- numeric(length(sets$size))
  keep <- which(!is.na(group))
  if (length(keep) == 0L) {
    return(out)
  }
  pair <- paste(sets$id[keep], group[keep], sep = "\r")
  u <- match(pair, unique(pair))
  count <- tabulate(u)
  set_of <- sets$id[keep][match(seq_along(count), u)]
  share <- rowsum((count / sets$size[set_of])^2, set_of)
  out[as.integer(rownames(share))] <- share[, 1L]
  out
}

# .fb_predict_new_field_counts() --- for each set, the pairs of its rows in
# one new trial binned by lattice distance, weighted 1 / |P|^2. Bin
# dr + 1 + n_row * dc, for row distance dr and column distance dc.
#
# @noRd
# @keywords internal
.fb_predict_new_field_counts <- function(sets, field) {
  n_bin <- field$n_row * field$n_col
  out <- matrix(0, length(sets$size), n_bin)
  keep <- which(!is.na(field$block))
  if (length(keep) == 0L) {
    return(out)
  }
  groups <- split(keep, paste(sets$id[keep], field$block[keep], sep = "\r"))
  for (g in groups) {
    s <- sets$id[g[1L]]
    dr <- abs(outer(field$row[g], field$row[g], "-"))
    dc <- abs(outer(field$col[g], field$col[g], "-"))
    out[s, ] <- out[s, ] +
      tabulate(dr + 1L + field$n_row * dc, n_bin) / sets$size[s]^2
  }
  out
}

# .fb_predict_new_draw_column() --- one column of the variance draws, by
# component name.
#
# @noRd
# @keywords internal
.fb_predict_new_draw_column <- function(vdraws, name) {
  if (!name %in% names(vdraws)) {
    stop(
      "flexyBayes internal error: the variance draws carry no `",
      name,
      "`; they carry ",
      paste(names(vdraws), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  as.numeric(vdraws[[name]])
}

# .fb_predict_new_moments() --- the mean and variance of every set mean,
# conditional on each draw. `eta` holds the rows' draws of the fixed part plus
# the effects at seen levels (rows x draws); `vdraws` the variance draws on
# the scale of fb_variance_draws() (standard deviations and correlations),
# one row per draw. Returns two sets x draws matrices.
#
# @noRd
# @keywords internal
.fb_predict_new_moments <- function(plan, sets, eta, vdraws) {
  n_draw <- ncol(eta)
  if (nrow(vdraws) != n_draw) {
    stop(
      "flexyBayes internal error: ",
      nrow(vdraws),
      " variance draws for ",
      n_draw,
      " draws of the linear predictor.",
      call. = FALSE
    )
  }
  mean <- rowsum(eta, sets$id, reorder = TRUE) / sets$size
  dimnames(mean) <- NULL
  variance <- matrix(0, length(sets$size), n_draw)
  add <- function(q, name) {
    if (any(q != 0)) {
      sd <- .fb_predict_new_draw_column(vdraws, name)
      variance <<- variance + outer(q, sd^2)
    }
  }
  for (b in plan$iid) {
    add(.fb_predict_new_share2(sets, b$level), b$component)
  }
  for (comp in unique(plan$residual$component)) {
    rows <- plan$residual$component == comp
    add(tabulate(sets$id[rows], length(sets$size)) / sets$size^2, comp)
  }
  f <- plan$field
  if (!is.null(f) && any(!is.na(f$block))) {
    counts <- .fb_predict_new_field_counts(sets, f)
    rho_r <- .fb_predict_new_draw_column(vdraws, f$rho[["row"]])
    rho_c <- .fb_predict_new_draw_column(vdraws, f$rho[["col"]])
    sd <- .fb_predict_new_draw_column(vdraws, f$component)
    power <- vapply(
      seq_len(n_draw),
      function(d) {
        as.vector(outer(
          rho_r[d]^(seq_len(f$n_row) - 1L),
          rho_c[d]^(seq_len(f$n_col) - 1L)
        ))
      },
      numeric(f$n_row * f$n_col)
    )
    variance <- variance + (counts %*% power) * rep(sd^2, each = nrow(counts))
  }
  list(mean = mean, variance = variance)
}

# .fb_predict_new_mixture_quantile() --- quantile p of each row's equal-weight
# mixture of normals, by bisection on the mixture CDF.
#
# @noRd
# @keywords internal
.fb_predict_new_mixture_quantile <- function(p, mean, sd) {
  lo <- apply(mean - 10 * sd, 1L, min)
  hi <- apply(mean + 10 * sd, 1L, max)
  for (i in seq_len(60L)) {
    mid <- (lo + hi) / 2
    below <- rowMeans(stats::pnorm(mid, mean, sd)) < p
    lo[below] <- mid[below]
    hi[!below] <- mid[!below]
  }
  (lo + hi) / 2
}

# .fb_predict_new_summary() --- the predictive mean, SD and central interval
# of each set mean from the conditional moments.
#
# @noRd
# @keywords internal
.fb_predict_new_summary <- function(moments, level) {
  m <- moments$mean
  v <- moments$variance
  estimate <- rowMeans(m)
  # The mixture variance: mean within-draw variance plus the spread of the
  # draw means, both over draws with equal weight.
  sd <- sqrt(rowMeans(v) + rowMeans((m - estimate)^2))
  alpha <- (1 - level) / 2
  data.frame(
    estimate = estimate,
    sd = sd,
    lower = .fb_predict_new_mixture_quantile(alpha, m, sqrt(v)),
    upper = .fb_predict_new_mixture_quantile(1 - alpha, m, sqrt(v))
  )
}

# .fb_predict_new_simulate() --- joint draws of every set mean: the
# conditional mean of each draw plus effects simulated at the new levels, the
# field of each new trial, and a residual per row, averaged over each set.
#
# @noRd
# @keywords internal
.fb_predict_new_simulate <- function(plan, sets, mean, vdraws) {
  n_draw <- ncol(mean)
  u <- matrix(0, n_draw, plan$n)
  for (b in plan$iid) {
    rows <- which(!is.na(b$level))
    if (length(rows) == 0L) {
      next
    }
    sd <- .fb_predict_new_draw_column(vdraws, b$component)
    z <- matrix(stats::rnorm(n_draw * max(b$level[rows])), n_draw) * sd
    u[, rows] <- u[, rows] + z[, b$level[rows], drop = FALSE]
  }
  sd_res <- vapply(
    plan$residual$component,
    function(comp) .fb_predict_new_draw_column(vdraws, comp),
    numeric(n_draw)
  )
  u <- u + matrix(stats::rnorm(n_draw * plan$n), n_draw) * sd_res
  f <- plan$field
  if (!is.null(f) && any(!is.na(f$block))) {
    rho_r <- .fb_predict_new_draw_column(vdraws, f$rho[["row"]])
    rho_c <- .fb_predict_new_draw_column(vdraws, f$rho[["col"]])
    sd <- .fb_predict_new_draw_column(vdraws, f$component)
    for (blk in unique(stats::na.omit(f$block))) {
      rows <- which(f$block == blk)
      node <- cbind(f$row[rows], f$col[rows])
      for (d in seq_len(n_draw)) {
        chol_r <- chol(.fb_ar1_correlation(f$n_row, rho_r[d]))
        chol_c <- chol(.fb_ar1_correlation(f$n_col, rho_c[d]))
        z <- matrix(stats::rnorm(f$n_row * f$n_col), f$n_row)
        # crossprod(chol_r, z) %*% chol_c has covariance R_row (x) R_col over
        # the (row, column) nodes.
        field <- crossprod(chol_r, z) %*% chol_c
        u[d, rows] <- u[d, rows] + sd[d] * field[node]
      }
    }
  }
  t(mean) + t(rowsum(t(u), sets$id, reorder = TRUE) / sets$size)
}


# ---------------------------------------------------------------- #
# The two engines                                                   #
# ---------------------------------------------------------------- #

# .fb_predict_new_fixed() --- the fixed-effect design of the rows, their
# offset and their estimability.
#
# @noRd
# @keywords internal
.fb_predict_new_fixed <- function(object, newdata) {
  off <- .fb_offset_values(.fb_fit_ir(object), newdata, newdata = TRUE)
  trms <- stats::delete.response(stats::terms(formula(object)))
  bhat <- coef(object)
  X <- .fb_fixef_model_matrix(trms, newdata, names(bhat), .fb_fit_data(object))
  list(
    X = X,
    offset = if (is.null(off)) numeric(nrow(newdata)) else off,
    status = .fb_fixed_alias_status(object, X)
  )
}

# .fb_predict_new_result() --- the returned table.
#
# @noRd
# @keywords internal
.fb_predict_new_result <- function(
  tab,
  draws,
  sets,
  row_status,
  level,
  engine
) {
  out <- cbind(sets$labels, n_rows = sets$size, tab)
  if (!is.null(row_status)) {
    aliased <- tabulate(sets$id[row_status == "Aliased"], length(sets$size))
    status <- ifelse(aliased > 0L, "Aliased", "Estimable")
    out[status == "Aliased", c("estimate", "sd", "lower", "upper")] <- NA_real_
    draws[, status == "Aliased"] <- NA_real_
    out$status <- status
  }
  structure(
    out,
    class = c("fb_predict_new", "data.frame"),
    draws = draws,
    level = level,
    engine = engine,
    mean_by = sets$vars
  )
}

# .fb_predict_new_inla() --- the INLA route.
#
# @noRd
# @keywords internal
.fb_predict_new_inla <- function(object, newdata, mean_by, n_draws, level) {
  .fb_predict_new_args(newdata, n_draws, level)
  .fb_predict_new_family(object)
  n_draws <- as.integer(n_draws)
  sets <- .fb_predict_new_sets(newdata, mean_by)
  plan <- .fb_predict_new_plan(object, newdata, positions = TRUE)
  fixed <- .fb_predict_new_fixed(object, newdata)

  draws <- .fb_inla_config_sample(object, n_draws)
  beta <- draws$latent[
    paste0(rownames(object$inla$summary.fixed), ":1"),
    ,
    drop = FALSE
  ]
  eta <- fixed$X %*% beta + fixed$offset
  for (s in plan$seen) {
    rows <- which(!is.na(s$position))
    if (length(rows) > 0L) {
      eta[rows, ] <- eta[rows, , drop = FALSE] +
        draws$latent[paste0(s$tag, ":", s$position[rows]), , drop = FALSE]
    }
  }
  rm(draws)
  vdraws <- fb_variance_draws(object, n_draws)

  moments <- .fb_predict_new_moments(plan, sets, eta, vdraws)
  tab <- .fb_predict_new_summary(moments, level)
  sim <- .fb_predict_new_simulate(plan, sets, moments$mean, vdraws)
  .fb_predict_new_result(tab, sim, sets, fixed$status, level, "inla")
}

# .fb_predict_new_brms() --- the brms route: brms's own posterior predictive
# at new levels, each drawn from the group-level SDs, averaged over each set.
#
# @noRd
# @keywords internal
.fb_predict_new_brms <- function(object, newdata, mean_by, n_draws, level) {
  .fb_predict_new_args(newdata, n_draws, level)
  .fb_predict_new_family(object)
  sets <- .fb_predict_new_sets(newdata, mean_by)
  .fb_predict_new_plan(object, newdata, positions = FALSE)
  # Refuses by name a newdata without the columns an offset reads.
  .fb_offset_values(.fb_fit_ir(object), newdata, newdata = TRUE)
  .check_installed(
    "brms",
    "Package 'brms' is required for predict() on a flexybayes_brms object."
  )
  total <- brms::ndraws(object$brms)
  draws <- brms::posterior_predict(
    object$brms,
    newdata = newdata,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian",
    ndraws = min(as.integer(n_draws), total)
  )
  sim <- t(rowsum(t(draws), sets$id, reorder = TRUE) / sets$size)
  dimnames(sim) <- NULL
  alpha <- (1 - level) / 2
  tab <- data.frame(
    estimate = colMeans(sim),
    sd = apply(sim, 2L, stats::sd),
    lower = apply(sim, 2L, stats::quantile, probs = alpha, names = FALSE),
    upper = apply(sim, 2L, stats::quantile, probs = 1 - alpha, names = FALSE)
  )
  .fb_predict_new_result(tab, sim, sets, NULL, level, "brms")
}


# ---------------------------------------------------------------- #
# Printing                                                          #
# ---------------------------------------------------------------- #

#' @rdname predict.flexybayes_inla
#' @param x An `fb_predict_new` table, printed with a banner naming what its
#'   numbers are.
#' @param digits Number of significant digits for the printed table.
#' @export
print.fb_predict_new <- function(x, digits = 4L, ...) {
  by <- attr(x, "mean_by")
  cat(sprintf(
    "Posterior predictive: %s  (%g%% predictive interval)\n",
    if (length(by)) {
      paste0("mean over the rows of each ", paste(by, collapse = " x "))
    } else {
      "one row at a time"
    },
    100 * (attr(x, "level") %||% 0.95)
  ))
  cat(
    "Levels the fit did not see are drawn from their fitted variances; ",
    "the residual is included.\n",
    sep = ""
  )
  body <- as.data.frame(x)
  num <- vapply(body, is.numeric, logical(1L))
  body[num] <- lapply(body[num], signif, digits = digits)
  print(body, row.names = FALSE)
  invisible(x)
}

# .fb_predict_new_route() --- the arguments that do not combine with the
# new-level path.
#
# @noRd
# @keywords internal
.fb_predict_new_route <- function(
  allow_new_levels,
  mean_by,
  classify,
  type,
  se.fit,
  re_formula = NULL,
  summary = TRUE
) {
  if (!isTRUE(allow_new_levels) && !isFALSE(allow_new_levels)) {
    stop("predict(): `allow_new_levels` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(allow_new_levels)) {
    if (!is.null(mean_by)) {
      stop(
        "predict(): `mean_by` averages the posterior predictive, which ",
        "allow_new_levels = TRUE returns.",
        call. = FALSE
      )
    }
    return(invisible(FALSE))
  }
  clash <- c(
    if (!is.null(classify)) "classify",
    if (!identical(type, "response")) "type = \"link\"",
    if (isTRUE(se.fit)) "se.fit = TRUE",
    if (!is.null(re_formula)) "re_formula",
    if (!isTRUE(summary)) "summary = FALSE"
  )
  if (length(clash) > 0L) {
    stop(
      "predict(): allow_new_levels = TRUE returns the posterior predictive ",
      "on the response scale, with its standard deviation and its draws ",
      "(attr(, \"draws\")), and does not combine with ",
      paste(clash, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
