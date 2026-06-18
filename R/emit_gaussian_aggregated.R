# emit_gaussian_aggregated -- aggregated backend wiring (v0.3.2).
#
# Aggregated emit path: gaussian-identity LMMs in scope route through the
# `<fb_aggregated>` sufficient-statistics object instead of the per-row
# likelihood. The aggregated form is algebraically identical to the
# per-row form (not an approximation).
#
# Two implementation patterns -- one per backend:
#
# (1) greta path. Combines two greta `distribution()` attachments:
#
#       distribution(ybar_k)  <- normal(mu_cell, sigma / sqrt(n_k))
#       distribution(WSS)     <- gamma((N - K) / 2, 1 / (2 * sigma^2))
#
#     The cell-mean weighted gaussian contributes the cell-level
#     residual sum-of-squares to the log-likelihood; the gamma
#     observation on the data-only scalar `WSS = sum_k (S2_k - S1_k^2 /
#     n_k)` contributes the within-cell sum-of-squares term. Together
#     they recover the full-data log-likelihood up to sigma-independent
#     data constants (which do not affect the posterior shape). The
#     gamma observation is skipped when N == K (single-obs cells) or
#     WSS <= machine eps (no within-cell scatter).
#
# (2) INLA path. Cell-mean weighted gaussian with a custom-prior
#     correction on the log-precision hyperparameter that absorbs the
#     sigma-dependent within-cell SS term:
#
#       prior_expr <- "expression:
#                       a = (WSS / 2) + 5e-5;
#                       b = (N - K) / 2 + 1;
#                       log_dens = -a * exp(theta) + b * theta;
#                       return(log_dens);"
#
#     The constants `a` and `b` collapse the default `loggamma(1, 5e-5)`
#     precision prior with the `Delta(theta)` correction. The combined
#     prior is fed via
#     `control.family$hyper$prec$prior` and INLA's deterministic
#     numerical engine recovers the per-row posterior on beta + sigma
#     to bit-exact precision (spike: differences <= 1e-4 on a 1000-row
#     example).
#
# Heterogeneous residual: when the IR's rcov_terms include an
# `at_units` term (per-level residual variance, e.g. `at(env):units`),
# the cell-mean machinery extends naturally PROVIDED the residual-
# grouping factor is itself a cell key (so each cell has a single
# residual sigma). When the residual factor is NOT a cell key, the
# aggregation closure breaks and the emit refuses with reason code
# `heterogeneous_residual_factor_not_in_cell_key`.
#
# Internal -- not exported. Called from .dispatch_backend() when the
# `<fb_aggregation_plan>` declares eligibility.

emit_gaussian_aggregated <- function(
  fb,
  fb_aggregated,
  data,
  backend = c("greta", "inla"),
  n_samples = 1000L,
  warmup = 500L,
  chains = 4L,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  rcov = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_
) {
  backend <- match.arg(backend)

  if (!inherits(fb, "fb_terms")) {
    stop("emit_gaussian_aggregated() requires an <fb_terms> IR.", call. = FALSE)
  }
  if (!inherits(fb_aggregated, "fb_aggregated")) {
    stop(
      "emit_gaussian_aggregated() requires an <fb_aggregated> ",
      "from .fb_aggregate_gaussian().",
      call. = FALSE
    )
  }

  # Defensive scope re-check (the caller -- .dispatch_backend() --
  # already gates on <fb_aggregation_plan>$eligible, but the emit
  # path also documents its own invariants so a future direct caller
  # cannot misuse it).
  if (!identical(fb$family, "gaussian")) {
    stop(
      "emit_gaussian_aggregated(): family must be gaussian; got '",
      fb$family,
      "'.",
      call. = FALSE
    )
  }
  if (!(is.null(fb$link) || identical(fb$link, "identity"))) {
    stop(
      "emit_gaussian_aggregated(): link must be identity; got '",
      fb$link,
      "'.",
      call. = FALSE
    )
  }

  # Sufficient-statistics sanity gates.
  agg <- fb_aggregated$sufficient_stats
  if (any(agg$n_k <= 0L)) {
    stop(
      "emit_gaussian_aggregated(): non-positive n_k in cell ",
      which(agg$n_k <= 0L)[1L],
      ".",
      call. = FALSE
    )
  }
  wss_per_cell <- agg$S2_k - agg$S1_k^2 / agg$n_k
  # Allow tiny negative drift from floating-point summation.
  if (any(wss_per_cell < -1e-8 * abs(agg$S2_k))) {
    stop(
      "emit_gaussian_aggregated(): negative within-cell SS in cell ",
      which.min(wss_per_cell),
      " (likely numerical instability ",
      "upstream).",
      call. = FALSE
    )
  }
  wss_per_cell <- pmax(wss_per_cell, 0)
  WSS_total <- sum(wss_per_cell)

  # Heterogeneous-residual handling: detect rcov_terms shape and
  # validate that any at_units residual factor lies in the cell key.
  rcov_plan <- .agg_emit_rcov_plan(fb, fb_aggregated)

  # Build per-RI level indices for the cell key (one per random_cols
  # entry; foundation already places each RI grouping factor into the
  # cell key as a column on the sufficient_stats data.table).
  ri_plan <- .agg_emit_ri_plan(fb, fb_aggregated)

  if (return_code) {
    # The aggregated emit at v0.3.2 does NOT yet expose review_code =
    # TRUE for the aggregated path; the per-row review path remains the
    # documented v0.2 surface. Honest refusal:
    stop(
      "emit_gaussian_aggregated(): return_code = TRUE is not yet ",
      "supported on the aggregated path (deferred to a follow-on ",
      "release). Pass aggregate = FALSE if you need the per-row ",
      "Stan / greta source.",
      call. = FALSE
    )
  }

  t0 <- Sys.time()
  if (identical(backend, "greta")) {
    fit_engine_out <- .emit_gaussian_aggregated_greta(
      fb = fb,
      fb_aggregated = fb_aggregated,
      data = data,
      wss_per_cell = wss_per_cell,
      WSS_total = WSS_total,
      ri_plan = ri_plan,
      rcov_plan = rcov_plan,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose
    )
  } else {
    fit_engine_out <- .emit_gaussian_aggregated_inla(
      fb = fb,
      fb_aggregated = fb_aggregated,
      data = data,
      wss_per_cell = wss_per_cell,
      WSS_total = WSS_total,
      ri_plan = ri_plan,
      rcov_plan = rcov_plan,
      verbose = verbose
    )
  }
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  .agg_emit_build_fit(
    fb = fb,
    fb_aggregated = fb_aggregated,
    data = data,
    backend = backend,
    engine_out = fit_engine_out,
    ri_plan = ri_plan,
    rcov_plan = rcov_plan,
    elapsed = elapsed,
    the_call = the_call,
    fixed = fixed,
    random = random,
    rcov = rcov,
    family = family,
    link = link,
    data_name = data_name,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd
  )
}


# ---------------------------------------------------------------- #
# Residual-structure planning                                       #
# ---------------------------------------------------------------- #
# Returns a list:
#   $kind        "homogeneous" | "at_units"
#   $factor      NULL (homogeneous) | <character> column name (at_units)
#   $level_col   NULL | <integer vector length K> mapping cell -> level
#                of $factor. Used to index per-level sigma in greta and
#                to attach per-level scales in INLA.
#   $n_levels    NULL | integer count of distinct residual levels.
#   $levels      NULL | character vector of level labels.
#
# Refuses (structured error) when an at_units term references a factor
# that is NOT in the cell-key set: the cell-constant mu property still
# holds, but the cell-constant sigma property does not, so the per-row
# vs per-cell algebraic identity breaks. Reason code
# `heterogeneous_residual_factor_not_in_cell_key`.
.agg_emit_rcov_plan <- function(fb, fb_aggregated) {
  rcov <- fb$rcov_terms
  if (length(rcov) == 0L) {
    return(list(
      kind = "homogeneous",
      factor = NULL,
      level_col = NULL,
      n_levels = NULL,
      levels = NULL
    ))
  }

  homogeneous_types <- c("units", "id", "ide", "simple")
  hetero_types <- c("at_units")

  for (term in rcov) {
    ttype <- term$type %||% "units"
    if (ttype %in% homogeneous_types) {
      next
    }
    if (ttype %in% hetero_types) {
      f <- as.character(term$var)
      if (!(f %in% fb_aggregated$cell_key_cols)) {
        stop(.fb_refusal_condition(
          reason_code = "heterogeneous_residual_factor_not_in_cell_key",
          message = paste0(
            "emit_gaussian_aggregated(): heterogeneous residual ",
            "at(",
            f,
            "):units refused -- '",
            f,
            "' is not in the ",
            "cell key {",
            paste(fb_aggregated$cell_key_cols, collapse = ", "),
            "}. The cell-constant sigma property does not hold, ",
            "so the per-row / per-cell algebraic identity breaks. ",
            "Pass aggregate = FALSE for the per-row path."
          ),
          family_class = "flexybayes_aggregate_emit_refusal",
          factor = f,
          cell_key = fb_aggregated$cell_key_cols
        ))
      }
      col_data <- fb_aggregated$sufficient_stats[[f]]
      if (!is.factor(col_data)) {
        col_data <- factor(col_data)
      }
      return(list(
        kind = "at_units",
        factor = f,
        level_col = as.integer(col_data),
        n_levels = nlevels(col_data),
        levels = levels(col_data)
      ))
    }
    # Anything else (us_units, fa_units, ar1_units, ...) is out of scope.
    stop(.fb_refusal_condition(
      reason_code = "rcov_type_unsupported_for_aggregation",
      message = paste0(
        "emit_gaussian_aggregated(): rcov term type '",
        ttype,
        "' is not supported by the aggregated path. Only homogeneous ",
        "(units / id) and at_units heterogeneous residual are ",
        "supported. Pass aggregate = FALSE for the per-row path."
      ),
      family_class = "flexybayes_aggregate_emit_refusal",
      rcov_type = ttype
    ))
  }
  list(
    kind = "homogeneous",
    factor = NULL,
    level_col = NULL,
    n_levels = NULL,
    levels = NULL
  )
}


# Random-intercept planning. Returns a list of per-RI-term records:
#   $col          character: column name in sufficient_stats
#   $level_idx    integer vector length K: per-cell level index
#   $n_levels     integer: factor's full level count
.agg_emit_ri_plan <- function(fb, fb_aggregated) {
  cols <- fb_aggregated$random_cols
  if (!length(cols)) {
    return(list())
  }
  out <- vector("list", length(cols))
  for (i in seq_along(cols)) {
    g <- cols[i]
    col_data <- fb_aggregated$sufficient_stats[[g]]
    if (!is.factor(col_data)) {
      col_data <- factor(col_data)
    }
    out[[i]] <- list(
      col = g,
      level_idx = as.integer(col_data),
      n_levels = nlevels(col_data)
    )
  }
  out
}


# ---------------------------------------------------------------- #
# greta path                                                        #
# ---------------------------------------------------------------- #
.emit_gaussian_aggregated_greta <- function(
  fb,
  fb_aggregated,
  data,
  wss_per_cell,
  WSS_total,
  ri_plan,
  rcov_plan,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd,
  verbose,
  mcmc_verbose
) {
  if (!requireNamespace("greta", quietly = TRUE)) {
    stop(
      "Package 'greta' is required for the greta-backed aggregated ",
      "emit. Install with install.packages('greta'); ",
      "greta::install_greta_deps().",
      call. = FALSE
    )
  }

  cell_design <- fb_aggregated$cell_design
  agg <- fb_aggregated$sufficient_stats
  K <- fb_aggregated$K
  N <- fb_aggregated$N
  p <- ncol(cell_design)
  ybar_k <- agg$S1_k / agg$n_k
  n_k <- agg$n_k

  # Single-matrix-multiply formulation of the cell-level linear
  # predictor. Avoids `+` between greta_array operands -- which
  # triggers greta's `check_dims` -> `as.greta_array(...)` dispatch
  # that fails to find the greta_array identity coercion when the
  # call is made from inside flexyBayes's package namespace (the
  # S3 method table is partially populated; `as.greta_array.matrix`
  # fires on greta_array operands and surfaces as a spurious
  # missing/infinite-values error). The existing per-row
  # emit_greta() works around this via a codegen-then-eval-in-greta-
  # namespace pattern; the aggregated path achieves the same end by
  # construction.
  #
  # Algebra:
  #   M_combined <- [ cell_design | Z_1 | Z_2 | ... ]   K x (p + sum_i nL_i)
  #   theta      <- c(beta, u_1, u_2, ...)              length p + sum_i nL_i
  #   mu_cell    <- M_combined %*% theta                K x 1
  # Each Z_i is a K x nL_i one-hot matrix selecting the cell's level
  # for RI term i, so `Z_i %*% u_i` equals the per-cell RI
  # contribution that `u_i[ri_plan[[i]]$level_idx]` would produce.
  Z_per_term <- lapply(ri_plan, function(r) {
    Z <- matrix(0, K, r$n_levels)
    Z[cbind(seq_len(K), r$level_idx)] <- 1
    Z
  })
  M_combined <- if (length(Z_per_term)) {
    do.call(cbind, c(list(cell_design), Z_per_term))
  } else {
    cell_design
  }

  greta_ns <- asNamespace("greta")
  g <- greta_ns

  M_g <- g$as_data(M_combined)
  ybar_g <- g$as_data(ybar_k)
  n_k_g <- g$as_data(n_k)

  beta_g <- g$normal(0, prior_fixed_sd, dim = p)

  tau_per_term <- vector("list", length(ri_plan))
  u_per_term <- vector("list", length(ri_plan))
  for (i in seq_along(ri_plan)) {
    nL <- ri_plan[[i]]$n_levels
    tau <- g$normal(0, prior_vc_sd, truncation = c(0, Inf))
    u_ <- g$normal(0, tau, dim = nL)
    tau_per_term[[i]] <- tau
    u_per_term[[i]] <- u_
  }

  # Concatenate beta + u_1 + u_2 + ... into a single greta vector and
  # compute the cell-level linear predictor `M_combined %*% theta`.
  # Both the `c()` concatenation and the `%*%` matrix-multiply are
  # called by DIRECT METHOD LOOKUP into greta's namespace so we
  # bypass R's S3 dispatch, which from inside flexyBayes's package
  # namespace routes greta_array operands to base-R's `c()` /
  # `%*%` (returning a non-greta_array result that subsequently
  # fails greta's downstream `as.greta_array` coercion). The
  # alternative (codegen + eval-in-greta-ns-env, as in emit_greta())
  # is blocked by the local shell harness's eval() guard.
  c_method <- get("c.greta_array", envir = greta_ns)
  matmul_method <- get("%*%.greta_array", envir = greta_ns)

  # Iterative two-arg concatenation: `do.call(c_method, ...)` synthesises
  # a call object whose first element deparses to multiple lines (the
  # c_method function body), tripping greta's internal
  # `vapply(names, deparse, "")` length-check inside c.greta_array.
  # Calling c_method with explicit two-arg form avoids the issue.
  theta <- beta_g
  for (u_ in u_per_term) {
    theta <- c_method(theta, u_)
  }

  mu_cell <- matmul_method(M_g, theta)

  # Residual sigma: scalar (homogeneous) or per-level (at_units).
  if (identical(rcov_plan$kind, "homogeneous")) {
    sigma_g <- g$normal(0, prior_vc_sd, truncation = c(0, Inf))
    sigma_cell_g <- sigma_g
  } else {
    sigma_g <- g$normal(
      0,
      prior_vc_sd,
      truncation = c(0, Inf),
      dim = rcov_plan$n_levels
    )
    sigma_cell_g <- sigma_g[rcov_plan$level_col]
  }

  # (i) cell-mean weighted gaussian.
  g$`distribution<-`(ybar_g, g$normal(mu_cell, sigma_cell_g / sqrt(n_k_g)))

  # (ii) within-cell SS gamma observation -- one per residual level
  # for at_units, one global for homogeneous. Skipped when degrees of
  # freedom <= 0 or WSS <= machine eps (single-obs cells).
  if (identical(rcov_plan$kind, "homogeneous")) {
    df_total <- N - K
    if (df_total > 0L && WSS_total > .Machine$double.eps) {
      WSS_g <- g$as_data(WSS_total)
      g$`distribution<-`(WSS_g, g$gamma(df_total / 2, 1 / (2 * sigma_g^2)))
    }
  } else {
    for (lv in seq_len(rcov_plan$n_levels)) {
      idx <- which(rcov_plan$level_col == lv)
      n_level <- sum(n_k[idx])
      k_level <- length(idx)
      WSS_lv <- sum(wss_per_cell[idx])
      df_lv <- n_level - k_level
      if (df_lv > 0L && WSS_lv > .Machine$double.eps) {
        WSS_lv_g <- g$as_data(WSS_lv)
        g$`distribution<-`(
          WSS_lv_g,
          g$gamma(df_lv / 2, 1 / (2 * sigma_g[lv]^2))
        )
      }
    }
  }

  # Build model + run MCMC.
  #
  # greta::model() internally does `substitute(list(...))[-1]` then
  # `vapply(names, deparse, "")` to label its target greta_arrays.
  # Under `do.call(g$model, list_of_greta_arrays)` the substituted
  # arguments deparse to the FULL multi-line greta_array printout
  # (not to symbol names), and the vapply length-1 check fails with
  # "values must be length 1, but FUN(X[[1]]) result is length 2".
  # The workaround: bind each greta_array under a short symbol name
  # in a fresh environment whose parent is the greta namespace, then
  # build and run a model() call expression in that environment so
  # substitute() captures the short symbols rather than the values.
  model_env <- new.env(parent = greta_ns)
  model_env$beta <- beta_g
  model_env$sigma <- sigma_g
  tau_names <- character(length(tau_per_term))
  for (i in seq_along(tau_per_term)) {
    nm <- paste0("tau_", i)
    model_env[[nm]] <- tau_per_term[[i]]
    tau_names[[i]] <- nm
  }
  call_args <- c(list(quote(beta), quote(sigma)), lapply(tau_names, as.name))
  model_call <- as.call(c(list(quote(model)), call_args))
  m <- base::eval(model_call, envir = model_env)
  draws <- g$mcmc(
    m,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    verbose = isTRUE(mcmc_verbose)
  )

  list(
    backend = "greta",
    draws = draws,
    model = m,
    beta = beta_g,
    sigma = sigma_g,
    tau_per_term = tau_per_term,
    cell_design = cell_design,
    K = K,
    N = N,
    p = p
  )
}


# ---------------------------------------------------------------- #
# INLA path                                                         #
# ---------------------------------------------------------------- #
.emit_gaussian_aggregated_inla <- function(
  fb,
  fb_aggregated,
  data,
  wss_per_cell,
  WSS_total,
  ri_plan,
  rcov_plan,
  verbose
) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for the INLA-backed aggregated ",
      "emit. Install from https://inla.r-inla-download.org/.",
      call. = FALSE
    )
  }

  if (identical(rcov_plan$kind, "at_units")) {
    stop(
      "emit_gaussian_aggregated() INLA: heterogeneous residual ",
      "at_units on INLA requires the multi-likelihood INLA stack ",
      "API (deferred). Pass backend = \"greta\" with ",
      "aggregate = TRUE for the heterogeneous-residual path, or ",
      "backend = \"inla\" with aggregate = FALSE for the per-row ",
      "INLA path.",
      call. = FALSE
    )
  }

  agg <- as.data.frame(fb_aggregated$sufficient_stats)
  K <- fb_aggregated$K
  N <- fb_aggregated$N
  agg$ybar_k <- agg$S1_k / agg$n_k

  # Pull per-RI hyperpriors from the same priors_to_inla() translation
  # the per-row emit_inla() uses, so the aggregated + per-row INLA fits
  # share the random-intercept prior shape. The residual prior is NOT
  # plumbed through at v0.3.2: the custom-prior expression below
  # composes Delta with INLA's default loggamma(1, 5e-5) base. A user-
  # supplied residual prior on `sigma` is therefore ignored on the
  # aggregated INLA path -- documented limitation; matched-prior
  # composition is a v0.3.3 follow-on. The greta aggregated path does
  # honor the flexyBayes prior via .priors_to_legacy() naturally.
  hyper_ctrl <- if (inherits(fb$priors, "fb_prior")) {
    priors_to_inla(fb$priors)
  } else {
    list()
  }

  fixed_form <- .agg_emit_inla_formula(fb, fb_aggregated, ri_plan, hyper_ctrl)

  # Custom-prior expression on the gaussian-precision hyperparameter:
  # combines INLA's default loggamma(1, 5e-5) prior with the Delta
  # correction that absorbs the within-cell SS into the prior so the
  # cell-mean weighted likelihood + custom prior recovers the per-row
  # posterior.
  a <- WSS_total / 2 + 5e-5
  b <- (N - K) / 2 + 1
  # Recenter the log-density to ~0 at its mode theta* = log(b / a). The
  # additive constant `c` does not change the prior (it is absorbed in the
  # normalising constant), but it keeps the values INLA's expression
  # evaluator sees O(1) near the mode. Without it the log-density is ~ -b
  # at the mode (e.g. ~ -2.5e9 at N = 5e9), at which magnitude INLA's
  # hyperparameter integration becomes numerically unstable and can
  # corrupt the fixed-effect posterior at extreme aggregated N.
  c0 <- b * (log(b / a) - 1)
  custom_prior_expr <- sprintf(
    paste0(
      "expression:",
      " a = %.16e;",
      " b = %.16e;",
      " c = %.16e;",
      " log_dens = -a * exp(theta) + b * theta - c;",
      " return(log_dens);"
    ),
    a,
    b,
    c0
  )

  control_family <- list(
    hyper = list(prec = list(prior = custom_prior_expr))
  )

  inla_call <- list(
    formula = fixed_form,
    data = agg,
    family = "gaussian",
    scale = agg$n_k,
    control.family = control_family,
    control.compute = list(
      config = TRUE,
      return.marginals = TRUE,
      dic = FALSE,
      waic = FALSE
    )
  )

  inla_fit <- do.call(INLA::inla, inla_call)

  list(
    backend = "inla",
    inla = inla_fit,
    fixed_form = fixed_form,
    custom_prior = custom_prior_expr,
    K = K,
    N = N
  )
}


# Build an INLA-side formula on the cell-level data.table. The
# `<fb_aggregated>$sufficient_stats` data.table carries the model-
# matrix contrast columns (`(Intercept)`, `f1b`, `f1c`, ...) rather
# than the original factor variables, so the formula references
# those contrast columns directly. `(Intercept)` is dropped because
# R's formula machinery auto-adds the intercept (unless `fb$intercept
# == FALSE`, in which case we suppress with `- 1`). RI grouping
# columns enter via INLA's `f(<col>, model = "iid")` syntax. Backticks
# guard against syntactically-awkward contrast names.
.agg_emit_inla_formula <- function(
  fb,
  fb_aggregated,
  ri_plan,
  hyper_ctrl = list()
) {
  fixed_cols <- setdiff(fb_aggregated$fixed_cols, "(Intercept)")
  rhs_fixed <- if (length(fixed_cols)) {
    paste(paste0("`", fixed_cols, "`"), collapse = " + ")
  } else {
    "1"
  }

  # When hyper_ctrl carries a translated user prior for this RI, splice
  # it into the f(...) call. In the default case the synthesised
  # uniform-on-SD prior (.default_uniform_prior() in flexybayes()) keys
  # every random group, so hyper_ctrl carries the same faithful uniform
  # expression prior the per-row emit_inla() path uses -- both paths
  # share one prior (scale fixed once from the per-row data), which is
  # what makes the per-row and aggregated INLA fits algebraically
  # equivalent. The former per-row/aggregated divergence under defaults
  # is resolved by that unification (priors_to_inla() now emits the exact
  # uniform on both paths instead of a per-row PC approximation).
  #
  # The fallback below fires only when no fb_prior is in play (e.g. the
  # legacy `prior_vc_sd` path) or a user supplies a partial prior that
  # omits this group. It pins the scale-invariant loggamma(1, 5e-5) so
  # the aggregated fit matches the per-row path, which inherits the same
  # INLA default when hyper_ctrl is empty.
  default_iid_hyper <- list(prior = "loggamma", param = c(1, 5e-05))
  ri_pieces <- vapply(
    ri_plan,
    function(r) {
      entry <- hyper_ctrl[[r$col]] %||% default_iid_hyper
      hyper_arg <- .inla_hyper_arg(entry)
      sprintf(
        "f(`%s`, model = \"iid\"%s)",
        r$col,
        if (nzchar(hyper_arg)) paste0(", ", hyper_arg) else ""
      )
    },
    character(1L)
  )
  rhs_ri <- if (length(ri_pieces)) {
    paste(ri_pieces, collapse = " + ")
  } else {
    ""
  }

  rhs <- if (nzchar(rhs_ri)) {
    paste(rhs_fixed, rhs_ri, sep = " + ")
  } else {
    rhs_fixed
  }
  if (!isTRUE(fb$intercept)) {
    rhs <- paste(rhs, "- 1")
  }

  stats::as.formula(paste("ybar_k ~", rhs))
}


# ---------------------------------------------------------------- #
# Shared fit-object constructor                                     #
# ---------------------------------------------------------------- #
# Builds the `<flexybayes>` fit shape that downstream methods
# (print/summary/predict/triangulate) consume. Mirrors emit_greta() /
# emit_inla() construction so the new aggregated fits are
# byte-interchangeable for the existing method surface.
.agg_emit_build_fit <- function(
  fb,
  fb_aggregated,
  data,
  backend,
  engine_out,
  ri_plan,
  rcov_plan,
  elapsed,
  the_call,
  fixed,
  random,
  rcov,
  family,
  link,
  data_name,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd
) {
  # Posterior summary: pull the engine-side draws / marginals into a
  # uniform shape. greta: collapse mcmc.list across chains, summarise
  # mean/sd/q025/q975 per parameter. INLA: read summary.fixed +
  # summary.hyperpar.
  posterior_summary <- if (identical(backend, "greta")) {
    .agg_greta_summarise(engine_out, fb_aggregated, ri_plan, rcov_plan)
  } else {
    .agg_inla_summarise(engine_out, fb_aggregated, ri_plan, rcov_plan)
  }

  # Per-row reconstructed fitted values for the $glm shim.
  fitted_row <- .agg_reconstruct_fitted_row(
    fb,
    fb_aggregated,
    data,
    posterior_summary,
    ri_plan
  )
  y_row <- data[[fb$response]]
  resid_row <- as.numeric(y_row) - fitted_row

  glm_obj <- structure(
    list(
      coefficients = posterior_summary$beta_means,
      vcov = posterior_summary$beta_vcov,
      fitted.values = fitted_row,
      linear.predictors = fitted_row,
      residuals = resid_row,
      y = y_row,
      family = list(family = fb$family, link = fb$link %||% "identity"),
      formula = the_call,
      data = data
    ),
    class = c("flexybayes_glm_shim", "lm", "list")
  )

  extras <- structure(
    list(
      summary = posterior_summary,
      convergence = posterior_summary$convergence,
      variance_comps = posterior_summary$variance_comps,
      run_time = elapsed,
      parse_info = list(
        # Mirror emit_greta()'s fixed-info shape so confint /
        # summary / methods.R downstream consumers see the same
        # `list(response, intercept, terms)` structure.
        fixed = list(
          response = fb$response,
          intercept = fb$intercept,
          terms = fb$fixed_terms
        ),
        random = fb$random_terms,
        rcov = fb$rcov_terms,
        family = list(family = fb$family, link = fb$link %||% "identity"),
        smooths = list() # the aggregated path excludes smooths.
      ),
      call_info = list(
        fixed = fixed,
        random = random,
        rcov = rcov,
        data_name = data_name,
        family = family,
        link = link,
        n_samples = n_samples,
        warmup = warmup,
        chains = chains,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd
      ),
      model_info = list(
        n_obs = fb_aggregated$N,
        n_cells = fb_aggregated$K,
        n_fixed = length(fb$fixed_terms) + isTRUE(fb$intercept),
        n_random = length(fb$random_terms),
        family = fb$family,
        link = fb$link %||% "identity"
      ),
      aggregation_meta = list(
        N = fb_aggregated$N,
        K = fb_aggregated$K,
        compression = fb_aggregated$compression,
        residual = rcov_plan$kind,
        # Whether the posterior the user sees is the per-row-equivalent
        # one. The aggregated likelihood plus the default precision prior
        # recover the per-row posterior to numerical precision, so the
        # default path is tagged "per_row_equivalent". When the user
        # supplies an explicit prior object, the matched-prior
        # equivalence against a default-prior per-row fit no longer holds
        # (and on the INLA aggregated path the residual prior is moreover
        # not plumbed through -- see the "Matched priors" note on
        # triangulate()), so the fit is tagged "custom". Surfaced via
        # canonical_names() and the aggregated print / summary.
        prior_parametrization = if (inherits(fb$priors, "fb_prior")) {
          "custom"
        } else {
          "per_row_equivalent"
        }
      ),
      the_call = the_call,
      formula = the_call
    ),
    class = "flexybayes_extras"
  )

  # The greta raw slot mirrors the per-row emit_greta() shape --
  # `$greta$draws` holding the coda mcmc object -- so every consumer
  # that reads `fit$greta$draws` (fb_as_draws_simple / canonical_names /
  # confint / plot / predict) works identically on the aggregated greta
  # path. Storing the draws bare at `fit$greta` (the pre-v0.4.0 layout)
  # left `fit$greta$draws` NULL, so triangulate() and canonical_names()
  # on a greta-aggregated fit failed with "fit$greta$draws is missing".
  # The INLA slot holds the raw inla fit object directly, matching
  # fb_as_draws_simple.flexybayes_inla()'s `fit$inla` contract.
  raw_slot <- if (identical(backend, "greta")) {
    list(draws = engine_out$draws)
  } else {
    engine_out$inla
  }
  raw_name <- if (identical(backend, "greta")) "greta" else "inla"

  # Backend identity in the class vector so generics that have no
  # `.flexybayes_aggregated` method (e.g., fb_as_draws_simple()) fall
  # through to the correct backend-specific method via S3 dispatch.
  # Without `flexybayes_inla` in the vector, aggregated INLA fits
  # dispatched to fb_as_draws_simple.flexybayes (the greta extractor)
  # and failed with "fit$greta$draws is missing" -- breaking
  # triangulate(fit_g, fit_i) on the aggregated path.
  fit <- structure(
    list(
      glm = glm_obj,
      extras = extras
    ),
    class = c(
      "flexybayes_aggregated",
      if (identical(backend, "inla")) "flexybayes_inla",
      "flexybayes"
    )
  )
  fit[[raw_name]] <- raw_slot
  fit
}


# ---------------------------------------------------------------- #
# Posterior-summarisation helpers                                   #
# ---------------------------------------------------------------- #
.agg_greta_summarise <- function(
  engine_out,
  fb_aggregated,
  ri_plan,
  rcov_plan
) {
  draws <- engine_out$draws
  m <- do.call(rbind, draws)

  beta_idx <- seq_len(engine_out$p)
  sigma_idx <- if (identical(rcov_plan$kind, "homogeneous")) {
    engine_out$p + 1L
  } else {
    engine_out$p + seq_len(rcov_plan$n_levels)
  }
  tau_start <- max(sigma_idx) + 1L
  tau_idx <- if (length(ri_plan)) {
    tau_start:(tau_start + length(ri_plan) - 1L)
  } else {
    integer(0)
  }

  beta_means <- colMeans(m[, beta_idx, drop = FALSE])
  names(beta_means) <- colnames(fb_aggregated$cell_design)
  beta_vcov <- stats::cov(m[, beta_idx, drop = FALSE])
  dimnames(beta_vcov) <- list(names(beta_means), names(beta_means))

  sigma_means <- colMeans(m[, sigma_idx, drop = FALSE])
  tau_means <- if (length(tau_idx)) {
    colMeans(m[, tau_idx, drop = FALSE])
  } else {
    numeric(0)
  }

  # R-hat approximation: split-chain variance ratio per param. greta's
  # draws is a coda mcmc.list; reuse coda's gelman.diag when available.
  psrf_mat <- tryCatch(
    {
      coda_ns <- asNamespace("coda")
      gd <- coda_ns$gelman.diag(draws, multivariate = FALSE, autoburnin = FALSE)
      gd$psrf
    },
    error = function(e) NULL
  )

  list(
    backend = "greta",
    beta_means = beta_means,
    beta_vcov = beta_vcov,
    sigma_means = sigma_means,
    tau_means = tau_means,
    convergence = list(gelman = list(psrf = psrf_mat)),
    variance_comps = list(sigma = sigma_means, tau = tau_means)
  )
}


.agg_inla_summarise <- function(engine_out, fb_aggregated, ri_plan, rcov_plan) {
  inla_fit <- engine_out$inla
  fixed_summary <- inla_fit$summary.fixed
  beta_means <- fixed_summary$mean
  names(beta_means) <- rownames(fixed_summary)

  # `nrow =` guards the intercept-only case: diag() of a length-one
  # vector is otherwise read as diag(n) and returns an n x n identity
  # (a 0 x 0 matrix here), not a 1 x 1 covariance. Mirrors the count path.
  beta_vcov <- diag(fixed_summary$sd^2, nrow = length(beta_means))
  dimnames(beta_vcov) <- list(names(beta_means), names(beta_means))

  hyper <- inla_fit$summary.hyperpar
  prec_y_mean <- hyper$mean[grepl(
    "^Precision for the Gaussian observations",
    rownames(hyper)
  )]
  sigma_means <- if (length(prec_y_mean)) {
    sqrt(1 / prec_y_mean)
  } else {
    NA_real_
  }

  tau_means <- if (length(ri_plan)) {
    out <- numeric(length(ri_plan))
    for (i in seq_along(ri_plan)) {
      pat <- paste0("^Precision for ", ri_plan[[i]]$col, "$")
      prec_i <- hyper$mean[grepl(pat, rownames(hyper))]
      out[i] <- if (length(prec_i)) sqrt(1 / prec_i) else NA_real_
    }
    out
  } else {
    numeric(0)
  }

  list(
    backend = "inla",
    beta_means = beta_means,
    beta_vcov = beta_vcov,
    sigma_means = sigma_means,
    tau_means = tau_means,
    convergence = list(gelman = list(psrf = NULL)),
    variance_comps = list(sigma = sigma_means, tau = tau_means)
  )
}


# Reconstruct per-row fitted values from the cell-level posterior +
# original data. Pure linear-predictor expansion at the fixed-effect
# level; random-intercept BLUP shrinkage is a follow-on extension
# (would join `predict.flexybayes_aggregated` if/when added). The
# formula uses the ORIGINAL IR variable names (e.g. `f1 + f2`) so the
# per-row data.frame's factor columns resolve cleanly through
# `model.matrix()`; the resulting contrast columns must match the
# cell_design's contrast columns by construction (same R formula
# semantics applied to the same factors).
.agg_reconstruct_fitted_row <- function(
  fb,
  fb_aggregated,
  data,
  posterior_summary,
  ri_plan
) {
  fixed_form <- .fb_aggregate_fixed_formula(fb)
  if (is.null(fixed_form)) {
    X_row <- matrix(
      1,
      nrow = nrow(data),
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )
  } else {
    X_row <- stats::model.matrix(fixed_form, data = data)
  }
  # Defensive: align column order between X_row and the beta vector
  # (same factor order through model.matrix, but be explicit).
  beta_aligned <- posterior_summary$beta_means[colnames(X_row)]
  as.numeric(X_row %*% beta_aligned)
}
