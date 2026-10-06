# fb_refit.R: the model of a fit, refitted on new data with what a package
# built on flexyBayes adds to it.

# Arguments of flexybayes() that state the model; fb_refit() takes them from
# the fit.
.FB_REFIT_MODEL_ARGS <- c(
  "random", "residual", "family", "link", "known_matrices", "prior",
  "prior_fixed_sd", "prior_vc_sd", "syntax"
)


#' Refit a model on new data, with added fixed columns
#'
#' Refits the model of `fit` on `data`, with the same random structure,
#' family, priors, engine and engine settings, known variance or covariance,
#' observation weights and offsets. A package built on flexyBayes uses it to
#' refit a user's model on a fold of the data, with candidate columns added
#' to the fixed part, or with brms prior rows and Stan code of its own.
#'
#' @section What the refit keeps:
#' * The fixed formula is the one `fit` was called with, or `fixed`, with
#'   the terms of `add_fixed` appended.
#' * A prior written with [fb_prior()], `prior_fixed_sd` or `prior_vc_sd` is
#'   passed on as written. The package default prior is built again over the
#'   refitted model, as [update()] does, so on the same data every term keeps
#'   the prior it had and an added term gets the default.
#' * The engine is the one `fit` was fitted with, `"inla"` or `"brms"`, even
#'   when the original call asked for `"auto"`. Its recorded sampler
#'   settings, `seed`, `control`, `na_action`, `fixed_alias` and `aggregate`
#'   are passed on.
#' * Weights and a known variance belong to rows. A refit on the fit's own
#'   data keeps them. A known variance is read from the column of that name
#'   in the fit's data, so rows of `fit$data` carry it. For any other `data`,
#'   pass `weights = ` or `known_variance = ` with one value per row;
#'   without them a weighted fit is refused, as is a fit with a known
#'   variance whose `data` lacks the column.
#'
#' @section Extending the brms model:
#' `extra_prior`, `stanvars` and `formula_hook` are brms objects, and the
#' refit then runs on brms: `fit` must be a brms fit, or `backend = "brms"`
#' must be given. The refit still starts from the brms formula and prior
#' rows the package builds for the model.
#'
#' `formula_hook` is called as `formula_hook(formula, prior)` with that
#' formula and the `brmsprior` of those rows, and returns
#' `list(formula = , prior = )`, the formula and prior rows to fit. A
#' non-linear formula usually moves the package's rows onto one of its
#' parts, which is why the hook receives them. `extra_prior` rows are then
#' merged in, replacing any row with the same `class`, `coef`, `group`,
#' `resp`, `dpar` and `nlpar`, and `stanvars` is passed to [brms::brm()].
#' The fitted object reads under the linear model's names.
#'
#' @section Readout:
#' `readout` names what the caller reads from an INLA refit. `"varcomp"`
#' (variance components), `"fitted"` (fitted values) and `"draws"`
#' (posterior draws, for example on held-out rows) ask INLA for that
#' readout only: no DIC, WAIC or marginal densities, posterior
#' configurations only for `"draws"`, and no second fit to test an AR1
#' field. The posterior mode, and the readout taken from it, are those of
#' the full refit. A brms refit is always returned in full.
#'
#' @param fit A flexyBayes fit.
#' @param data A data frame with the columns the model uses and any added
#'   ones. `NULL`, the default, refits on the data the fit carries.
#' @param add_fixed Fixed-effect terms to append to the fixed formula: a
#'   one-sided formula such as `~ z1 + z2`, or a character vector of term
#'   labels. `NULL` adds none.
#' @param extra_prior A `brmsprior` from [brms::set_prior()] with rows to add
#'   to the prior, or `NULL`.
#' @param stanvars A `stanvars` object from [brms::stanvar()], or `NULL`.
#' @param formula_hook A function of `formula` and `prior` returning
#'   `list(formula = , prior = )`, or `NULL`. See "Extending the brms model".
#' @param readout `"full"`, the default, or `"varcomp"`, `"fitted"` or
#'   `"draws"`. See "Readout".
#' @param fixed A two-sided formula that replaces the fit's fixed formula,
#'   for a refit whose fixed part differs by more than added terms. `NULL`,
#'   the default, keeps the fit's.
#' @param ... Named arguments of [flexybayes()] that change how this refit
#'   runs, for example `seed`, `backend`, `n_samples`, `control`,
#'   `na_action`, `verbose`, `return_code`, `weights` or `known_variance`.
#'   On the brms route of "Extending the brms model", arguments that
#'   [flexybayes()] does not have go to [brms::brm()]. Arguments that state
#'   the model (`random`, `residual`, `family`, `link`, `known_matrices`,
#'   `prior`, `prior_fixed_sd`, `prior_vc_sd`) are refused; [update()]
#'   refits a changed model.
#' @returns A fitted flexyBayes object of the engine the refit ran on, or,
#'   with `return_code = TRUE`, the code that engine would run.
#' @seealso [update.flexybayes()] to refit with any argument changed.
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(site = factor(rep(1:12, each = 10)), x = rnorm(120),
#'                   z = rnorm(120))
#'   d$y <- rnorm(12)[d$site] + 0.5 * d$x + rnorm(120)
#'   fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
#'                     verbose = FALSE)
#'   # The same model with z added, on the first nine sites.
#'   fold <- fb_refit(fit, data = d[d$site %in% 1:9, ], add_fixed = ~z,
#'                    verbose = FALSE)
#'   coef(fold)
#' }
#' }
#' @export
fb_refit <- function(
  fit,
  data = NULL,
  add_fixed = NULL,
  extra_prior = NULL,
  stanvars = NULL,
  formula_hook = NULL,
  readout = c("full", "varcomp", "fitted", "draws"),
  fixed = NULL,
  ...
) {
  if (!inherits(fit, "flexybayes")) {
    stop("fb_refit(): `fit` must be a flexyBayes fit.", call. = FALSE)
  }
  readout <- match.arg(readout)
  dots <- list(...)
  if (length(dots) > 0L &&
        (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop(
      "fb_refit(): every argument in `...` must be named, as an argument ",
      "of flexybayes().",
      call. = FALSE
    )
  }
  stated <- intersect(names(dots), .FB_REFIT_MODEL_ARGS)
  if (length(stated) > 0L) {
    stop(
      "fb_refit(): ", paste0("`", stated, "`", collapse = ", "),
      " would change the model the fit states. Use update() to refit a ",
      "changed model.",
      call. = FALSE
    )
  }
  if (!is.null(formula_hook) && !is.function(formula_hook)) {
    stop("fb_refit(): `formula_hook` must be a function or NULL.",
         call. = FALSE)
  }

  info <- fit$extras$call_info %||% list()
  ir <- .fb_fit_ir(fit)
  fixed <- .fb_refit_fixed(
    fixed %||% info$fixed,
    .fb_refit_add_labels(add_fixed)
  )
  rows <- .fb_refit_rows(fit, info, ir, data, dots)

  hooked <- !is.null(extra_prior) || !is.null(stanvars) ||
    !is.null(formula_hook)
  if (hooked) {
    return(.fb_refit_brms(
      fit, info, ir, fixed, rows, extra_prior, stanvars, formula_hook, dots
    ))
  }

  priors <- ir$priors
  if (!is.null(attr(priors, "fb_prior_default_basis"))) {
    priors <- NULL
  }
  args <- list(
    fixed = fixed,
    random = info$random,
    residual = info$residual,
    data = rows$data,
    family = info$family %||% ir$family,
    link = info$link,
    known_matrices = info$known_matrices %||% list(),
    prior = priors,
    seed = info$seed,
    control = info$control,
    backend = .fb_refit_engine(fit, info)
  )
  # A scalar prior argument is passed only when the fit's call gave it, so a
  # fit on the package default refits on the default.
  for (which in c("fixed_sd", "vc_sd")) {
    if (.fb_prior_scalar_supplied(ir, which)) {
      args[[paste0("prior_", which)]] <- .fb_prior_scalar_value(ir, which)
    }
  }
  for (nm in c("na_action", "fixed_alias", "aggregate", "n_samples",
               "warmup", "chains")) {
    if (!is.null(info[[nm]])) {
      args[[nm]] <- info[[nm]]
    }
  }
  args$weights <- rows$weights
  args$known_variance <- rows$known_variance
  for (nm in names(dots)) {
    args[nm] <- list(dots[[nm]])
  }

  plan <- if (identical(readout, "full")) NULL else .fb_fold_refit_plan(readout)
  .fb_with_fold_refit(plan, do.call(flexybayes, args))
}


# The engine a fit ran on. The call record holds the routing request, which
# may be "auto".
.fb_refit_engine <- function(fit, info) {
  if (inherits(fit, "flexybayes_brms")) {
    return("brms")
  }
  if (inherits(fit, "flexybayes_inla")) {
    return("inla")
  }
  info$backend %||% "auto"
}


.fb_refit_add_labels <- function(add_fixed) {
  if (is.null(add_fixed)) {
    return(character(0))
  }
  if (inherits(add_fixed, "formula") && length(add_fixed) == 2L) {
    return(attr(stats::terms(add_fixed), "term.labels"))
  }
  if (is.character(add_fixed) && !anyNA(add_fixed) && all(nzchar(add_fixed))) {
    return(add_fixed)
  }
  stop(
    "fb_refit(): `add_fixed` must be a one-sided formula, such as ",
    "`~ z1 + z2`, or a character vector of term labels.",
    call. = FALSE
  )
}


# `fixed` with `add` appended to its right-hand side; a brmsformula keeps its
# other parts.
.fb_refit_fixed <- function(fixed, add) {
  if (is.null(fixed)) {
    stop(.fb_refusal_condition(
      reason_code = "refit_not_reconstructable",
      message = paste0(
        "fb_refit(): the fit does not record the fixed formula it was ",
        "called with, so its model cannot be rebuilt. Pass `fixed`."
      )
    ))
  }
  if (length(add) == 0L) {
    return(fixed)
  }
  if (inherits(fixed, "brmsformula")) {
    fixed$formula <- .fb_refit_fixed(fixed$formula, add)
    return(fixed)
  }
  if (!inherits(fixed, "formula") || length(fixed) != 3L) {
    stop(
      "fb_refit(): the fixed formula is not a two-sided formula, so ",
      "`add_fixed` cannot be appended to it.",
      call. = FALSE
    )
  }
  lhs <- paste(deparse(fixed[[2L]]), collapse = " ")
  rhs <- paste(deparse(fixed[[3L]]), collapse = " ")
  stats::as.formula(
    paste0(lhs, " ~ ", rhs, " + ", paste(add, collapse = " + ")),
    env = environment(fixed)
  )
}


# The rows of the refit: its data, and the weights and known variance that
# go with them.
.fb_refit_rows <- function(fit, info, ir, data, dots) {
  own <- is.null(data)
  data <- data %||% fit$glm$data %||% fit$data
  if (is.null(data)) {
    stop(.fb_refusal_condition(
      reason_code = "refit_not_reconstructable",
      message = paste0(
        "fb_refit(): this fit does not keep the data it was fitted to. ",
        "Pass `data`."
      )
    ))
  }
  data <- as.data.frame(data)

  weights <- NULL
  if ("weights" %in% names(dots)) {
    weights <- dots$weights
  } else {
    w <- .fb_ir_weights(ir) %||% info$weights
    if (.fb_weights_nonconstant(w)) {
      if (!own) {
        stop(.fb_refusal_condition(
          reason_code = "refit_not_reconstructable",
          message = paste0(
            "fb_refit(): the fit carries observation weights, one per row ",
            "of its own data, and new `data` was given without them. Pass ",
            "`weights = ` with one weight per row of `data`."
          )
        ))
      }
      weights <- w
    }
  }

  known_variance <- NULL
  if ("known_variance" %in% names(dots)) {
    known_variance <- dots$known_variance
  } else if (!is.null(.fb_known_variance_term(ir))) {
    known_variance <- .fb_known_variance_reissue(ir, data)
    if (is.null(known_variance)) {
      .fb_known_variance_refuse(paste0(
        "fb_refit(): the fit carries a known variance per row, and `data` ",
        "has no `", .FB_KNOWN_VARIANCE_COL, "` column to read it from. ",
        "Pass rows of the fit's own data (`fit$data`), or ",
        "`known_variance = ` with one value per row of `data`."
      ))
    }
  }
  data <- data[setdiff(names(data), .fb_known_variance_columns())]
  list(data = data, weights = weights, known_variance = known_variance)
}


# The refit on brms with extra prior rows, Stan code or a formula hook. The
# model is built as flexybayes() builds it and handed to emit_brms(), which
# takes these three; flexybayes() takes none of them.
.fb_refit_brms <- function(fit, info, ir, fixed, rows, extra_prior, stanvars,
                           formula_hook, dots) {
  engine <- if ("backend" %in% names(dots)) {
    dots$backend
  } else {
    .fb_refit_engine(fit, info)
  }
  if (!identical(.resolve_engine_string(engine), "brms")) {
    stop(.fb_refusal_condition(
      reason_code = "refit_hooks_require_brms",
      message = paste0(
        "fb_refit(): `extra_prior`, `stanvars` and `formula_hook` are brms ",
        "objects, and this refit would run on \"",
        .resolve_engine_string(engine), "\". Pass `backend = \"brms\"`."
      )
    ))
  }
  opts <- .fb_engine_opts(engine) %||% list()
  dots$backend <- NULL
  setting <- function(name) {
    if (name %in% names(dots)) dots[[name]] else info[[name]]
  }
  prior_fixed_sd <- info$prior_fixed_sd %||% 100
  prior_vc_sd <- info$prior_vc_sd %||% 1

  data <- rows$data
  fb <- .build_ir_polymorphic(
    fixed = fixed,
    random = info$random,
    residual = info$residual,
    data = data,
    family = info$family %||% ir$family,
    link = info$link %||% ir$link,
    weights = rows$weights,
    known_matrices = info$known_matrices %||% list(),
    prior = NULL,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    syntax = "auto"
  )
  # The fit's resolved prior is carried as it is: rebuilt here, the default
  # would be this release's and not the one the fit ran under.
  if (!is.null(ir$priors)) {
    fb$priors <- ir$priors
  }
  if (!is.null(ir$prior_scalars)) {
    fb$prior_scalars <- ir$prior_scalars
  }
  known <- .fb_known_variance_term(ir)
  if (!is.null(rows$known_variance)) {
    known <- known %||% list(type = "known_variance",
                             column = .FB_KNOWN_VARIANCE_COL)
    known$values <- as.numeric(rows$known_variance)
    data[[.FB_KNOWN_VARIANCE_COL]] <- known$values
    fb$addition_terms <- c(fb$addition_terms %||% list(), list(known))
  }

  passed <- c("seed", "control", "n_samples", "warmup", "chains",
              "na_action", "weights", "known_variance")
  args <- c(
    list(
      fb = fb,
      data = data,
      known_matrices = info$known_matrices %||% list(),
      weights = rows$weights,
      n_samples = opts$n_samples %||% setting("n_samples") %||% 1000L,
      warmup = opts$warmup %||% setting("warmup") %||% 500L,
      chains = opts$chains %||% setting("chains") %||% 4L,
      seed = setting("seed"),
      control = setting("control"),
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      extra_prior = extra_prior,
      stanvars = stanvars,
      formula_hook = formula_hook,
      fixed = fixed,
      random = info$random,
      residual = info$residual,
      family = info$family %||% ir$family,
      link = info$link %||% ir$link,
      data_name = info$data_name %||% NA_character_,
      na_action = setting("na_action"),
      requested_backend = "brms"
    ),
    dots[setdiff(names(dots), passed)]
  )
  out <- do.call(emit_brms, args)
  if (!is.list(out)) {
    return(out)
  }
  # The parsed model and the data, on the slots an INLA fit carries them in,
  # so readers of either engine's fit find them.
  out$fb <- fb
  out$data <- data
  out$extras$fb_terms <- fb
  out
}


# Calls a formula hook and checks what it returns.
.fb_refit_apply_formula_hook <- function(formula_hook, formula, prior) {
  out <- formula_hook(formula, prior)
  if (!is.list(out) || inherits(out, "brmsformula") ||
        !all(c("formula", "prior") %in% names(out))) {
    stop(
      "fb_refit(): `formula_hook` must return a list with elements ",
      "`formula` and `prior`.",
      call. = FALSE
    )
  }
  out[c("formula", "prior")]
}
