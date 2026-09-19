# fb_prior -- PC-canonical hybrid prior DSL.
#
# Implements the cross-engine prior interlingua. v0.1 minimum subset:
# user-facing constructor + S3 class + structured spec list +
# print method. Cross-engine translation tables for brms / INLA emit
# are in priors_to_inla() (this file) and priors_to_brms()
# (R/priors_to_brms.R).
#
# Targets supported in v0.1:
#   sigma             - residual standard deviation
#   sd(group = name)  - random-effect standard deviation
#   b("name")         - fixed-effect coefficient
#   cor(group = name) - correlation parameter
#   smooth(var)       - smoother variance (rw1/rw2)
#
# Distribution families supported in v0.1:
#   pc(upper, prob)            - PC prior (penalised complexity);
#                                exponential on sigma with
#                                rate = -log(prob) / upper.
#   half_normal(scale)         - on sd scale.
#   half_cauchy(scale)         - on sd scale.
#   student_t(df, scale)       - on sd scale (df > 0).
#   normal(mean, sd)           - on coefficient scale.
#   exponential(rate)          - on sd scale.
#   lkj(eta)                   - LKJ correlation prior.
#   uniform(lower, upper)      - on sd scale (lower >= 0, upper > lower).
#                                Sits outside the PC interlingua so cross-
#                                engine compilation is engine-direct, not
#                                via the PC translation table.
#
# All specifications are stored verbatim with their parsed args;
# the emit_*() functions translate to backend-specific syntax at
# fit time per the cross-engine translation table.

# ---------------------------------------------------------------- #
# Constructor                                                      #
# ---------------------------------------------------------------- #

#' Specify priors via the PC-canonical hybrid DSL
#'
#' The flexyBayes prior DSL (domain-specific language) lives on the
#' standard-deviation scale (never precision / variance) and accepts
#' two canonical idioms: distributional (`half_normal(scale = 1)`)
#' and tail-quantile / PC (penalised complexity)
#' (`pc(upper = 1, prob = 0.01)` meaning `Pr(sigma > 1) = 0.01`).
#' The PC idiom is the cross-engine interlingua -- it survives every
#' backend because it is a probability statement, not a
#' distributional name.
#'
#' v0.1 supports the targets and distributions listed above in the
#' file header. Calls outside the supported set raise a structured
#' error with the supported list.
#'
#' Argument matching is strict. Each distribution has a fixed parameter
#' list, and a call is matched against it by base R's own rules, so
#' `half_normal(1)` and `half_normal(scale = 1)` are the same
#' specification. An argument name the distribution does not have, a
#' duplicated argument, a missing required argument, a non-scalar or
#' non-finite value, or a value outside its domain (a non-positive
#' scale, a tail probability outside `(0, 1)`) is refused at
#' construction with a condition carrying a `flexybayes_refusal_*`
#' class -- see [fb_refusals()]. Both halves of a PC prior are
#' required: `pc(upper = 1)` states no probability, and is refused
#' rather than completed with a default the caller never wrote.
#'
#' @param ... One or more two-sided formulas of the form
#'   `target ~ distribution(args)`, one per parameter being given a
#'   prior.
#'
#'   Examples:
#'
#'   * `sigma ~ pc(upper = 2, prob = 0.05)`
#'   * `sd(group = "subject") ~ half_normal(scale = 1)`
#'   * `b("treatment") ~ student_t(df = 4, scale = 2.5)`
#'   * `cor(group = "subject") ~ lkj(eta = 2)`
#'   * `sd(group = "subject") ~ uniform(lower = 0, upper = 5)`
#'
#' Supported distribution families: `pc`, `half_normal`, `half_cauchy`,
#' `student_t`, `normal`, `exponential`, `lkj`, `cauchy`, `gamma`,
#' `uniform`. Note that `uniform()` on a variance component sits
#' outside the PC-canonical interlingua, but both active backends
#' represent it faithfully on the SD scale: the INLA backend via an
#' expression-prior on the log-precision, and brms as a bounded
#' `uniform()` prior on each simple random-effect (and residual) SD.
#'
#' @section What each backend can carry:
#'
#' `fb_prior()` is engine-neutral: it records what was asked for. Whether
#' a given (target, distribution) pair can be *carried* depends on the
#' backend, and the two differ.
#'
#' * **brms** takes any of the ten distributions on `sigma` and
#'   `sd(group = )` -- the parameter is bounded below at zero and brms
#'   renormalises a two-sided density over that support, which is why
#'   `normal(0, s)` and `half_normal(s)` are the same prior there -- and
#'   `normal` / `student_t` on `b()`. It carries no `cor()` or
#'   `smooth()` prior, because no model this package emits on brms has
#'   either parameter.
#' * **INLA** takes `half_normal`, `uniform`, `pc` and `exponential` on
#'   `sigma`, `sd(group = )` and `smooth()` (`exponential` is the PC
#'   prior written in its distributional form), and `normal` on `b()`,
#'   which arrives through `control.fixed`.
#'
#' A specification the selected backend cannot carry is refused at the
#' fit, naming the row and both remedies -- re-express it in a
#' distribution that engine does carry, or switch backend. Before 0.9.2
#' the INLA route discarded such a row silently while
#' [prior_summary()] printed it as applied, so a prior-sensitivity
#' analysis could vary a scale and compare the engine's default with
#' itself. `prior_summary()` lists each prior with the engine slot it
#' reached, as the engine reports it after fitting, and names each stated
#' or default prior that reached no slot, with the reason.
#'
#' A prior naming a coefficient or a variance component the model does
#' not have is likewise refused at the fit, against the engine's own
#' parameter list and in the vocabulary the caller wrote.
#'
#' @section The default prior:
#'
#' A fit that states neither `prior` nor `prior_vc_sd` runs under the
#' package default: `half_normal(scale = s)` on the residual standard
#' deviation and on the standard deviation of every variance component
#' the default reaches, on both engines. `s` is `5 * sd(y)` on an identity
#' link, `3` on a log link and `5` on a logit link; on an uncorrelated
#' random slope on a numeric covariate `x` it is `s / sd(x)`. The default
#' has no upper bound. [prior_summary()] prints the scale a fit ran under.
#'
#' @returns An `fb_prior` object, an S3 class inheriting from list,
#'   whose `$specs` element carries the parsed `target` and `spec`
#'   pairs.
#'
#' @examples
#' p <- fb_prior(
#'   sigma                 ~ pc(upper = 2, prob = 0.05),
#'   sd(group = "subject") ~ half_normal(scale = 1),
#'   b("treatment")        ~ student_t(df = 4, scale = 2.5)
#' )
#' p
#' @export
fb_prior <- function(...) {
  args <- list(...)
  if (length(args) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "prior_spec_empty",
      message = "`fb_prior()` requires at least one specification."
    ))
  }

  specs <- vector("list", length(args))
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (!inherits(a, "formula")) {
      stop(.fb_refusal_condition(
        reason_code = "prior_spec_not_formula",
        message = paste0(
          "Each `fb_prior()` argument must be a two-sided formula ",
          "`target ~ distribution(...)`. Argument ",
          i,
          " is: ",
          .fb_prior_deparse(a),
          "."
        )
      ))
    }
    if (length(a) != 3L) {
      stop(.fb_refusal_condition(
        reason_code = "prior_spec_not_two_sided",
        message = paste0(
          "`fb_prior()` formula ",
          i,
          " must be two-sided ",
          "(target on left, distribution on right). Got: ",
          .fb_prior_deparse(a),
          "."
        )
      ))
    }
    specs[[i]] <- list(
      target = .parse_prior_target(a[[2]]),
      spec = .parse_prior_distribution(a[[3]], envir = environment(a))
    )
  }

  obj <- list(specs = specs)
  class(obj) <- c("fb_prior", "list")
  obj
}

is_fb_prior <- function(x) inherits(x, "fb_prior")

# ---------------------------------------------------------------- #
# Parsing helpers                                                  #
# ---------------------------------------------------------------- #

# Parse target side of a prior spec formula:
#   sigma                 -> list(type = "sigma")
#   sd(group = "subject") -> list(type = "sd", group = "subject")
#   b("treatment")        -> list(type = "b", name = "treatment")
#   cor(group = "subject") -> list(type = "cor", group = "subject")
#   smooth("time")        -> list(type = "smooth", var = "time")
.parse_prior_target <- function(expr) {
  if (is.name(expr)) {
    nm <- as.character(expr)
    if (nm == "sigma") {
      return(list(type = "sigma"))
    }
    # A bare name other than `sigma` used to be carried as an opaque
    # target and then dropped at emit time, so `phi ~ half_normal(...)`
    # was accepted and never applied. There is no parameter it can name.
    .stop_prior_target_unsupported(expr)
  }
  if (is.call(expr)) {
    fn <- as.character(expr[[1]])
    args_list <- as.list(expr[-1])

    if (fn == "sigma") {
      return(list(type = "sigma"))
    }

    if (fn == "sd") {
      group <- .extract_string_arg(args_list, "group")
      if (is.null(group)) {
        stop(.fb_refusal_condition(
          reason_code = "prior_target_argument_missing",
          message = paste0(
            "sd() prior target requires `group = \"...\"` argument."
          )
        ))
      }
      return(list(type = "sd", group = group))
    }

    if (fn == "b") {
      if (length(args_list) < 1L) {
        stop(.fb_refusal_condition(
          reason_code = "prior_target_argument_missing",
          message = paste0(
            "b() prior target requires a name string, e.g., ",
            "b(\"treatment\")."
          )
        ))
      }
      # One coefficient, named by a string or a bare symbol. Anything else
      # used to reach as.character() on the unevaluated argument, which
      # deparses a call head-first: b(c("x", "z")) became the three names
      # c("c", "x", "z") and b(starts_with("x")) became c("starts_with", "x").
      # Neither was refused; the phantom head survived fb_plan() and the
      # length-3 target then failed inside the emitter with the base error
      # "the condition has length > 1" -- an R internals message about a
      # formula the user wrote. Vector and tidy-select targets are a real
      # feature and need a typed term resolver, not a coercion; until that
      # exists the refusal names what was written and what b() takes.
      raw <- args_list[[1]]
      ok <- is.symbol(raw) ||
        (is.character(raw) && length(raw) == 1L && !is.na(raw))
      if (!ok) {
        stop(.fb_refusal_condition(
          reason_code = "prior_target_not_a_single_name",
          message = paste0(
            "`b(", paste(deparse(raw), collapse = ""), ")` does not name one ",
            "coefficient. b() takes a single coefficient name, written as a ",
            "string -- b(\"treatment\") -- or as a bare symbol -- ",
            "b(treatment).",
            "\nA vector or a tidy-select helper is not yet supported: give ",
            "each coefficient its own row in fb_prior()."
          )
        ))
      }
      return(list(type = "b", name = as.character(raw)))
    }

    if (fn == "cor") {
      group <- .extract_string_arg(args_list, "group")
      if (is.null(group)) {
        stop(.fb_refusal_condition(
          reason_code = "prior_target_argument_missing",
          message = paste0(
            "cor() prior target requires `group = \"...\"` argument."
          )
        ))
      }
      return(list(type = "cor", group = group))
    }

    if (fn == "smooth") {
      if (length(args_list) < 1L) {
        stop(.fb_refusal_condition(
          reason_code = "prior_target_argument_missing",
          message = paste0(
            "smooth() prior target requires a variable name, e.g., ",
            "smooth(\"time\")."
          )
        ))
      }
      var <- as.character(args_list[[1]])
      basis <- if (!is.null(args_list$basis)) {
        as.character(args_list$basis)
      } else {
        "rw2"
      }
      return(list(type = "smooth", var = var, basis = basis))
    }
  }
  .stop_prior_target_unsupported(expr)
}

# Deparse a language object to one line -- a long call deparses to a
# character vector, and a refusal message must be a single string.
.fb_prior_deparse <- function(expr) {
  paste(deparse(expr), collapse = " ")
}

.stop_prior_target_unsupported <- function(expr) {
  stop(.fb_refusal_condition(
    reason_code = "prior_target_unsupported",
    message = paste0(
      "Unsupported prior target: ",
      .fb_prior_deparse(expr),
      ". Supported: sigma, sd(group = ...), b(\"name\"), ",
      "cor(group = ...), smooth(\"var\")."
    )
  ))
}

# Helper -- extract a named string argument from a parsed call.
.extract_string_arg <- function(args_list, name) {
  if (!is.null(args_list[[name]])) {
    return(as.character(args_list[[name]]))
  }
  unnamed <- args_list[
    !nzchar(names(args_list) %||% rep("", length(args_list)))
  ]
  if (length(unnamed) >= 1L) {
    return(as.character(unnamed[[1]]))
  }
  NULL
}

# Parse distribution side of a prior spec formula. Returns:
#   list(family = "pc", args = list(upper = 1, prob = 0.01))

# The parameter contract for every supported prior family -- the one
# place the DSL's argument names, its required arguments, and the domain
# each argument lives on are written down. `.name_prior_args()` matches a
# user call against `params` so a positional call binds exactly as R
# would, and `.check_prior_args()` enforces `required` and `domain`
# afterwards. The names match what the emit paths read (legacy reads `sd`
# for normal, `scale` for half_normal; INLA and brms read `upper` /
# `prob` for pc), so a named and a positional call are the same
# specification downstream.
#
# Before 0.9.2 this table held parameter names only, an unmatchable call
# was returned unchanged, and every emit-side read carried a `%||%`
# default -- so `half_normal(sd = 1.5)` was accepted and fitted under
# `half_normal(scale = 1)`. The required and domain columns exist so the
# refusal happens where the user can see it.
#
# Domains:
#   positive     finite and strictly greater than zero (scales, rates,
#                shapes, degrees of freedom, a PC upper bound)
#   nonnegative  finite and >= 0 (the uniform lower bound, which lives on
#                the non-negative SD scale)
#   real         finite, any sign (means and locations)
#   unit_open    finite and strictly inside (0, 1) (tail probabilities)
.fb_prior_family_table <- function() {
  list(
    normal = list(
      params = c("mean", "sd"),
      required = "sd",
      domain = c(mean = "real", sd = "positive")
    ),
    half_normal = list(
      params = "scale",
      required = "scale",
      domain = c(scale = "positive")
    ),
    half_cauchy = list(
      params = "scale",
      required = "scale",
      domain = c(scale = "positive")
    ),
    cauchy = list(
      params = c("location", "scale"),
      required = "scale",
      domain = c(location = "real", scale = "positive")
    ),
    student_t = list(
      params = c("df", "location", "scale"),
      required = c("df", "scale"),
      domain = c(df = "positive", location = "real", scale = "positive")
    ),
    exponential = list(
      params = "rate",
      required = "rate",
      domain = c(rate = "positive")
    ),
    gamma = list(
      params = c("shape", "rate"),
      required = c("shape", "rate"),
      domain = c(shape = "positive", rate = "positive")
    ),
    lkj = list(
      params = "eta",
      required = "eta",
      domain = c(eta = "positive")
    ),
    pc = list(
      params = c("upper", "prob"),
      required = c("upper", "prob"),
      domain = c(upper = "positive", prob = "unit_open")
    ),
    uniform = list(
      params = c("lower", "upper"),
      required = "upper",
      domain = c(lower = "nonnegative", upper = "positive")
    )
  )
}

# The supported distribution names, read off the one table rather than
# re-listed. A family without a parameter contract cannot be emitted.
.fb_prior_supported_families <- function() {
  names(.fb_prior_family_table())
}

# Rewrite a prior-distribution call so positional arguments carry their
# canonical names, using base R's own argument matching via match.call()
# against a stub whose formals are the family's parameters. An
# unmatchable call is a refusal, not a fallback: the argument names are
# the mini-language's public surface, and a name the family does not have
# is a specification the package cannot honour.
.name_prior_args <- function(expr, fn) {
  entry <- .fb_prior_family_table()[[fn]]
  if (is.null(entry)) {
    return(expr)
  }
  params <- entry$params
  stub <- function() NULL
  fm <- stats::setNames(rep(list(quote(expr = )), length(params)), params)
  formals(stub) <- fm
  matched <- tryCatch(
    match.call(stub, expr),
    error = function(e) conditionMessage(e)
  )
  if (is.character(matched)) {
    .stop_prior_call_unmatchable(fn, params, expr, matched)
  }
  matched
}

# Turn base R's own argument-matching failure into a typed refusal that
# names the family's parameters. Two shapes reach here: a duplicated
# formal, and a name (or a positional overflow) the family does not have.
.stop_prior_call_unmatchable <- function(fn, params, expr, detail) {
  supported <- paste0(fn, "(", paste(params, collapse = ", "), ")")
  if (grepl("matched by multiple", detail, fixed = TRUE)) {
    stop(.fb_refusal_condition(
      reason_code = "prior_argument_duplicated",
      message = paste0(
        "Prior specification ",
        .fb_prior_deparse(expr),
        " binds the same ",
        "hyperparameter more than once (",
        detail,
        "). Supply each ",
        "argument of ",
        supported,
        " once."
      )
    ))
  }
  stop(.fb_refusal_condition(
    reason_code = "prior_argument_unknown",
    message = paste0(
      "Prior specification ",
      .fb_prior_deparse(expr),
      " carries an ",
      "argument `",
      fn,
      "()` does not have (",
      detail,
      "). Supported: ",
      supported,
      "."
    )
  ))
}

# Enforce the family's required arguments and the domain of each one,
# after evaluation, on the canonically-named argument list.
.check_prior_args <- function(fn, args, expr) {
  entry <- .fb_prior_family_table()[[fn]]
  if (is.null(entry)) {
    return(invisible(NULL))
  }
  supported <- paste0(fn, "(", paste(entry$params, collapse = ", "), ")")
  nms <- names(args) %||% rep("", length(args))

  unknown <- setdiff(nms[nzchar(nms)], entry$params)
  if (length(unknown)) {
    stop(.fb_refusal_condition(
      reason_code = "prior_argument_unknown",
      message = paste0(
        "Prior specification ",
        .fb_prior_deparse(expr),
        " carries the ",
        "argument",
        if (length(unknown) > 1L) "s" else "",
        " ",
        paste0("`", unknown, "`", collapse = ", "),
        ", which `",
        fn,
        "()` does not have. Supported: ",
        supported,
        "."
      )
    ))
  }

  missing_args <- setdiff(entry$required, nms)
  if (length(missing_args)) {
    stop(.fb_refusal_condition(
      reason_code = "prior_argument_missing",
      message = paste0(
        "Prior specification ",
        .fb_prior_deparse(expr),
        " omits ",
        paste0("`", missing_args, "`", collapse = ", "),
        ", which `",
        fn,
        "()` requires. Supported: ",
        supported,
        if (identical(fn, "pc")) {
          paste0(
            ". A PC prior is the pair: `pc(upper = U, prob = p)` states ",
            "Pr(sigma > U) = p, and half of that statement fixes no prior."
          )
        } else {
          "."
        }
      )
    ))
  }

  for (nm in names(args)) {
    if (!nzchar(nm)) {
      next
    }
    .check_prior_hyperparameter(
      value = args[[nm]],
      name = nm,
      domain = entry$domain[[nm]] %||% "real",
      fn = fn,
      expr = expr
    )
  }
  invisible(NULL)
}

# One hyperparameter: a finite numeric scalar, then the family's domain.
.check_prior_hyperparameter <- function(value, name, domain, fn, expr) {
  label <- paste0("`", name, "` in ", .fb_prior_deparse(expr))
  if (
    length(value) == 1L &&
      !is.numeric(value) &&
      is.atomic(value) &&
      is.na(value)
  ) {
    stop(.fb_refusal_condition(
      reason_code = "prior_hyperparameter_not_scalar",
      message = paste0(label, " must be a finite number. Got NA.")
    ))
  }
  if (!is.numeric(value) || length(value) != 1L) {
    stop(.fb_refusal_condition(
      reason_code = "prior_hyperparameter_not_scalar",
      message = paste0(
        label,
        " must be a single number. Got ",
        if (is.numeric(value)) {
          paste0("a numeric vector of length ", length(value))
        } else if (is.language(value)) {
          paste0("the unevaluated expression `", .fb_prior_deparse(value), "`")
        } else {
          paste0("an object of class ", paste(class(value), collapse = "/"))
        },
        "."
      )
    ))
  }
  if (!is.finite(value)) {
    stop(.fb_refusal_condition(
      reason_code = "prior_hyperparameter_not_scalar",
      message = paste0(
        label,
        " must be a finite number. Got ",
        format(value),
        "."
      )
    ))
  }
  ok <- switch(
    domain,
    "positive" = value > 0,
    "nonnegative" = value >= 0,
    "unit_open" = value > 0 && value < 1,
    "real" = TRUE,
    TRUE
  )
  if (isTRUE(ok)) {
    return(invisible(NULL))
  }
  requirement <- switch(
    domain,
    "positive" = paste0(
      "be strictly positive -- zero or a negative value names no ",
      "distribution on this parameter"
    ),
    "nonnegative" = "be >= 0 (the SD scale is non-negative)",
    "unit_open" = "be a probability strictly inside (0, 1)",
    "lie in its domain"
  )
  stop(.fb_refusal_condition(
    reason_code = "prior_hyperparameter_out_of_domain",
    message = paste0(
      label,
      " must ",
      requirement,
      ". Got ",
      format(value),
      "."
    )
  ))
}

.parse_prior_distribution <- function(expr, envir = baseenv()) {
  if (!is.call(expr)) {
    stop(.fb_refusal_condition(
      reason_code = "prior_distribution_not_a_call",
      message = paste0(
        "Prior distribution must be a call (e.g., ",
        "`pc(upper = 1, prob = 0.01)`). Got: ",
        .fb_prior_deparse(expr),
        "."
      )
    ))
  }

  fn <- as.character(expr[[1]])
  supported <- .fb_prior_supported_families()

  if (!fn %in% supported) {
    stop(.fb_refusal_condition(
      reason_code = "prior_distribution_unknown",
      message = paste0(
        "Unsupported prior distribution: ",
        fn,
        ". Supported distributions: ",
        paste(supported, collapse = ", "),
        "."
      )
    ))
  }

  # Name positional arguments by the family's canonical parameter order
  # *before* evaluating, so that `normal(0, 50)` and `normal(0, sd = 50)`
  # carry identical names downstream. Every emit path (INLA, brms)
  # reads these arguments by name; without this step a
  # positional call silently drops to the default scale on the by-name
  # paths. Matching is delegated to base R via match.call() on a stub
  # with the canonical formals, so named and positional binding follow
  # the usual rules -- and a call base R cannot match is refused here
  # rather than carried unnamed into the emit paths.
  expr <- .name_prior_args(expr, fn)

  # Evaluate distribution arguments in the formula's calling
  # environment so users can pass data-driven scales like
  # `pc(upper = 2.5 * sd(y), prob = 0.05)`. An expression that does not
  # evaluate is left in place and refused by the scalar check below,
  # which names the argument.
  args <- as.list(expr[-1])
  args <- lapply(args, function(a) {
    tryCatch(eval(a, envir = envir), error = function(e) a)
  })

  # uniform carries one argument the table cannot express: `lower`
  # defaults to 0, and `upper` must exceed it. The default is applied
  # before the shared check so a bare `uniform(upper = 5)` is complete.
  if (fn == "uniform" && is.null(args$lower)) {
    args$lower <- 0
  }

  # Required arguments, scalar-ness, and the per-argument domain.
  .check_prior_args(fn, args, expr)

  if (fn == "uniform" && args$upper <= args$lower) {
    stop(.fb_refusal_condition(
      reason_code = "prior_hyperparameter_out_of_domain",
      message = paste0(
        "uniform() upper must be > lower. Got lower = ",
        args$lower,
        ", upper = ",
        args$upper,
        "."
      )
    ))
  }

  list(family = fn, args = args)
}

# Local null-coalescing op (R < 4.4 compat -- DESCRIPTION targets >= 4.1).
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------- #
# Default-prior construction (Simpson-2017 PC) -- internal         #
# ---------------------------------------------------------------- #

# Pick a sensible scale U for the PC default from the response and
# (family, link). Identity/Gaussian: U = 2.5 * sd(y). Log-link
# (Poisson, NegBin, Gamma): U = 2.5 * sd(log(y + 0.5)) so the scale
# matches the random effect's working scale. Logit-link (binomial,
# beta): U = 2.5 (logit-scale random effects are typically O(1) so
# response-derived scale is meaningless). Returns (scale, basis)
# where basis records which branch fired -- captured so the
# deprecation message can explain itself.
.default_prior_scale <- function(data, response, family, link = NULL) {
  fam <- if (inherits(family, "family")) {
    family$family
  } else {
    as.character(family)
  }
  fam <- tolower(fam %||% "gaussian")
  link_name <- if (!is.null(link)) {
    tolower(link)
  } else if (inherits(family, "family")) {
    tolower(family$link)
  } else {
    NULL
  }
  is_log_link <- isTRUE(link_name == "log") ||
    fam %in%
      c(
        "poisson",
        "negative_binomial",
        "negbinom",
        "nbinomial",
        "gamma",
        "hurdle_gamma"
      )
  is_logit_link <- isTRUE(link_name == "logit") ||
    fam %in% c("binomial", "binary", "beta")

  if (is_logit_link) {
    return(list(scale = 2.5, basis = "logit_default"))
  }

  y <- if (is.data.frame(data)) {
    .fb_response_values(list(response = response), data)
  } else {
    NULL
  }

  if (is_log_link) {
    if (!is.numeric(y)) {
      return(list(scale = 1.0, basis = "log_default_nonnumeric"))
    }
    yp <- as.numeric(y)
    yp <- yp[is.finite(yp) & yp >= 0]
    if (length(yp) < 2L) {
      return(list(scale = 1.0, basis = "log_default_too_small"))
    }
    sdv <- stats::sd(log(yp + 0.5))
    if (!is.finite(sdv) || sdv <= 0) {
      sdv <- 1.0
    }
    return(list(scale = 2.5 * sdv, basis = "log_link_sd"))
  }

  if (!is.numeric(y)) {
    return(list(scale = 1.0, basis = "identity_default"))
  }
  yv <- as.numeric(y)
  yv <- yv[is.finite(yv)]
  if (length(yv) < 2L) {
    return(list(scale = 1.0, basis = "identity_too_small"))
  }
  sdv <- stats::sd(yv)
  if (!is.finite(sdv) || sdv <= 0) {
    sdv <- 1.0
  }
  list(scale = 2.5 * sdv, basis = "identity_sd")
}

# Construct a PC prior fb_prior() for a model. PC spec applied to
# sigma + every named random group in `random_groups`. Retained in
# v0.1.x as an explicit-choice constructor (e.g., used in the advanced
# priors vignette appendix and PC-default tests); no longer the v0.1
# default. The bounded-uniform default below supersedes the former
# PC default.
.default_pc_prior <- function(
  data,
  response,
  family,
  link = NULL,
  random_groups = character(0),
  alpha = 0.05
) {
  scl <- .default_prior_scale(data, response, family, link)
  specs <- list(list(
    target = list(type = "sigma"),
    spec = list(family = "pc", args = list(upper = scl$scale, prob = alpha))
  ))
  for (g in unique(random_groups)) {
    if (!nzchar(g)) {
      next
    }
    specs[[length(specs) + 1L]] <- list(
      target = list(type = "sd", group = g),
      spec = list(family = "pc", args = list(upper = scl$scale, prob = alpha))
    )
  }
  obj <- list(specs = specs)
  class(obj) <- c("fb_prior", "list")
  attr(obj, "fb_prior_default_basis") <- scl$basis
  attr(obj, "fb_prior_default_scale") <- scl$scale
  obj
}

# The default variance-component prior is half_normal(scale = s) on each
# standard deviation, with s = U / .FB_DEFAULT_SCALE_DIVISOR, where U is the
# bound the uniform default of 0.9.0 to 0.12.x used: 5 sd(y) on an identity
# link, 3 on a log link, 5 on a logit link. The uniform was replaced because
# INLA's reported interval for a weakly identified component can run past a
# bound, which the support audit refuses (agridat::john.alpha,
# `rep + rep:block`: 8 of 8 repeats); a half-normal has no bound to cross.
# The divisor is the "matched to the current bound" rule adopted in unit
# U-PRIOR-HN-1 against the anchored validation tests.
.FB_DEFAULT_SCALE_DIVISOR <- 1

# The multiplier that turns the response SD into the default's scale on the
# identity scale: s = .FB_DEFAULT_SD_MULTIPLIER * sd(y). It is named because
# two other places depend on the *relationship* and not just the number: with
# a heterogeneous residual the scalar sigma prior is retargeted onto the
# log-sigma coefficients (.brms_retarget_sigma_for_heterogeneous_residual()
# in R/priors_to_brms.R), and the INLA Wishart scale (.fb_inla_response_var()
# in R/emit_inla.R) is written on the response variance; both recover sd(y)
# by dividing s by this multiplier. Written out twice, the constants drifted
# apart once, and the retarget inherited the superseded PC default's 2.5
# while the uniform default had moved to 5.
.FB_DEFAULT_SD_MULTIPLIER <- 5 / .FB_DEFAULT_SCALE_DIVISOR

# Pick the half-normal scale s of the default per (family, link). Returns
# (scale, basis):
# - Gaussian / identity link: s = 5 * sd(y) / d; basis =
#   "identity_sd_half_normal"
# - Log link (Poisson, NegBin, Gamma): s = 3 / d; basis = "log_link_half_normal"
#   (3 on the log scale is a 20x ratio of group means).
# - Logit link (Binomial, Beta): s = 5 / d; basis = "logit_half_normal"
# where d is .FB_DEFAULT_SCALE_DIVISOR. The 5 * sd(y) is a flexyBayes
# heuristic, wide enough to be weakly informative; Gelman (2006, Bayesian
# Analysis 1(3):515-534) discusses half-t and half-Cauchy alternatives for
# very small J, available via fb_prior().
#
# `resp_sd`, when given, is the response standard deviation read from
# somewhere other than `data` -- the streamed route holds only a head sample
# of the rows and reads it off the full sufficient statistics instead.
.default_half_normal_scale <- function(data, response, family, link = NULL,
                                       resp_sd = NULL) {
  link_class <- .fb_default_link_class(family, link)

  if (identical(link_class, "logit")) {
    return(list(scale = 5 / .FB_DEFAULT_SCALE_DIVISOR,
                basis = "logit_half_normal"))
  }

  if (identical(link_class, "log")) {
    return(list(scale = 3 / .FB_DEFAULT_SCALE_DIVISOR,
                basis = "log_link_half_normal"))
  }

  if (!is.null(resp_sd)) {
    sdv <- resp_sd
    if (!is.finite(sdv) || sdv <= 0) {
      sdv <- 1.0
    }
    return(list(
      scale = .FB_DEFAULT_SD_MULTIPLIER * sdv,
      basis = "identity_sd_half_normal"
    ))
  }

  # A transformed response, log(y) or I(10 * y), is not a column: it is read
  # through the same resolver the rest of the fit uses, so the scale is
  # written on the quantity the residual is written on.
  y <- if (is.data.frame(data)) {
    .fb_response_values(list(response = response), data)
  } else {
    NULL
  }

  if (!is.numeric(y)) {
    return(list(scale = 1.0, basis = "identity_default_half_normal"))
  }
  yv <- as.numeric(y)
  yv <- yv[is.finite(yv)]
  if (length(yv) < 2L) {
    return(list(scale = 1.0, basis = "identity_too_small_half_normal"))
  }
  sdv <- stats::sd(yv)
  if (!is.finite(sdv) || sdv <= 0) {
    sdv <- 1.0
  }
  list(
    scale = .FB_DEFAULT_SD_MULTIPLIER * sdv,
    basis = "identity_sd_half_normal"
  )
}

# Which scale a default prior is written on: "logit" and "log" for the links
# whose linear predictor is unit-free, "identity" for everything whose
# linear predictor carries the response's units. Shared by the
# variance-component default above and the fixed-effect default on INLA
# (R/inla_default_priors.R), so the two cannot classify a family differently.
.fb_default_link_class <- function(family, link = NULL) {
  fam <- if (inherits(family, "family")) {
    family$family
  } else {
    as.character(family)
  }
  fam <- tolower(fam %||% "gaussian")
  link_name <- if (!is.null(link)) {
    tolower(link)
  } else if (inherits(family, "family")) {
    tolower(family$link)
  } else {
    NULL
  }
  is_log_link <- isTRUE(link_name == "log") ||
    fam %in%
      c(
        "poisson",
        "negative_binomial",
        "negbinom",
        "nbinomial",
        "gamma",
        "hurdle_gamma"
      )
  is_logit_link <- isTRUE(link_name == "logit") ||
    fam %in% c("binomial", "binary", "beta")

  if (is_logit_link) {
    return("logit")
  }
  if (is_log_link) {
    return("log")
  }
  "identity"
}

# Construct the default variance-component prior fb_prior() for a model:
# sigma ~ half_normal(scale = s) and every named random group's
# SD ~ half_normal(scale = s), with s from .default_half_normal_scale().
#
# The default also fires on vm() and ped() structured-cov random terms.
# The two forms share the simple-RI scale derivation (5 * sd(response)
# heuristic; identity / log / logit family-aware branches via
# .default_half_normal_scale()) because their natural variance-matrix
# scaling is identity in the SD-on-the-Cholesky-scaled-random-effect
# interpretation that codegen .code_random's vm / ped branches use
# (u_<tag> = L_G %*% (z * sigma_<tag>); the Cholesky absorbs G's scale and
# sigma is purely the residual-scale multiplier). At 0.9.0 the nested
# interaction intercept, the heterogeneous-variance diag / idh / at group
# and the unstructured us group joined them -- see
# .fb_default_prior_targets() for which term types reach the default and
# which are left to the engine.
#
# `vm_ped_groups` is a character vector of `term$var` strings (one
# per vm / ped random term) and is keyed identically to `random_groups`
# in priors_to_legacy() / priors_to_inla() so the per-component
# entries surface through the existing dispatch.
.default_half_normal_prior <- function(
  data,
  response,
  family,
  link = NULL,
  random_groups = character(0),
  vm_ped_groups = character(0),
  resp_sd = NULL
) {
  scl <- .default_half_normal_scale(
    data, response, family, link, resp_sd = resp_sd
  )
  specs <- list(list(
    target = list(type = "sigma"),
    spec = list(family = "half_normal", args = list(scale = scl$scale))
  ))
  for (g in unique(random_groups)) {
    if (!nzchar(g)) {
      next
    }
    specs[[length(specs) + 1L]] <- list(
      target = list(type = "sd", group = g),
      spec = list(family = "half_normal", args = list(scale = scl$scale))
    )
  }
  # vm() + ped() branch. Add a separate spec per structured group
  # so prior_summary() surfaces them distinctly from the simple-RI
  # groups; carry a `_default_form` tag (vm / ped) so the
  # surface can explain which structured form the spec came from.
  for (g in unique(vm_ped_groups)) {
    if (!nzchar(g)) {
      next
    }
    sp <- list(
      target = list(type = "sd", group = g),
      spec = list(family = "half_normal", args = list(scale = scl$scale))
    )
    attr(sp, "_default_form") <- "vm_or_ped"
    specs[[length(specs) + 1L]] <- sp
  }
  obj <- list(specs = specs)
  class(obj) <- c("fb_prior", "list")
  attr(obj, "fb_prior_default_basis") <- scl$basis
  attr(obj, "fb_prior_default_scale") <- scl$scale
  if (length(vm_ped_groups) > 0L) {
    attr(obj, "fb_prior_default_vm_ped_groups") <- unique(vm_ped_groups)
  }
  obj
}

# .fb_default_prior_targets() --- which random terms the shared default
# reaches, and which are left to whatever the engine does on its own.
#
# One walk over the IR's random terms, consumed twice: flexybayes() reads
# `$shared` / `$vm_ped` to build the default half-normal prior, and the
# model fingerprint (R/model_fingerprint.R) reads `$engine` to record the
# parameters that carry no prior this package chose. Both consumers must
# agree about which is which -- a term whose prior is the shared default in
# one place and an engine default in the other would make the fit's prior
# record, which prior_summary() prints, disagree with the prior the fit
# carried, so the two views are derived here rather than written out twice.
#
# `$shared` and `$vm_ped` are group names in the key the prior DSL uses
# (`sd(group = <name>)`); `$engine` is a data.frame of *canonical parameter
# names* (the brms-shaped names of R/canonical_names.R) with the reason
# each is left alone.
#
# The uncorrelated-slope form contributes two entries, the intercept
# variance under the grouping factor's own name and the slope variance
# under the synthesised `<slope_var>_<group>` key that
# .priors_to_brms_specs() unwraps back into brms's (class = "sd",
# coef = <slope_var>, group = <group>) row.
.fb_default_prior_targets <- function(fb) {
  shared <- character(0)
  vm_ped <- character(0)
  eng_param <- character(0)
  eng_reason <- character(0)

  add_engine <- function(param, reason) {
    eng_param <<- c(eng_param, param)
    eng_reason <<- c(eng_reason, reason)
  }

  for (t in fb$random_terms %||% list()) {
    ttype <- t$type %||% "<unknown>"

    if (ttype %in% c("simple", "ide", "id")) {
      if (!is.null(t$var) && nzchar(t$var)) {
        shared <- c(shared, t$var)
      }
      next
    }

    if (identical(ttype, "simple_slope_uncor")) {
      if (isTRUE(t$with_intercept) && !is.null(t$var) && nzchar(t$var)) {
        shared <- c(shared, t$var)
      }
      if (!is.null(t$slope_var) && nzchar(t$slope_var)) {
        shared <- c(shared, paste0(t$slope_var, "_", t$var))
      }
      next
    }

    # (x | g) -- the same two standard deviations, plus the correlation
    # between them. brms takes an SD prior on both; the correlation keeps
    # brms's LKJ, which the prior DSL does not carry.
    if (identical(ttype, "simple_slope_cor")) {
      if (!is.null(t$var) && nzchar(t$var)) {
        shared <- c(shared, t$var)
      }
      if (!is.null(t$slope_var) && nzchar(t$slope_var)) {
        shared <- c(shared, paste0(t$slope_var, "_", t$var))
        add_engine(
          paste0("cor_", t$var),
          paste0(
            "(", t$slope_var, " | ", t$var, ") intercept-slope correlation ",
            "keeps brms's default LKJ prior; flexyBayes sets no correlation ",
            "prior"
          )
        )
      }
      next
    }

    if (ttype %in% c("vm", "ped")) {
      if (!is.null(t$var) && nzchar(t$var)) {
        vm_ped <- c(vm_ped, t$var)
      }
      next
    }

    # Interaction random intercept, A:B. brms groups by the interaction
    # itself, so the group name is the colon-joined pair and one SD carries
    # the whole term.
    if (identical(ttype, "nested")) {
      grp <- paste0(t$outer, ":", t$inner)
      shared <- c(shared, grp)
      next
    }

    # diag(f):g / idh(f):g / at(f):g -- one free SD per level of the outer
    # factor, all of them in brms's `sd` class for the inner group, so a
    # single group-keyed spec covers every level.
    if (identical(ttype, "at_simple")) {
      if (!is.null(t$inner) && nzchar(t$inner)) {
        shared <- c(shared, t$inner)
      }
      next
    }

    # us(f):g -- the same per-level SDs, plus a correlation block. The SDs
    # take the shared default; the correlations keep brms's own LKJ, which
    # has no counterpart in the prior DSL and no INLA analogue.
    # fa(f, k):g -- the common factors and the per-level specific variances
    # are ordinary iid fields keyed by the inner group, so the shared SD
    # default reaches them. It has to: left on INLA's default log-gamma
    # precision, the specific variance at the level with the largest loading
    # collapses to zero (measured: a true 0.25 came back at 0.007). The
    # LOADINGS are copy-scaling parameters with no SD prior to set, so they
    # keep the engine's own and are recorded as such.
    if (identical(ttype, "fa_gxe")) {
      if (!is.null(t$inner) && nzchar(t$inner)) {
        shared <- c(shared, t$inner)
        add_engine(
          paste0("loadings_", t$inner),
          paste0(
            "fa(", t$outer, ", ", t$k, "):", t$inner,
            " loadings are INLA copy-scaling parameters; flexyBayes sets ",
            "no prior on them"
          )
        )
      }
      next
    }

    if (identical(ttype, "us_gxe")) {
      if (!is.null(t$inner) && nzchar(t$inner)) {
        shared <- c(shared, t$inner)
        add_engine(
          paste0("cor_", t$inner),
          paste0(
            "us(",
            t$outer,
            "):",
            t$inner,
            " level correlations keep brms's ",
            "default LKJ prior; flexyBayes sets no correlation prior"
          )
        )
      }
      next
    }

    # Everything below is outside the walker. Recording the parameter is
    # the point: an unrecorded prior and a matched one are indistinguishable
    # downstream, and the comparison then runs on terms the two engines were
    # never asked the same question about.
    if (identical(ttype, "combo")) {
      add_engine(
        paste0("sd_", paste(as.character(t$vars), collapse = ":")),
        "multi-way random interaction is outside the default-prior walker"
      )
      next
    }

    if (ttype %in% c("ar1", "ar1_spatial")) {
      add_engine(
        "sd_spatial",
        "AR1 latent field keeps INLA's own default hyperprior"
      )
      add_engine(
        "rho_row",
        "AR1 latent field keeps INLA's own default hyperprior"
      )
      if (identical(ttype, "ar1_spatial")) {
        add_engine(
          "rho_col",
          "AR1 latent field keeps INLA's own default hyperprior"
        )
      }
      next
    }

    lbl <- t$var %||% t$label %||% ttype
    add_engine(
      paste0("sd_", lbl),
      paste0(
        "random term type \"",
        ttype,
        "\" is outside the default-prior walker"
      )
    )
  }

  list(
    shared = unique(shared),
    vm_ped = unique(vm_ped),
    engine = data.frame(
      parameter = eng_param,
      reason = eng_reason,
      stringsAsFactors = FALSE
    )
  )
}

# One-time announcement of the default variance-component prior.
# Silenceable via options(flexyBayes.silence_default_prior_note = TRUE).
# The once-flag is held in the package-internal emit-state env (see
# R/emit_state.R); v0.3.9 migrated the flag out of the options() namespace
# where unrelated callers could consume it before the intended emission.
.default_prior_note_once <- function(scale) {
  if (isTRUE(getOption("flexyBayes.silence_default_prior_note", FALSE))) {
    return(invisible())
  }
  if (.emit_state_get("default_prior_note")) {
    return(invisible())
  }
  message(
    "flexyBayes: by default the residual and every random-effect group ",
    "have a half_normal(scale = ",
    format(scale, digits = 3),
    ") prior on their standard deviation, with the scale divided by sd(x) ",
    "for a random slope on x. See ?fb_prior to set other priors; ",
    "options(flexyBayes.silence_default_prior_note = TRUE) hides this note."
  )
  .emit_state_set("default_prior_note", TRUE)
}

# Faithful INLA expression prior for a uniform(lower, upper) prior on a
# variance component on the SD scale. INLA parameterises the precision
# hyperparameter internally as theta = log(precision), so sigma =
# exp(-theta / 2). Transforming the flat density p(sigma) = 1 / (U - L)
# on (L, U) through that change of variables gives
#
#   log p(theta) = -log(U - L) - log(2) - theta / 2
#
# on the support L < sigma < U, i.e. -2 log(U) <= theta <= -2 log(L). The
# upper-theta bound vanishes when L = 0 (sigma may approach 0). Outside the
# support a large negative log-density (-1e10, INLA's finite stand-in for
# -Inf) enforces the bound. The string is emitted on a single line so it
# splices unmodified through .inla_hyper_arg() and survives the
# formula-as-text -> as.formula() -> INLA::inla() round trip.
#
# This is the exact representation of the package's uniform-on-SD default
# (and any user-supplied uniform() prior), replacing the former lossy
# PC-prior approximation that concentrated mass at sigma = 0 and so
# disagreed with brms's flat uniform on small group counts.
.inla_uniform_sd_expr <- function(upper, lower = 0) {
  if (lower <= 0) {
    sprintf(
      paste0(
        "expression: U=%.16g; lb=-2*log(U);",
        " ld=-log(U)-log(2)-theta/2;",
        " return( theta<lb ? -1.0e10 : ld );"
      ),
      upper
    )
  } else {
    sprintf(
      paste0(
        "expression: L=%.16g; U=%.16g; lo=-2*log(U); hi=-2*log(L);",
        " ld=-log(U-L)-log(2)-theta/2;",
        " return( (theta<lo)||(theta>hi) ? -1.0e10 : ld );"
      ),
      lower,
      upper
    )
  }
}

# Faithful INLA expression prior for a half_normal(scale = s) prior on a
# variance component on the SD scale. The half-normal density
# p(sigma) = (2 / (s sqrt(2 pi))) exp(-sigma^2 / (2 s^2)) on (0, Inf),
# transformed through sigma = exp(-theta / 2) (Jacobian (1/2)
# exp(-theta / 2)), gives
#
#   log p(theta) = -0.5 log(2 pi) - log(s) - sigma^2 / (2 s^2) - theta / 2
#
# (the half-normal's log(2) cancels the Jacobian's -log(2)). Single line
# for the same splicing reason as .inla_uniform_sd_expr(). Replaces the
# former PC-prior approximation of half_normal().
.inla_halfnormal_sd_expr <- function(scale) {
  sprintf(
    paste0(
      "expression: s=%.16g; sig=exp(-theta/2);",
      " return( -0.5*log(2*3.141592653589793) - log(s)",
      " - (sig*sig)/(2*s*s) - theta/2 );"
    ),
    scale
  )
}

# Faithful INLA expression prior for the legacy `lognormal(0, s)`
# variance-component prior on the SD scale -- the brms path emits this
# as a `lognormal` row for `prior_vc_sd`. INLA parameterises the
# hyperparameter as theta = log(precision), so
# sigma = exp(-theta / 2) and log sigma = -theta / 2. Transforming
# p(sigma) = 1 / (sigma s sqrt(2 pi)) exp(-(log sigma)^2 / (2 s^2))
# through that change of variables (Jacobian sigma / 2) gives
#
#   log p(theta) = -log(2) - log(s) - 0.5 log(2 pi) - theta^2 / (8 s^2)
#
# on the whole real line, so the prior needs no support guard. Single line
# for the same splicing reason as the two above.
#
# Before 0.9.2 the INLA route consumed no legacy scalar at all: a fit
# passed `prior_vc_sd = 3` ran under INLA's own log-gamma precision
# default while `prior_summary()` printed the lognormal in its header.
#
# A location mu on log sigma, which the standardised INLA fit needs
# (R/inla_standardise.R), generalises the last term to
# -(theta / 2 + mu)^2 / (2 s^2); at mu = 0 the expression is unchanged.
.inla_lognormal_sd_expr <- function(scale, location = 0) {
  if (identical(location, 0)) {
    return(sprintf(
      paste0(
        "expression: s=%.16g;",
        " return( -log(2.0) - log(s) - 0.5*log(2*3.141592653589793)",
        " - (theta*theta)/(8.0*s*s) );"
      ),
      scale
    ))
  }
  sprintf(
    paste0(
      "expression: s=%.16g; mu=%.16g; u=theta/2.0+mu;",
      " return( -log(2.0) - log(s) - 0.5*log(2*3.141592653589793)",
      " - (u*u)/(2.0*s*s) );"
    ),
    scale,
    location
  )
}

# The legacy scalar variance-component prior, keyed the way
# priors_to_inla() keys an fb_prior: "sigma" for the residual, the
# grouping-factor tag for each variance component the default-prior
# walker reaches. One walk, so the two routes cannot disagree about
# which components carry a prior this package chose.
.priors_legacy_to_inla <- function(fb, vc_sd) {
  targets <- .fb_default_prior_targets(fb)
  keys <- unique(c("sigma", targets$shared, targets$vm_ped))
  keys <- keys[nzchar(keys)]
  entry <- list(
    prior = .inla_lognormal_sd_expr(vc_sd),
    meta = list(family = "lognormal", scale = vc_sd)
  )
  stats::setNames(rep(list(entry), length(keys)), keys)
}

# Did the caller write this scalar prior argument? flexybayes() records
# the answer on the IR; a direct emit_*() call carries no record, and the
# absence means "not supplied", which is what it meant before 0.9.2.
.fb_prior_scalar_supplied <- function(fb, which) {
  flags <- fb$prior_scalars$supplied
  if (is.null(flags) || is.null(flags[[which]])) {
    return(FALSE)
  }
  isTRUE(unname(flags[[which]]))
}

# The value of a scalar prior argument as the caller left it, for the
# surfaces that report what the fit ran under. NA when the fit carries no
# record (a direct emit_*() call).
.fb_prior_scalar_value <- function(fb, which) {
  value <- fb$prior_scalars[[which]]
  if (is.null(value)) NA_real_ else value
}

# ---------------------------------------------------------------- #
# Print method                                                     #
# ---------------------------------------------------------------- #

#' Print method for fb_prior
#'
#' @param x   An `fb_prior` object as returned by [fb_prior()].
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the specification list
#'   it prints.
#' @keywords internal
#' @export
print.fb_prior <- function(x, ...) {
  cat(
    "<fb_prior> ",
    length(x$specs),
    " specification",
    if (length(x$specs) != 1L) "s" else "",
    "\n",
    sep = ""
  )
  for (s in x$specs) {
    cat(
      "  ",
      .format_prior_target(s$target),
      " ~ ",
      .format_prior_distribution(s$spec),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

.format_prior_target <- function(target) {
  switch(
    target$type,
    "sigma" = "sigma",
    "sd" = paste0("sd(group = \"", target$group, "\")"),
    "b" = paste0("b(\"", target$name, "\")"),
    "cor" = paste0("cor(group = \"", target$group, "\")"),
    "smooth" = paste0(
      "smooth(\"",
      target$var,
      "\", basis = \"",
      target$basis,
      "\")"
    ),
    "name" = target$name,
    deparse(target)
  )
}

.format_prior_distribution <- function(spec) {
  arg_strs <- vapply(
    seq_along(spec$args),
    function(i) {
      nm <- names(spec$args)[i]
      val <- spec$args[[i]]
      if (is.numeric(val) && length(val) == 1L) {
        return(paste0(
          if (nzchar(nm %||% "")) paste0(nm, " = ") else "",
          format(val)
        ))
      }
      paste0(if (nzchar(nm %||% "")) paste0(nm, " = ") else "", deparse(val))
    },
    character(1)
  )
  paste0(spec$family, "(", paste(arg_strs, collapse = ", "), ")")
}

# ---------------------------------------------------------------- #
# Translation helpers (cross-engine emit hooks)                    #
# ---------------------------------------------------------------- #

# Translate fb_prior -> INLA hyperpar control list: a named list keyed
# by f()-term group / smoother variable, with "sigma" for the residual,
# mapping to the `hyper` body INLA takes on that term.
#
# Untranslatable rows refuse here rather than falling through. Before
# 0.9.2 this function handled `pc`, `half_normal` and `uniform` on
# `sigma` / `sd()` / `smooth()` and returned nothing at all for every
# other pair, while `prior_summary()` printed each of them as applied --
# a prior-sensitivity analysis on INLA could vary a `student_t` scale
# and compare the engine's default with itself (field-sweep FS-21).
# Which pairs translate is declared once, in
# `.fb_prior_translation_table()` (R/prior_translation.R).
#
# `b()` rows are translatable but do not arrive here: INLA states the
# fixed-effect prior in `control.fixed`, which
# `.priors_to_inla_control_fixed()` builds. They are skipped rather than
# refused.
priors_to_inla <- function(prior) {
  out <- list()
  if (!inherits(prior, "fb_prior")) {
    return(out)
  }
  .fb_check_prior_translation(prior, "inla")

  for (s in prior$specs) {
    if (!s$target$type %in% c("sigma", "sd", "smooth")) {
      # b() arrives through control.fixed; nothing else survives the
      # translation check above.
      next
    }
    key <- switch(
      s$target$type,
      "sigma" = "sigma",
      "sd" = s$target$group,
      "smooth" = s$target$var
    )
    if (s$spec$family == "pc") {
      u <- as.numeric(s$spec$args$upper %||% 1)
      a <- as.numeric(s$spec$args$prob %||% 0.01)
      out[[key]] <- list(prior = "pc.prec", param = c(u, a))
    } else if (s$spec$family == "exponential") {
      # The PC prior IS exponential on the SD scale: `pc.prec(U, alpha)`
      # states P(sigma > U) = alpha, i.e. Exp(-log(alpha) / U). Writing
      # U = 1 / rate and alpha = exp(-1) recovers the requested rate
      # exactly for any positive rate, and avoids the underflow that
      # U = 1, alpha = exp(-rate) would hit for a large rate.
      rate <- as.numeric(s$spec$args$rate)
      out[[key]] <- list(
        prior = "pc.prec",
        param = c(1 / rate, exp(-1)),
        meta = list(family = "exponential", rate = rate)
      )
    } else if (s$spec$family == "half_normal") {
      scale <- as.numeric(s$spec$args$scale %||% 1)
      # Exact half_normal(scale = s) via an INLA expression prior on the
      # SD scale (no longer the lossy PC approximation).
      out[[key]] <- list(
        prior = .inla_halfnormal_sd_expr(scale),
        meta = list(family = "half_normal", scale = scale)
      )
    } else if (s$spec$family == "uniform") {
      lower <- as.numeric(s$spec$args$lower %||% 0)
      upper <- as.numeric(s$spec$args$upper)
      # Exact uniform(lower, upper) on the SD scale via an INLA
      # expression prior (.inla_uniform_sd_expr). This supersedes the
      # former lossy PC-prior approximation: the PC prior concentrated
      # mass at sigma = 0 and so shrank a small-group variance component
      # more than brms's flat uniform did, producing a
      # cross-engine prior mismatch that surfaced as spurious
      # cross-engine disagreement on the variance components. The
      # expression prior is the faithful representation, so the two
      # engines now carry genuinely the same default prior and there is
      # no approximation to flag.
      out[[key]] <- list(
        prior = .inla_uniform_sd_expr(upper = upper, lower = lower),
        meta = list(family = "uniform", lower = lower, upper = upper)
      )
    }
  }
  out
}

# Translate the `b()` rows of an fb_prior into INLA's control.fixed.
#
# INLA states the fixed-effect prior as a mean and a precision, and both
# accept per-coefficient named lists with a `default` fallback (verified
# against INLA 25.10.19: `control.fixed = list(mean = list(x = 1,
# default = 0), prec = list(x = 0.0625, default = 0.001))` round-trips
# into `fit$.args$control.fixed`). The intercept has its own pair of
# arguments, `mean.intercept` / `prec.intercept`.
#
# `available` is the design-matrix column names the emit will build, so
# a `b()` row naming a coefficient the model does not carry refuses by
# name (field-sweep FS-20) instead of arriving as a `control.fixed`
# entry INLA silently ignores.
.priors_to_inla_control_fixed <- function(prior, available = NULL) {
  out <- list(mean = list(), prec = list())
  intercept <- list()
  if (!inherits(prior, "fb_prior")) {
    return(list())
  }
  for (s in prior$specs) {
    if (!identical(s$target$type, "b")) {
      next
    }
    nm <- s$target$name
    mean_v <- as.numeric(.named_or_positional(s$spec$args, "mean", 1L, 0))
    sd_v <- as.numeric(.named_or_positional(s$spec$args, "sd", 2L, 1))
    if (nm %in% c("Intercept", "(Intercept)")) {
      intercept <- list(
        mean.intercept = mean_v,
        prec.intercept = 1 / (sd_v^2)
      )
      next
    }
    if (!is.null(available) && !nm %in% available) {
      .fb_stop_prior_target_absent(
        target_label = paste0("b(\"", nm, "\")"),
        kind = "fixed-effect coefficient",
        available = available,
        engine = "inla"
      )
    }
    out$mean[[nm]] <- mean_v
    out$prec[[nm]] <- 1 / (sd_v^2)
  }
  if (!length(out$mean)) {
    out$mean <- NULL
    out$prec <- NULL
  }
  c(out[!vapply(out, is.null, logical(1))], intercept)
}
