# fb_from_brms: flexyBayes brms-format ingest path.
#
# Parses a brms-style formula
#   `response ~ fixed_effects + (1 | g1) + (1 | g2) + ...`
# via R's formula abstract syntax tree (AST), with no runtime brms
# dependency. brmsformula objects (built via brms::bf()) are accepted
# by unwrapping `formula$formula` when brms is loaded.
#
# Supported:
# - The families .resolve_family() admits, as a name or a family
#   object, with the default link, the object's own link, or an
#   explicit override.
# - Linear fixed effects (factor / continuous / interaction /
#   I() expressions), parsed by .parse_fixed as on the asreml
#   ingest path.
# - Random intercepts via `(1 | g)` and `(1 || g)`.
# - Uncorrelated random intercept + slope via `(x || g)` and the
#   equivalent `(1 + x || g)`. Each level of g carries an
#   independent intercept deviation AND an independent slope
#   deviation, both drawn from their own univariate-normal priors
#   (no correlation parameter). The explicit `(0 + x || g)` shape
#   suppresses the intercept block (slope-only). Semantics match
#   lme4 / brms exactly: `(x || g)` expands to
#   `(1 | g) + (0 + x | g)`.
# - Optional `weights` argument mapped to addition_terms with
#   type = "weights" (mirrors fb_from_asreml).
#
# Refused at ingest, with a message naming the term:
# - Smoothers `s()`, `t2()`.
# - Gaussian processes `gp()`.
# - Autocorrelation `ar()`, `ma()`, `arma()`, `cosy()`, `car()`.
# - Distributional regression and addition terms (`cens`, `trunc`,
#   `trials`, `se`, `mi`, `weights()` inside the formula). The
#   `weights` argument is supported instead.
# - Multivariate / mixture / categorical families, and every other
#   family .resolve_family() does not admit (refused by that gate
#   here, for a family name and a family object alike).

#' Ingest a brms-format formula into the flexyBayes IR
#'
#' Parses a brms / lme4-style two-sided formula
#' (`response ~ fixed + (1 | g)`) into a `fb_terms` object --
#' flexyBayes's backend-agnostic intermediate representation (IR). The IR
#' is what every engine emits from, so building it explicitly lets a
#' power user inspect the parsed model, cache it, or hand it to a fitting
#' verb. Most users never call this directly: [fb()] and [flexybayes()]
#' detect a brms-style `formula` argument and build the IR internally.
#'
#' Linear fixed effects, random intercepts (`(1 | g)`, `(1 || g)`) and
#' uncorrelated random intercept+slope (`(x || g)`) are supported;
#' correlated random slopes (`(x | g)`), smoothers, Gaussian processes
#' and autocorrelation terms refuse at ingest with a structured message.
#'
#' @param formula A base-R formula (recommended) or a `brmsformula`
#'   object. `brmsformula` support requires \pkg{brms} to be installed.
#' @param data A data.frame containing every referenced variable. May be
#'   `NULL` only on the advanced metadata-only path (see `carry_n_rows`).
#' @param family Character family name or a base-R `family()` object.
#' @param link Character link override, or `NULL` for the family default.
#' @param prior An [fb_prior()] object, or `NULL`.
#' @param weights Optional numeric weights vector of length `nrow(data)`,
#'   mapped to a single addition term with `type = "weights"`.
#' @param prior_fixed_sd Numeric SD for fixed-effect normal priors when
#'   `prior` is `NULL`.
#' @param prior_vc_sd Numeric hyperparameter for the variance-component
#'   priors when `prior` is `NULL`.
#' @param carry_n_rows Advanced: a positive integer enabling the
#'   metadata-only IR path (`data = NULL`). The formula's variables are
#'   realised as a one-row stand-in data frame and the row count is recorded as
#'   `carry_n_rows`, for stress-testing the preflight layer without
#'   materialising the full data. Leave `NULL` for ordinary use.
#' @param ... Reserved for future brms-ingest options (specials, autocor
#'   handling); currently unused.
#'
#' @return An `fb_terms` object with `source = "brms"`.
#'
#' @section Deprecated:
#' `fb_from_brms()` is deprecated as of 0.13.0. Use [fb_plan()] to inspect
#' the backend-neutral model plan without fitting it.
#'
#' @family flexyBayes ingest adapters
#' @seealso [fb()] and [flexybayes()] for the universal fitting entry;
#'   [fb_from_asreml()] for the other active dialect.
#' @examples
#' set.seed(1L)
#' df <- data.frame(
#'   y = rnorm(30),
#'   x = rnorm(30),
#'   g = factor(rep(letters[1:5], 6))
#' )
#' plan <- fb_plan(y ~ x + (1 | g), data = df)
#' class(plan)
#' @export
fb_from_brms <- function(
  formula,
  data,
  family = "gaussian",
  link = NULL,
  prior = NULL,
  weights = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  carry_n_rows = NULL,
  ...
) {
  .fb_warn_deprecated_export(
    "fb_from_brms",
    "Use `fb_plan()` to inspect the compiled model plan."
  )
  .fb_from_brms_impl(
    formula = formula,
    data = data,
    family = family,
    link = link,
    prior = prior,
    weights = weights,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    carry_n_rows = carry_n_rows,
    ...
  )
}

.fb_from_brms_impl <- function(
  formula,
  data,
  family = "gaussian",
  link = NULL,
  prior = NULL,
  weights = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  carry_n_rows = NULL,
  ...
) {
  if (inherits(formula, "brmsformula")) {
    .check_installed(
      "brms",
      "Package 'brms' is required to ingest brmsformula objects. ",
      "Install via: install.packages(\"brms\")"
    )
    .fb_brms_refuse_unfitted_parts(formula, family, link)
    formula <- formula$formula
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula or a brmsformula object.", call. = FALSE)
  }
  if (length(formula) < 3L) {
    stop("`formula` must be two-sided: response ~ predictors.", call. = FALSE)
  }

  # Metadata-only path: data = NULL + carry_n_rows lets
  # the caller build an IR for stress-testing the preflight layer
  # without materialising N rows of data. The formula's free variables
  # are realised as a one-row stand-in data.frame (all numeric).
  # Random-term level counts are stripped from the IR so .fb_preflight()
  # reads them from the <fb_dataset>'s dictionaries instead.
  if (is.null(data)) {
    if (is.null(carry_n_rows)) {
      stop(
        "`data` is NULL but `carry_n_rows` was not supplied. ",
        "Pass `carry_n_rows = N` for the metadata-only IR path, ",
        "or pass a real data.frame.",
        call. = FALSE
      )
    }
    if (
      !is.numeric(carry_n_rows) ||
        length(carry_n_rows) != 1L ||
        carry_n_rows < 1
    ) {
      stop(
        "`carry_n_rows` must be a positive numeric scalar; got: ",
        deparse(carry_n_rows),
        call. = FALSE
      )
    }
    data <- .fb_from_brms_standin_data(formula)
    fb <- .fb_from_brms_impl(
      formula = formula,
      data = data,
      family = family,
      link = link,
      prior = prior,
      weights = NULL, # weights cannot be metadata-only
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      carry_n_rows = NULL, # break the recursion
      ...
    )
    fb$data_summary$n <- as.integer(carry_n_rows)
    # Preflight reads level counts from <fb_dataset>$dictionaries when
    # var_n is NA.
    for (i in seq_along(fb$random_terms)) {
      fb$random_terms[[i]]$var_n <- NA_integer_
      fb$random_terms[[i]]$var_levels <- NULL
    }
    for (i in seq_along(fb$fixed_terms)) {
      if (identical(fb$fixed_terms[[i]]$type, "factor")) {
        fb$fixed_terms[[i]]$n_levels <- NA_integer_
        fb$fixed_terms[[i]]$levels <- NULL
      }
    }
    fb$source <- "brms_metadata_only"
    return(fb)
  }

  # When `data` is present and `carry_n_rows` is set, override the
  # cached row count on the returned IR (use case: stress-testing the
  # preflight ceiling on a synthetic claim that the data is bigger
  # than nrow(data)).
  override_n_rows <- carry_n_rows

  # LHS-pipe brms addition forms (`y | trials(n) ~ ...`,
  # `y | cens(c) ~ ...`, `y | weights(w) ~ ...`, `y | se(s) ~ ...`).
  # The walker does not parse addition forms. Without this refusal the
  # later names(data) lookup fails with "Response variable not found
  # in data".
  if (
    is.call(formula[[2L]]) &&
      identical(formula[[2L]][[1L]], as.name("|"))
  ) {
    stop(
      "brms ingest does not yet support LHS addition forms ",
      "(`y | trials(n) ~ ...`, `y | cens(c) ~ ...`, etc.). The ",
      "supported brms corpus is single-column responses only ",
      "(use `family = \"binomial\"` for Bernoulli with a 0/1 ",
      "response). Addition forms are queued for a future ",
      "expansion of the brms ingest layer. Got: ",
      deparse(formula[[2L]]),
      call. = FALSE
    )
  }

  # ------------------------- Family + link ----------------------- #
  # A family object goes through the same gate as its name. Reading
  # `family$family` straight off the object let this grammar fit families
  # the string spelling refuses (brms::lognormal() on either engine) and
  # send others to an untyped emitter error (brms::student()).
  # .resolve_family() reads the object's name, keeps the link the object
  # supplies, and refuses a `link =` that contradicts it.
  fl <- .resolve_family(family, link)
  fam_name <- fl$family
  link_name <- fl$link

  # ------------------------- Response --------------------------- #
  # A transformed response, log(y), is evaluated by .parse_fixed() below,
  # as the ASReml grammar does.
  response <- paste(deparse(formula[[2]]), collapse = "")

  # ------------------------- Walk RHS --------------------------- #
  walked <- .brms_walk_rhs(formula[[3]])

  # Refuse unsupported features at ingest.
  if (length(walked$ext_terms) > 0L) {
    # A univariate smoother gets the smoother refusal, which names
    # `random = ~ spl(x)` on INLA. The multivariate forms (`te`, `ti`,
    # `t2`, `gp`, `sos`) have no `spl()` counterpart and get the
    # ingest-feature refusal below.
    for (t in walked$ext_terms) {
      lbl <- paste(t$deparse, collapse = " ")
      head <- sub("\\(.*$", "", lbl)
      if (head %in% c("s", "spl")) {
        .stop_fixed_smoother_unsupported(lbl, head, data)
      }
    }
    classes <- vapply(
      walked$ext_terms,
      function(t) {
        # deparse() returns a multi-line character vector for a
        # sufficiently long call (e.g. s(x, k = 10, representation =
        # list(scheme = "low_rank_smooth", rank = 5L))); collapse to a
        # single string so this fail-fast refusal fires cleanly rather
        # than tripping the vapply length-1 constraint.
        paste0(t$type, " (", paste(t$deparse, collapse = " "), ")")
      },
      character(1)
    )
    stop(.fb_refusal_condition(
      reason_code = "brms_ingest_feature_unsupported",
      message = paste0(
        "The brms ingest layer has no lowering for: ",
        paste(classes, collapse = ", "),
        ".\n",
        "These features are queued for a later expansion of the ingest ",
        "layer."
      ),
      family_class = "flexybayes_brms_ingest_refusal"
    ))
  }

  # Build a fixed-only formula and re-parse via .parse_fixed (gives
  # us byte-identical fixed_terms output to the asreml ingest path)
  fixed_only_rhs <- if (length(walked$fixed_labels)) {
    paste(walked$fixed_labels, collapse = " + ")
  } else {
    "1"
  }
  fixed_form <- stats::as.formula(
    paste0(response, " ~ ", fixed_only_rhs)
  )
  fixed_info <- .parse_fixed(fixed_form, data)

  # Translate each (1 | g) / (x || g) / (1 + x || g) RE pair to a
  # parse_formula.R-style term descriptor. `(x | g)` correlated slopes
  # raise a structured refusal with a deferral pointer;
  # everything else outside the supported set raises the legacy
  # refusal.
  random_terms <- list()
  for (pair in walked$re_pairs) {
    rt <- .brms_re_to_fb(pair, data)
    if (identical(rt$type, "brms_re_factor_slope_unsupported")) {
      .stop_factor_slope_unsupported(
        grouping_factor = rt$rhs,
        slope_variable = rt$slope_var,
        correlated = rt$correlated,
        found = rt$deparse
      )
    }
    if (identical(rt$type, "brms_re_unsupported")) {
      stop(.fb_refusal_condition(
        reason_code = "brms_random_effect_form_unsupported",
        message = paste0(
          "The bar grammar has no lowering for this random-effect ",
          "specification.\n",
          "Found: ",
          rt$deparse,
          "\n",
          "Supported bar forms: (1 | g), (1 || g), (x || g), ",
          "(1 + x || g).\n",
          "Multi-variable slopes and the group-of-coefficients shorthand ",
          "are queued for\n  a later expansion of the brms ingest layer, ",
          "and have no ASReml surface that\n  expresses the same model ",
          "today."
        ),
        family_class = "flexybayes_brms_ingest_refusal"
      ))
    }
    random_terms[[length(random_terms) + 1L]] <- rt
  }

  # ------------------------- Weights --------------------------- #
  addition_terms <- if (!is.null(weights)) {
    .fb_check_weights(weights, data)
    list(list(type = "weights", values = weights))
  } else {
    list()
  }

  # ------------------------- Priors ---------------------------- #
  priors <- if (!is.null(prior)) {
    prior
  } else {
    list(fixed_sd = prior_fixed_sd, vc_sd = prior_vc_sd, legacy = TRUE)
  }

  data_summary <- list(
    n = nrow(data),
    known_matrices = character(0),
    response_sd = .fb_response_sd(
      list(response = response, offsets = fixed_info$offsets),
      data
    )
  )

  if (!is.null(override_n_rows)) {
    data_summary$n <- as.integer(override_n_rows)
  }

  new_fb_terms(
    response = response,
    response_is_expression = isTRUE(fixed_info$response_is_expression),
    family = fam_name,
    link = link_name,
    intercept = fixed_info$intercept,
    fixed_terms = fixed_info$terms,
    random_terms = random_terms,
    residual_terms = list(), # brms folds residual into family
    addition_terms = addition_terms,
    priors = priors,
    data_summary = data_summary,
    capabilities = character(),
    source = "brms",
    offsets = fixed_info$offsets
  )
}


# .fb_brms_refuse_unfitted_parts(): refuse the parts of a brms::bf()
# object that the ingest does not carry to an engine.
#
# The ingest reads `$formula` alone, so a non-linear declaration, a
# formula or a fixed value for a distributional parameter, and a family
# named inside bf() each refuse by name. `family` and `link` are the
# fit's own, and a family inside bf() that agrees with them passes.
# `caller` prefixes the message.
#
# @noRd
# @keywords internal
.fb_brms_refuse_unfitted_parts <- function(formula, family, link, caller = "") {
  refuse <- function(what) {
    stop(.fb_refusal_condition(
      reason_code = "brms_ingest_feature_unsupported",
      message = paste0(caller, what),
      family_class = "flexybayes_brms_ingest_refusal"
    ))
  }
  shown <- function(x) paste(deparse(x), collapse = " ")

  if (isTRUE(attr(formula$formula, "nl"))) {
    refuse(paste0(
      "brms::bf(nl = TRUE) declares a non-linear model (",
      shown(formula$formula),
      "). flexyBayes fits predictors that are ",
      "linear in their terms, so it refuses this model."
    ))
  }
  pforms <- formula$pforms %||% list()
  if (length(pforms) > 0L) {
    refuse(paste0(
      "brms::bf() gives a distributional parameter its own formula (",
      paste(vapply(pforms, shown, character(1L)), collapse = ", "),
      "). flexyBayes fits one value for each distributional parameter and ",
      "carries no such formula to either engine, so it refuses this model."
    ))
  }
  pfix <- formula$pfix %||% list()
  if (length(pfix) > 0L) {
    refuse(paste0(
      "brms::bf() fixes a distributional parameter (",
      paste(
        names(pfix),
        vapply(pfix, shown, character(1L)),
        sep = " = ",
        collapse = ", "
      ),
      "). flexyBayes estimates every distributional parameter and carries ",
      "no fixed value to either engine, so it refuses this model."
    ))
  }
  if (!is.null(formula$family)) {
    named <- .resolve_family(formula$family, NULL)
    fitted <- .resolve_family(family, link)
    if (
      !identical(c(named$family, named$link), c(fitted$family, fitted$link))
    ) {
      refuse(paste0(
        "brms::bf() names the family ",
        named$family,
        " (",
        named$link,
        " link), and the `family` argument gives ",
        fitted$family,
        " (",
        fitted$link,
        " link). flexyBayes reads the family from the ",
        "`family` argument, so it refuses this model. Pass the family ",
        "there and leave it out of brms::bf()."
      ))
    }
  }
  invisible(NULL)
}


# Synthesise a one-row, all-numeric stand-in
# data.frame from a brms-style formula. Used by the data = NULL +
# carry_n_rows ingest path so the brms walker can produce the IR
# shape. The realised level counts are stripped after the walk, and
# the preflight layer reads them from <fb_dataset>.
.fb_from_brms_standin_data <- function(formula) {
  all_vars <- all.vars(formula)
  if (!length(all_vars)) {
    stop(
      "formula references no variables; cannot synthesise ",
      "stand-in data.",
      call. = FALSE
    )
  }
  pl <- as.list(rep(1.0, length(all_vars)))
  names(pl) <- all_vars
  as.data.frame(pl, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------- #
# AST walking helpers                                              #
# ---------------------------------------------------------------- #

# Walk a brms-style formula RHS and bucket terms.
#
# Returns:
#   list(
#     fixed_labels = character vector of deparsed top-level terms,
#     re_pairs     = list of (lhs, rhs, cor) triplets for `(... | ...)`,
#     ext_terms    = list of fb_term descriptors for unsupported
#                    features (smoothers, GP, autocorrelation)
#   )
.brms_walk_rhs <- function(expr) {
  fixed_labels <- character(0)
  re_pairs <- list()
  ext_terms <- list()

  visit <- function(e) {
    # Strip parentheses around `(... | ...)` patterns first
    if (is.call(e) && identical(e[[1]], as.name("("))) {
      inner <- e[[2]]
      if (is.call(inner)) {
        op <- as.character(inner[[1]])
        if (op == "|") {
          re_pairs[[length(re_pairs) + 1L]] <<- list(
            lhs = inner[[2]],
            rhs = inner[[3]],
            cor = TRUE
          )
          return(invisible())
        }
        if (op == "||") {
          re_pairs[[length(re_pairs) + 1L]] <<- list(
            lhs = inner[[2]],
            rhs = inner[[3]],
            cor = FALSE
          )
          return(invisible())
        }
      }
      # Bare parentheses around a non-RE expression: descend
      visit(e[[2]])
      return(invisible())
    }

    # Recurse on `+` and `-` operators
    if (is.call(e)) {
      op <- as.character(e[[1]])
      if (op %in% c("+", "-")) {
        visit(e[[2]])
        if (length(e) >= 3L) {
          visit(e[[3]])
        }
        return(invisible())
      }

      # Unsupported brms specials go into ext_terms
      if (op %in% c("s", "t2")) {
        ext_terms[[length(ext_terms) + 1L]] <<- list(
          type = "smoother",
          fn = op,
          deparse = deparse(e)
        )
        return(invisible())
      }
      if (op == "gp") {
        ext_terms[[length(ext_terms) + 1L]] <<- list(
          type = "gp",
          deparse = deparse(e)
        )
        return(invisible())
      }
      if (op %in% c("ar", "ma", "arma", "cosy", "car")) {
        ext_terms[[length(ext_terms) + 1L]] <<- list(
          type = "autocorrelation",
          fn = op,
          deparse = deparse(e)
        )
        return(invisible())
      }
      if (
        op %in%
          c("cens", "trunc", "trials", "se", "mi", "weights", "subset", "rate")
      ) {
        ext_terms[[length(ext_terms) + 1L]] <<- list(
          type = "addition_form",
          fn = op,
          deparse = deparse(e)
        )
        return(invisible())
      }
    }

    # Default: a fixed-effect term, re-parsed later by .parse_fixed
    fixed_labels <<- c(fixed_labels, deparse(e))
    invisible()
  }

  visit(expr)

  list(
    fixed_labels = fixed_labels,
    re_pairs = re_pairs,
    ext_terms = ext_terms
  )
}

# Translate a single brms RE pair (lhs | rhs) to a parse_formula.R
# term descriptor.
#
# Supported:
#   (1 | g)        -> type = "simple"             (random intercept)
#   (1 || g)       -> type = "simple"             (random intercept)
#   (x || g)       -> type = "simple_slope_uncor" (intercept + slope,
#                                                  uncorrelated, with the
#                                                  intercept implicit as
#                                                  in lme4 / brms)
#   (1 + x || g)   -> type = "simple_slope_uncor" (with_intercept = TRUE;
#                                                  same model as above)
#   (0 + x || g)   -> type = "simple_slope_uncor" (with_intercept = FALSE;
#                                                  slope-only)
#   (x | g)        -> type = "simple_slope_cor" (intercept and slope share
#                     one 2-by-2 unstructured covariance)
#
# Anything else returns a marker descriptor whose type =
# "brms_re_unsupported" so the caller can reject with a clear
# message.
.brms_re_to_fb <- function(pair, data) {
  lhs <- pair$lhs
  rhs <- pair$rhs
  cor <- pair$cor

  group_name <- deparse(rhs)
  lhs_str <- deparse(lhs)

  is_one <- function(e) {
    if (
      is.numeric(e) && length(e) == 1L && isTRUE(all.equal(as.numeric(e), 1))
    ) {
      return(TRUE)
    }
    if (is.name(e) && as.character(e) == "1") {
      return(TRUE)
    }
    identical(deparse(e), "1")
  }

  is_simple_var <- function(e) {
    if (!is.name(e)) {
      return(FALSE)
    }
    nm <- as.character(e)
    nm %in%
      names(data) &&
      !is.factor(data[[nm]]) &&
      !is.character(data[[nm]]) &&
      is.numeric(data[[nm]])
  }

  # Random intercept: (1 | g) or (1 || g).
  if (is_one(lhs)) {
    # A grouping expression, not a column. `(1 | g:h)` is the bar
    # spelling of the ASReml `random = ~ g:h`, so it goes through the
    # same walker and gets the same descriptor, and the same gates.
    if (is.call(rhs) && identical(rhs[[1L]], as.name(":"))) {
      walked_group <- .walk(rhs)
      if (length(walked_group) == 1L) {
        return(.enrich(walked_group[[1L]], data = data))
      }
    }
    out <- list(type = "simple", var = group_name)
    if (group_name %in% names(data)) {
      f <- factor(data[[group_name]])
      out$var_n <- nlevels(f)
      out$var_levels <- levels(f)
    }
    return(out)
  }

  # Intercept + slope uncorrelated: (x || g). The bare-slope form is
  # sugar for (1 | g) + (0 + x | g), as in ?lme4::lFormula. An
  # explicit `0 +` suppresses the intercept (handled below).
  if (isFALSE(cor) && is_simple_var(lhs)) {
    return(.brms_re_to_fb_slope_uncor(
      slope_var = as.character(lhs),
      group_name = group_name,
      with_intercept = TRUE,
      data = data
    ))
  }

  # `lhs ± slope` forms: (1 + x || g), (0 + x || g), (x + 1 || g).
  # The `1` token signals the intercept block, and the `0` token
  # suppresses it. Anything else with the same shape, such as
  # (x + z || g) for two slopes, falls through to the generic
  # "brms_re_unsupported" refusal.
  is_zero <- function(e) {
    if (
      is.numeric(e) && length(e) == 1L && isTRUE(all.equal(as.numeric(e), 0))
    ) {
      return(TRUE)
    }
    if (is.name(e) && as.character(e) == "0") {
      return(TRUE)
    }
    identical(deparse(e), "0")
  }
  if (
    isFALSE(cor) &&
      is.call(lhs) &&
      identical(lhs[[1L]], as.name("+")) &&
      length(lhs) == 3L
  ) {
    left <- lhs[[2L]]
    right <- lhs[[3L]]
    sv <- NULL
    with_int <- NA
    if (is_one(left) && is_simple_var(right)) {
      sv <- as.character(right)
      with_int <- TRUE
    }
    if (is_one(right) && is_simple_var(left)) {
      sv <- as.character(left)
      with_int <- TRUE
    }
    if (is_zero(left) && is_simple_var(right)) {
      sv <- as.character(right)
      with_int <- FALSE
    }
    if (is_zero(right) && is_simple_var(left)) {
      sv <- as.character(left)
      with_int <- FALSE
    }
    if (!is.null(sv)) {
      return(.brms_re_to_fb_slope_uncor(
        slope_var = sv,
        group_name = group_name,
        with_intercept = with_int,
        data = data
      ))
    }
  }

  # Random slope on a factor: (f | g), (0 + f | g), (f || g). The bar
  # grammar has no lowering for it, so the refusal names the ASReml
  # spelling that fits: `random = ~ us(f):g` for the correlated form and
  # `random = ~ diag(f):g` for the uncorrelated one.
  factor_slope <- .brms_factor_slope_var(lhs, data)
  if (!is.na(factor_slope)) {
    return(list(
      type = "brms_re_factor_slope_unsupported",
      lhs = lhs_str,
      rhs = group_name,
      slope_var = factor_slope,
      correlated = isTRUE(cor),
      deparse = paste0(
        "(",
        lhs_str,
        if (isTRUE(cor)) " | " else " || ",
        group_name,
        ")"
      )
    ))
  }

  # Slope correlated: (x | g) or (1 + x | g). The intercept and the slope
  # share one 2-by-2 unstructured covariance.
  if (isTRUE(cor)) {
    slope_var <- NA_character_
    if (is_simple_var(lhs)) {
      slope_var <- as.character(lhs)
    } else if (
      is.call(lhs) && identical(lhs[[1L]], as.name("+")) && length(lhs) == 3L
    ) {
      left <- lhs[[2L]]
      right <- lhs[[3L]]
      if (is_one(left) && is_simple_var(right)) {
        slope_var <- as.character(right)
      }
      if (is_one(right) && is_simple_var(left)) slope_var <- as.character(left)
    }
    if (!is.na(slope_var)) {
      return(.brms_re_to_fb_slope_cor(
        slope_var = slope_var,
        group_name = group_name,
        data = data
      ))
    }
  }

  list(
    type = "brms_re_unsupported",
    lhs = lhs_str,
    rhs = group_name,
    cor = cor,
    deparse = paste0(
      "(",
      lhs_str,
      if (isTRUE(cor)) " | " else " || ",
      group_name,
      ")"
    )
  )
}

# Build an IR slot for an uncorrelated random-slope term. Used by
# both (x || g) and (1 + x || g). Carries the canonical-name pieces
# (slope_var, var) that codegen + canonical_names + emit_brms read.
.brms_re_to_fb_slope_uncor <- function(
  slope_var,
  group_name,
  with_intercept,
  data
) {
  out <- list(
    type = "simple_slope_uncor",
    var = group_name,
    slope_var = slope_var,
    with_intercept = isTRUE(with_intercept)
  )
  if (group_name %in% names(data)) {
    f <- factor(data[[group_name]])
    out$var_n <- nlevels(f)
    out$var_levels <- levels(f)
  }
  if (slope_var %in% names(data)) {
    out$slope_var_values <- as.numeric(data[[slope_var]])
  }
  out
}

# The correlated random slope `(x | g)`, carried as one term. The
# intercept and the slope deviations of a group share a 2-by-2
# unstructured covariance: two standard deviations and one correlation.
# brms fits it as `(x | g)` verbatim, INLA as an iid2d field over a
# two-block index.
.brms_re_to_fb_slope_cor <- function(slope_var, group_name, data) {
  out <- list(
    type = "simple_slope_cor",
    var = group_name,
    slope_var = slope_var,
    with_intercept = TRUE
  )
  if (group_name %in% names(data)) {
    f <- factor(data[[group_name]])
    out$var_n <- nlevels(f)
    out$var_levels <- levels(f)
  }
  if (slope_var %in% names(data)) {
    out$slope_var_values <- as.numeric(data[[slope_var]])
  }
  out
}

# The factor a bar-grammar random slope is written on, or NA when the
# left-hand side is not a factor slope. Recognises the bare form (`f`),
# the intercept-suppressed form (`0 + f`) and the intercept-carrying form
# (`1 + f`) in either operand order. All three have the same ASReml
# spelling.
.brms_factor_slope_var <- function(lhs, data) {
  is_factor_var <- function(e) {
    if (!is.name(e)) {
      return(FALSE)
    }
    nm <- as.character(e)
    nm %in% names(data) && (is.factor(data[[nm]]) || is.character(data[[nm]]))
  }
  if (is_factor_var(lhs)) {
    return(as.character(lhs))
  }
  if (is.call(lhs) && identical(lhs[[1L]], as.name("+")) && length(lhs) == 3L) {
    left <- lhs[[2L]]
    right <- lhs[[3L]]
    is_const <- function(e) {
      identical(deparse(e), "0") || identical(deparse(e), "1")
    }
    if (is_const(left) && is_factor_var(right)) {
      return(as.character(right))
    }
    if (is_const(right) && is_factor_var(left)) {
      return(as.character(left))
    }
  }
  NA_character_
}

# Raise the structured refusal for a bar-grammar factor slope, naming
# the ASReml surface that fits.
#
# The correlated form `(f | g)` is one free variance per level of `f`
# plus the full correlation block between levels: `us(f):g`. The
# uncorrelated form `(f || g)` drops the correlations and keeps the
# per-level variances: `diag(f):g`. Both fit on brms. On INLA `diag`
# fits and `us` refuses, and the message says so.
.stop_factor_slope_unsupported <- function(
  grouping_factor,
  slope_variable,
  correlated,
  found
) {
  asreml <- paste0(
    "random = ~ ",
    if (isTRUE(correlated)) "us(" else "diag(",
    slope_variable,
    "):",
    grouping_factor
  )
  other <- paste0(
    "random = ~ ",
    if (isTRUE(correlated)) "diag(" else "us(",
    slope_variable,
    "):",
    grouping_factor
  )
  msg <- paste0(
    "The bar grammar has no lowering for a random slope on a factor.\n",
    "Found: ",
    found,
    "   (`",
    slope_variable,
    "` is a factor)\n\n",
    "The same model in the ASReml grammar fits on brms:\n\n",
    "  flexybayes(<response> ~ ",
    slope_variable,
    ", ",
    asreml,
    ", data = <data>,\n",
    "             backend = \"brms\")\n\n",
    if (isTRUE(correlated)) {
      paste0(
        "`us(",
        slope_variable,
        "):",
        grouping_factor,
        "` is one free variance per level of `",
        slope_variable,
        "` within each\n  `",
        grouping_factor,
        "`, plus the correlations between those levels, which is what ",
        "the\n  correlated bar form asks for. Drop the correlations with ",
        "`",
        other,
        "`.\n"
      )
    } else {
      paste0(
        "`diag(",
        slope_variable,
        "):",
        grouping_factor,
        "` is one free variance per level of `",
        slope_variable,
        "` within each\n  `",
        grouping_factor,
        "`, with no correlation between levels, which is what the ",
        "uncorrelated\n  bar form asks for. Add the correlations with `",
        other,
        "`.\n"
      )
    },
    "On INLA, diag(",
    slope_variable,
    "):",
    grouping_factor,
    " fits and ",
    "us(",
    slope_variable,
    "):",
    grouping_factor,
    " is refused:\n  ",
    "flexyBayes fits estimated correlation matrices on brms only."
  )
  stop(.fb_refusal_condition(
    reason_code = "brms_factor_random_slope_unsupported",
    message = msg,
    family_class = "flexybayes_brms_ingest_refusal",
    grouping_factor = grouping_factor,
    slope_variable = slope_variable,
    correlated = isTRUE(correlated),
    asreml_surface = asreml
  ))
}
