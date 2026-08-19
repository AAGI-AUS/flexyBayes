# fb_from_brms -- flexyBayes brms-format ingest path.
#
# v0.1 MINIMUM SUBSET -- parses a brms-style formula
#   `response ~ fixed_effects + (1 | g1) + (1 | g2) + ...`
# via R's formula abstract syntax tree (AST). No runtime brms
# dependency for this minimum subset; brmsformula objects (built via
# brms::bf()) are accepted by unwrapping `formula$formula` when brms
# is loaded.
#
# Supported (v0.1):
# - Gaussian / binomial / poisson / Gamma / lognormal / beta family
#   names; default link or explicit override.
# - Linear fixed effects (factor / continuous / interaction /
#   I() expressions) -- re-parsed via the existing .parse_fixed
#   helper for byte-identity with the asreml ingest path.
# - Random intercepts via `(1 | g)` and `(1 || g)`.
# - Uncorrelated random intercept + slope via `(x || g)` and the
#   equivalent `(1 + x || g)`. Each level of g carries an
#   independent intercept deviation AND an independent slope
#   deviation, both drawn from their own univariate-normal priors
#   (no correlation parameter). The explicit `(0 + x || g)` shape
#   suppresses the intercept block (slope-only). Semantics match
#   lme4 / brms exactly: `(x || g)` expands to
#   `(1 | g) + (0 + x | g)`. Lifts the v0.2 refusal.
# - Optional `weights` argument mapped to addition_terms with
#   type = "weights" (mirrors fb_from_asreml).
#
# Refused (v0.1) -- fail-fast at ingest with a clear message:
# - Correlated random slopes `(x | g)` (deferred to the v0.3
#   structured-covariance representation; refuses with
#   condition class flexybayes_correlated_slope_unsupported).
# - Smoothers `s()`, `t2()`.
# - Gaussian processes `gp()`.
# - Autocorrelation `ar()`, `ma()`, `arma()`, `cosy()`, `car()`.
# - Distributional regression and addition terms (`cens`, `trunc`,
#   `trials`, `se`, `mi`, `weights()` inside formula; the
#   `weights` argument is supported instead).
# - Multivariate / hurdle / mixture / categorical families (caught
#   by lgm_gate post-construction or by the family-name check
#   here).
#
# Future iterations will use brms::brmsterms() for the full
# v0.1 ingest set and add full random-slope translation,
# smoothers via INLA's rw1/rw2, GP via SPDE, distributional
# regression refusal hooks, and the addition-form family.

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
#'   realised as a one-row placeholder and the row count is recorded as
#'   `carry_n_rows`, for stress-testing the preflight layer without
#'   materialising the full data. Leave `NULL` for ordinary use.
#' @param ... Reserved for future brms-ingest options (specials, autocor
#'   handling); currently unused.
#'
#' @return An `fb_terms` object with `source = "brms"`.
#'
#' @family flexyBayes ingest adapters
#' @seealso [fb()] and [flexybayes()] for the universal fitting entry;
#'   [fb_from_asreml()] and [fb_from_greta()] for the other dialects.
#' @examples
#' df <- data.frame(
#'   y = rnorm(30),
#'   x = rnorm(30),
#'   g = factor(rep(letters[1:5], 6))
#' )
#' ir <- fb_from_brms(y ~ x + (1 | g), data = df)
#' class(ir)
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
  if (inherits(formula, "brmsformula")) {
    .check_installed(
      "brms",
      "Package 'brms' is required to ingest brmsformula objects. ",
      "Install via: install.packages(\"brms\")"
    )
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
  # are realised as a one-row placeholder data.frame (all numeric);
  # random-term level counts are stripped from the IR so .fb_preflight()
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
    data <- .fb_from_brms_placeholder_data(formula)
    fb <- fb_from_brms(
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
    # Strip per-term cached level counts -- preflight reads from
    # <fb_dataset>$dictionaries when var_n is NA. This keeps the
    # IR explicit about not knowing the true group sizes.
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
  # The hand-rolled v0.1 walker does not parse the addition-form
  # surface; refuse at ingest with a corpus-naming message rather
  # than fall through to the cryptic "Response variable not found
  # in data" error from the later names(data) lookup.
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
  if (inherits(family, "family")) {
    fam_name <- family$family
    link_name <- if (!is.null(link)) link else family$link
  } else if (is.character(family) && length(family) == 1L) {
    fl <- .resolve_family(family, link)
    fam_name <- fl$family
    link_name <- fl$link
  } else {
    stop(
      "`family` must be a character(1) name or a `family()` object.",
      call. = FALSE
    )
  }

  # ------------------------- Response --------------------------- #
  response <- deparse(formula[[2]])
  if (!response %in% names(data)) {
    stop("Response variable '", response, "' not found in data.", call. = FALSE)
  }

  # ------------------------- Walk RHS --------------------------- #
  walked <- .brms_walk_rhs(formula[[3]])

  # Reject unsupported features fail-fast (per v0.1 minimum scope)
  if (length(walked$ext_terms) > 0L) {
    # A univariate smoother gets the smoother refusal, whichever grammar
    # it was written in: `y ~ s(x)` means the same thing on both
    # surfaces, and the useful sentence is the one naming
    # `random = ~ spl(x)` on the INLA backend rather than "queued for a
    # later expansion" (field-sweep FS-18). The multivariate forms
    # (`te`, `ti`, `t2`, `gp`, `sos`) have no `spl()` counterpart, so
    # they keep the ingest-feature refusal below, which claims nothing.
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
    if (identical(rt$type, "brms_re_correlated_slope_unsupported")) {
      .stop_correlated_slope_unsupported(
        grouping_factor = rt$rhs,
        slope_variable = rt$slope_var
      )
    }
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
    if (!is.numeric(weights) || length(weights) != nrow(data)) {
      stop(
        "`weights` must be a numeric vector of length nrow(data).",
        call. = FALSE
      )
    }
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
    known_matrices = character(0)
  )

  if (!is.null(override_n_rows)) {
    data_summary$n <- as.integer(override_n_rows)
  }

  new_fb_terms(
    response = response,
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
    source = "brms"
  )
}


# Synthesise a one-row, all-numeric placeholder
# data.frame from a brms-style formula. Used by the data = NULL +
# carry_n_rows ingest path so the existing brms walker can produce
# the IR shape; the realised level counts are stripped after the
# walk so the preflight layer reads them from <fb_dataset> instead.
.fb_from_brms_placeholder_data <- function(formula) {
  all_vars <- all.vars(formula)
  if (!length(all_vars)) {
    stop(
      "formula references no variables; cannot synthesise ",
      "placeholder data.",
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
      # Bare parentheses around a non-RE expression -- descend
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

      # Unsupported brms specials -- bucket into ext_terms
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

    # Default: a fixed-effect term -- record its deparse for later
    # re-parsing through .parse_fixed
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
# v0.2.6 supported:
#   (1 | g)        -> type = "simple"             (random intercept)
#   (1 || g)       -> type = "simple"             (random intercept)
#   (x || g)       -> type = "simple_slope_uncor" (intercept + slope,
#                                                  uncorrelated; the
#                                                  intercept is implicit
#                                                  per lme4 / brms semantics)
#   (1 + x || g)   -> type = "simple_slope_uncor" (with_intercept = TRUE;
#                                                  same model as above)
#   (0 + x || g)   -> type = "simple_slope_uncor" (with_intercept = FALSE;
#                                                  slope-only)
#   (x | g)        -> type = "brms_re_correlated_slope_unsupported"
#                     (correlated slope; deferred to v0.3)
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
    # A grouping *expression* rather than a column. `(1 | g:h)` is the
    # bar spelling of the ASReml `random = ~ g:h`, and the two surfaces
    # have to produce the same term descriptor or the structural gates
    # see two different models: the ASReml spelling was refused by
    # lgm_gate() with a typed, remedy-naming refusal while the bar
    # spelling was typed as a simple intercept on a group called "g:h",
    # passed the gate, and reached the INLA binary, which failed with a
    # message about the emitted formula rather than the user's model.
    # The interaction is classified by the same walker the ASReml side
    # uses, so the descriptors are identical -- including the
    # numeric-in-a-random-interaction guard `.enrich()` fires.
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

  # Intercept + slope uncorrelated: (x || g) (bare slope -- the
  # implicit "+ 1" matches lme4 / brms semantics: see
  # ?lme4::lFormula -- the bare-slope form is sugar for
  # (1 | g) + (0 + x | g)). The intercept block is included by
  # default; explicit 0 + suppresses it (handled below).
  if (isFALSE(cor) && is_simple_var(lhs)) {
    return(.brms_re_to_fb_slope_uncor(
      slope_var = as.character(lhs),
      group_name = group_name,
      with_intercept = TRUE,
      data = data
    ))
  }

  # `lhs ± slope` forms: (1 + x || g), (0 + x || g), (x + 1 || g).
  # The `1` token signals the intercept block; the `0` token
  # suppresses it. Anything else with the same shape (e.g.
  # (x + z || g) for two slopes) falls through to the generic
  # "brms_re_unsupported" refusal -- multi-slope (x + z || g) is
  # queued for a later expansion of the brms ingest layer.
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

  # Random slope on a FACTOR: (f | g), (0 + f | g), (f || g). The bar
  # grammar has no lowering for it, but the ASReml surface of the same
  # model fits -- `random = ~ us(f):g` for the correlated form and
  # `random = ~ diag(f):g` for the uncorrelated one -- so the refusal
  # names that surface verbatim rather than listing bar spellings none of
  # which expresses the model (field-sweep FS-6 / field finding U1).
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

  # Slope correlated: (x | g) or (1 + x | g). Refused with the
  # structured deferral message.
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
      return(list(
        type = "brms_re_correlated_slope_unsupported",
        lhs = lhs_str,
        rhs = group_name,
        slope_var = slope_var,
        deparse = paste0("(", lhs_str, " | ", group_name, ")")
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

# The factor a bar-grammar random slope is written on, or NA when the
# left-hand side is not a factor slope. Recognises the bare form (`f`),
# the intercept-suppressed form (`0 + f`) and the intercept-carrying form
# (`1 + f`) in either operand order -- all three are unstructured or
# diagonal per-level effects of the same factor within the grouping
# factor, and all three have the same ASReml surface.
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
# plus the full correlation block between levels -- `us(f):g`. The
# uncorrelated form `(f || g)` drops the correlations and keeps the
# per-level variances -- `diag(f):g`. Both fit on brms; both refuse on
# INLA with the gate's own structured-covariance message, so the refusal
# says which backend carries them.
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
        "`, plus the correlations between those levels -- which is what ",
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
        "`, with no correlation between levels -- which is what the ",
        "uncorrelated\n  bar form asks for. Add the correlations with `",
        other,
        "`.\n"
      )
    },
    "backend = \"inla\" refuses both: INLA's f() machinery does not ",
    "represent a\n  structured covariance across factor levels."
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

# Raise the structured refusal for (x | g).
# The custom condition class flexybayes_correlated_slope_unsupported
# carries the grouping factor, slope variable, deferral target, and
# workaround so downstream tooling can pattern-match on the slots
# without parsing free text.
.stop_correlated_slope_unsupported <- function(
  grouping_factor,
  slope_variable
) {
  uncorrelated <- paste0("(", slope_variable, " || ", grouping_factor, ")")
  msg <- paste0(
    "Correlated random slopes (x | g) are not yet supported.\n",
    "Got: (",
    slope_variable,
    " | ",
    grouping_factor,
    ")\n\n",
    "Fit ",
    uncorrelated,
    " instead. That is the uncorrelated random\n",
    "intercept and slope: a per-",
    grouping_factor,
    " intercept deviation\n",
    "and a per-",
    grouping_factor,
    " slope on ",
    slope_variable,
    ", each with its own\n",
    "variance and no correlation parameter between them. It emits on brms\n",
    "and is the right model whenever the intercept-slope correlation is not\n",
    "itself of inferential interest.\n\n",
    "Do NOT reach for the asreml crossing ~ ",
    grouping_factor,
    " + ",
    grouping_factor,
    ":",
    slope_variable,
    ".\n",
    "With a numeric ",
    slope_variable,
    " that is one independent deviation per\n",
    grouping_factor,
    "-by-",
    slope_variable,
    " cell, not a per-",
    grouping_factor,
    " slope -- a\n",
    "different model with a different parameter count.\n\n",
    "If the correlation parameter itself is needed, it is deferred to a\n",
    "future release (structured-covariance representation). No active\n",
    "backend represents the correlated form today."
  )
  cond <- structure(
    class = c("flexybayes_correlated_slope_unsupported", "error", "condition"),
    list(
      message = msg,
      call = NULL,
      grouping_factor = grouping_factor,
      slope_variable = slope_variable,
      deferral_target = "a future release",
      workaround = uncorrelated
    )
  )
  stop(cond)
}
