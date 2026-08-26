# fb -- flexyBayes user-facing entry points and ingest dispatchers.
#
# - fb_from_asreml(): asreml-format ingest path,
#   wraps the existing parse_formula.R helpers (.parse_fixed,
#   .parse_formula, .resolve_family) and emits an fb_terms object
#   (the intermediate representation, IR; see R/fb_terms.R).
# - fb(): literal alias for flexybayes(). The
#   brms-format ingest (former fb() body) is deferred to v0.2 under
#   the working name fb_brms().

# ---------------------------------------------------------------- #
# Asreml-format ingest -> fb_terms                                 #
# ---------------------------------------------------------------- #

#' Ingest an ASReml-format model specification into the flexyBayes IR
#'
#' Parses an ASReml-style `fixed` / `random` / `residual` specification into
#' a `fb_terms` object -- flexyBayes's backend-agnostic intermediate
#' representation (IR) of a model. The IR is what every engine emits
#' from, so building it explicitly lets a power user inspect the parsed
#' model, cache it, or hand it to a fitting verb. Most users never call
#' this directly: [flexybayes()] and [fb()] build the IR internally from
#' the same arguments. Argument names and defaults mirror [flexybayes()]
#' one-for-one, so this is a drop-in for that function's parsing step.
#'
#' @param fixed Two-sided formula: `response ~ fixed_effects`.
#' @param random One-sided formula `~ random_terms` (ASReml syntax), or
#'   `NULL`.
#' @param residual One-sided formula `~ residual_structure`, or `NULL`. `NULL`
#'   defaults to iid residuals (`list(list(type = "units"))`), matching
#'   [flexybayes()].
#' @param rcov Defunct in flexyBayes 0.9.0. The ASReml 3 name for the
#'   residual-structure argument, renamed to `residual` in ASReml 4;
#'   supplying it now raises an error. Use `residual` instead.
#' @param data A data.frame containing every variable referenced.
#' @param family Character family name (`gaussian`, `binomial`, `binary`,
#'   `poisson`, `negative_binomial`, `negbinom`, `gamma`, `beta`).
#' @param link Character link override, or `NULL` for the family default.
#' @param weights Optional numeric vector of length `nrow(data)`. When
#'   non-`NULL`, mapped to a single addition term with `type = "weights"`.
#' @param known_matrices Named list of known matrices (e.g.
#'   `list(Gmat = G)`); names are recorded in
#'   `data_summary$known_matrices`.
#' @param prior An [fb_prior()] object, or `NULL`. `NULL` falls back to
#'   the scalar priors `prior_fixed_sd` / `prior_vc_sd`.
#' @param prior_fixed_sd Numeric SD for fixed-effect normal priors when
#'   `prior` is `NULL`.
#' @param prior_vc_sd Numeric hyperparameter for the variance-component
#'   priors when `prior` is `NULL`.
#'
#' @return An `fb_terms` object with `source = "asreml"`.
#'
#' @family flexyBayes ingest adapters
#' @seealso [flexybayes()] and [fb()] for the universal fitting entry;
#'   [fb_from_brms()] for the other ingest dialect.
#' @examples
#' df <- data.frame(
#'   yield = rnorm(20),
#'   geno  = factor(rep(letters[1:4], 5)),
#'   env   = factor(rep(c("a", "b"), 10))
#' )
#' ir <- fb_from_asreml(yield ~ env, random = ~ geno, data = df)
#' class(ir)
#' @export
fb_from_asreml <- function(
  fixed,
  random = NULL,
  residual = NULL,
  data,
  family = "gaussian",
  link = NULL,
  weights = NULL,
  known_matrices = list(),
  prior = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  rcov = lifecycle::deprecated()
) {
  # `rcov` is the ASReml 3 spelling of the residual-structure argument,
  # renamed to `residual` in ASReml 4 and in flexyBayes 0.9.0. Kept only
  # as a tripwire that raises a guiding error if supplied.
  if (lifecycle::is_present(rcov)) {
    lifecycle::deprecate_stop(
      when = "0.9.0",
      what = "fb_from_asreml(rcov)",
      with = "fb_from_asreml(residual)",
      details = paste0(
        "ASReml-R renamed this argument from `rcov` (ASReml 3) to ",
        "`residual` (ASReml 4)."
      )
    )
  }

  # Delegate parsing to existing helpers verbatim -- no behaviour
  # change relative to flexybayes().
  fixed_info <- .parse_fixed(fixed, data)
  random_terms <- if (!is.null(random)) {
    .parse_formula(random, data)
  } else {
    list()
  }
  residual_terms <- if (!is.null(residual)) {
    .parse_formula(residual, data)
  } else {
    list(list(type = "units"))
  }
  fam_link <- .resolve_family(family, link)

  # weights -> addition_terms (per IR design sketch Sec.2)
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

  # Priors: forward-compat -- accept an fb_prior object when supplied;
  # otherwise capture legacy scalar priors. Marked legacy = TRUE so
  # downstream code can detect and dispatch accordingly.
  priors <- if (!is.null(prior)) {
    prior
  } else {
    list(fixed_sd = prior_fixed_sd, vc_sd = prior_vc_sd, legacy = TRUE)
  }

  data_summary <- list(
    n = nrow(data),
    known_matrices = names(known_matrices)
  )

  new_fb_terms(
    response = fixed_info$response,
    family = fam_link$family,
    link = fam_link$link,
    intercept = fixed_info$intercept,
    fixed_terms = fixed_info$terms,
    random_terms = random_terms,
    residual_terms = residual_terms,
    addition_terms = addition_terms,
    priors = priors,
    data_summary = data_summary,
    capabilities = character(),
    source = "asreml"
  )
}

# ---------------------------------------------------------------- #
# Grammar polymorphism for the universal entry                    #
# ---------------------------------------------------------------- #

# .detect_grammar() -- deterministic call-shape dispatch for the
# universal entry (fb() / flexybayes()). The rule is by SHAPE, never by
# guessing intent:
#   - a forced `syntax` (anything but "auto") wins outright;
#   - a brmsformula object  -> "brms";
#   - a formula whose terms carry lme4 / brms bar-grouping ("(... | g)")
#     -> "brms" (ASReml `fixed` formulae never contain a bar);
#   - any other formula      -> "asreml" (the long-standing default; a
#     bar-free `y ~ x` is identical under either grammar, so the safe
#     default matches historical behaviour).
# A non-formula `spec` falls through to "asreml" so the
# established fb_from_asreml() validation owns the error message
# (behaviour-preserving for malformed input).
.detect_grammar <- function(
  spec,
  random = NULL,
  residual = NULL,
  syntax = "auto"
) {
  if (!identical(syntax, "auto")) {
    return(syntax)
  }
  if (inherits(spec, "brmsformula")) {
    return("brms")
  }
  if (inherits(spec, "formula")) {
    has_bar <- any(grepl("|", as.character(spec), fixed = TRUE))
    if (has_bar) return("brms")
  }
  "asreml"
}

# .build_ir_polymorphic() -- the universal entry's ingest step. Detects
# the grammar from the call shape and routes to the matching exported
# adapter, producing the backend-agnostic fb_terms IR. ASReml ingest is
# byte-identical to the historical direct fb_from_asreml() call;
# brms-grammar ingest was added later.
.build_ir_polymorphic <- function(
  fixed,
  random,
  residual,
  data,
  family,
  link,
  weights,
  known_matrices,
  prior,
  prior_fixed_sd,
  prior_vc_sd,
  syntax = "auto"
) {
  grammar <- .detect_grammar(fixed, random, residual, syntax)
  switch(
    grammar,
    asreml = fb_from_asreml(
      fixed = fixed,
      random = random,
      residual = residual,
      data = data,
      family = family,
      link = link,
      weights = weights,
      known_matrices = known_matrices,
      prior = prior,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd
    ),
    brms = {
      if (!is.null(random) || !is.null(residual)) {
        stop(.fb_refusal_condition(
          reason_code = "grammar_brms_with_asreml_terms",
          message = paste0(
            "A brms-style formula (with `(... | g)` grouping) cannot be ",
            "combined with `random` / `residual` (ASReml grammar). Put every ",
            "grouping term inside the formula, or use the ASReml ",
            "`fixed` / `random` / `residual` form throughout."
          ),
          family_class = "flexybayes_grammar_brms_with_asreml_terms"
        ))
      }
      if (length(known_matrices) > 0L) {
        stop(.fb_refusal_condition(
          reason_code = "grammar_brms_known_matrices_unsupported",
          message = paste0(
            "`known_matrices` is not supported with brms-grammar ingest ",
            "via fb() / flexybayes(). Use the ASReml form ",
            "(flexybayes(fixed = , random = , known_matrices = ) or ",
            "fb_from_asreml()) for known-matrix carriers."
          ),
          family_class = "flexybayes_grammar_brms_known_matrices_unsupported"
        ))
      }
      fb_from_brms(
        formula = fixed,
        data = data,
        family = family,
        link = link,
        prior = prior,
        weights = weights,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd
      )
    }
  )
}

# Note: the `fb` alias for `flexybayes()` is defined at the bottom
# of R/flexybayes.R so that the alias target is in scope at source
# time (R sources files alphabetically; this file loads before
# flexybayes.R).
