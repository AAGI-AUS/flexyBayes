# Internal utility helpers for flexyBayes
# Not exported

# Add a line of code to the code generation context
.add <- function(ctx, ...) {
  ctx$code <- c(ctx$code, paste0(...))
  ctx
}

# Register a parameter name in the context
.add_param <- function(ctx, ...) {
  ctx$params <- c(ctx$params, c(...))
  ctx
}

# Add a term to the linear predictor
.add_pred <- function(ctx, term_str) {
  ctx$predictor <- c(ctx$predictor, term_str)
  ctx
}

# Safe deparse of formula sub-expression
.dep <- function(expr, i) {
  if (length(expr) >= i) deparse(expr[[i]]) else NA_character_
}

# Safe deparse of a named argument inside a call expression. Returns
# NA_character_ when the argument is absent. Used by the formula parser
# to extract structured-covariance named arguments (chol, precision,
# ...) from vm() / ped() special-term calls without disturbing the
# existing
# positional-V backward-compat path.
.dep_named <- function(expr, name) {
  nms <- names(expr)
  if (is.null(nms)) {
    return(NA_character_)
  }
  idx <- which(nms == name)
  if (length(idx) == 0L) {
    return(NA_character_)
  }
  deparse(expr[[idx[[1L]]]])
}

# Normalise the `family` argument to the character spelling the rest of
# the package works in.
#
# `flexybayes()` documents `family` as a string, but `stats::binomial()`
# and friends are what an R user reaches for first, and the roxygen
# imports name them. Until 0.9.0 a family object reached `tolower()` and
# the scalar `if` below it, and the user saw the raw base-R error "the
# condition has length > 1" -- untyped, and about the wrong thing.
#
# A family object carries both halves of the specification, so its
# `$link` is honoured. That is faithful rather than convenient:
# `Gamma()` names the inverse link, not the log link `family = "gamma"`
# defaults to, and reading it as anything else would fit a model the
# call did not ask for. Supplying `link =` as well is refused rather
# than resolved, because either resolution silently discards half of
# what the call said.
.normalise_family_argument <- function(family, link) {
  # A generator passed unevaluated (`family = binomial`) is the other
  # idiom base R accepts; evaluate it and carry on with the object.
  if (is.function(family)) {
    family <- tryCatch(family(), error = function(e) e)
    if (inherits(family, "error")) {
      stop(.fb_refusal_condition(
        reason_code = "family_argument_not_recognised",
        message = paste0(
          "`family` was a function that did not evaluate to a family ",
          "object: ",
          conditionMessage(family),
          ". Pass a family name as a string (for example ",
          "family = \"binomial\"), or a `stats` family object such as ",
          "binomial()."
        ),
        family_class = "flexybayes_family_argument_refusal"
      ))
    }
  }

  if (inherits(family, "family")) {
    fam <- family$family
    lnk <- family$link
    if (
      !is.character(fam) ||
        length(fam) != 1L ||
        !is.character(lnk) ||
        length(lnk) != 1L
    ) {
      stop(.fb_refusal_condition(
        reason_code = "family_argument_not_recognised",
        message = paste0(
          "`family` carries the class \"family\" but not a usable ",
          "$family / $link pair, so the model it names cannot be read. ",
          "Pass a family name as a string (for example ",
          "family = \"binomial\"), or a well-formed `stats` family ",
          "object such as binomial()."
        ),
        family_class = "flexybayes_family_argument_refusal"
      ))
    }
    if (!is.null(link) && !identical(tolower(link), tolower(lnk))) {
      stop(.fb_refusal_condition(
        reason_code = "family_argument_not_recognised",
        message = paste0(
          "`family = ",
          fam,
          "(link = \"",
          lnk,
          "\")` and `link = \"",
          link,
          "\"` name different links. flexyBayes refuses rather than ",
          "choosing one: drop `link =`, or pass the family as a string ",
          "(family = \"",
          tolower(fam),
          "\", link = \"",
          link,
          "\")."
        ),
        family_class = "flexybayes_family_argument_refusal"
      ))
    }
    return(list(family = fam, link = lnk))
  }

  if (!is.character(family) || length(family) != 1L || !nzchar(family)) {
    stop(.fb_refusal_condition(
      reason_code = "family_argument_not_recognised",
      message = paste0(
        "`family` must be a single family name or a `stats` family ",
        "object; got an object of class ",
        paste(class(family), collapse = "/"),
        " of length ",
        length(family),
        ". Pass, for example, family = \"binomial\" or ",
        "family = binomial()."
      ),
      family_class = "flexybayes_family_argument_refusal"
    ))
  }

  list(family = family, link = link)
}

# Resolve family and link function
#
# @param family Character family name, or a `stats` family object such
#   as `binomial()`; a family object supplies its own link.
# @param link Character or NULL: override default link
# @return List with family and link
.resolve_family <- function(family, link) {
  normalised <- .normalise_family_argument(family, link)
  family <- normalised$family
  link <- normalised$link

  defaults <- list(
    gaussian = "identity",
    binomial = "logit",
    binary = "logit",
    poisson = "log",
    negative_binomial = "log",
    negbinom = "log",
    gamma = "log",
    beta = "logit",
    # brms-native, brms-only: `brms::brmsfamily("hurdle_gamma")` declares
    # dpars mu, shape, hu on brms 2.23.0, and the entry allowlist used to
    # refuse a family the engine behind it carries (field-sweep FS-4 /
    # field finding C2). INLA's likelihood roster has no counterpart, so
    # lgm_gate() refuses `backend = "inla"` for it on the family row.
    hurdle_gamma = "log"
  )
  fam <- tolower(family)
  # The generalised extreme value and Dirichlet families are not GLM-link
  # mixed models -- block maxima have no mean-link, and a Dirichlet response
  # is a simplex, not a scalar -- so they are fitted through their own
  # dedicated entry points rather than the formula emit path. Point the user
  # there explicitly instead of refusing them as merely "unsupported".
  if (fam %in% c("gen_extreme_value", "gev", "extreme_value")) {
    stop(
      "Family \"",
      family,
      "\" (generalised extreme value) is fitted via ",
      "the dedicated `fb_gev()` entry point, not the `flexybayes()` ",
      "formula path: block maxima have no mean-link parameterisation. See ",
      "`?fb_gev`.",
      call. = FALSE
    )
  }
  if (fam == "dirichlet") {
    stop(
      "Family \"dirichlet\" is fitted via the dedicated `fb_dirichlet()` ",
      "entry point, not the `flexybayes()` formula path: a Dirichlet ",
      "response is a simplex (proportions summing to one), not a scalar ",
      "with a mean-link. See `?fb_dirichlet`.",
      call. = FALSE
    )
  }
  if (!fam %in% names(defaults)) {
    stop(.fb_refusal_condition(
      reason_code = "unsupported_family",
      message = paste0(
        "Unsupported family \"",
        family,
        "\". flexyBayes supports: ",
        paste(names(defaults), collapse = ", "),
        ". Other families (including survival / time-to-event and ",
        "multivariate responses) are planned future work; see ",
        "fb_refusals(). Block-maxima (generalised extreme value) and ",
        "compositional (Dirichlet) data have dedicated fitters: ",
        "see `?fb_gev` and `?fb_dirichlet`.",
        .fb_family_boundary_note(fam)
      ),
      family_class = "flexybayes_unsupported_family"
    ))
  }
  lnk <- if (!is.null(link)) tolower(link) else defaults[[fam]]
  list(family = fam, link = lnk)
}

# The documented-boundary note for a family a user is likely to reach
# for and not find. Appended to the unsupported-family refusal so the
# message says which layer the boundary is at, rather than leaving the
# user to infer that flexyBayes is merely behind its engines.
#
# The three families here were the ones the field engagement named
# (finding C2). Their status was checked against the installed engines
# rather than recalled -- `brms::brmsfamily()` refuses all three on brms
# 2.23.0, so blocking them is not a narrowing over that engine; INLA's
# roster (`names(INLA::inla.models()$likelihood)`, 107 entries at INLA
# 25.10.19) does carry `tweedie`, so that one is a genuine flexyBayes
# boundary and says so. Getting this backwards changes what the fix is:
# the fourth family the field named, `hurdle_gamma`, IS brms-native and
# is now admitted above.
.fb_family_boundary_note <- function(fam) {
  notes <- list(
    tweedie = paste0(
      "\n\nBoundary note. INLA's likelihood roster does carry `tweedie`, ",
      "so this is a\nflexyBayes boundary rather than an engine one: the ",
      "package has no INLA emit for\nit (no link / power parameter, no ",
      "validated fit). brms carries no Tweedie family\nat all. Tracked ",
      "as a feature request in inst/KNOWN_ISSUES.md; for compound-Poisson\n",
      "gamma data today, fit the zero / positive parts separately and ",
      "recombine on the\nposterior, or use family = \"hurdle_gamma\" with ",
      "backend = \"brms\"."
    ),
    zero_inflated_gamma = paste0(
      "\n\nBoundary note. Neither active engine carries this family: ",
      "`brms::brmsfamily(\n\"zero_inflated_gamma\")` refuses on brms ",
      "2.23.0, and INLA's roster has no\ncounterpart. The nearest ",
      "implemented alternative is family = \"hurdle_gamma\" with\n",
      "backend = \"brms\" -- a hurdle and a zero-inflated gamma are the ",
      "same model when the\npositive part has no mass at zero, which a ",
      "gamma does not."
    ),
    compound_poisson = paste0(
      "\n\nBoundary note. `brms::brmsfamily(\"compound_poisson\")` ",
      "refuses on brms 2.23.0, so\nthis is not a narrowing over that ",
      "engine. INLA's roster carries `tweedie`, which\nis the ",
      "compound-Poisson gamma; flexyBayes has no emit for it. See the ",
      "`tweedie`\nentry in inst/KNOWN_ISSUES.md."
    )
  )
  notes[[fam]] %||% ""
}

# Map flexyBayes family string to stats::family object
.get_stats_family <- function(fam_link) {
  switch(
    fam_link$family,
    "gaussian" = gaussian(link = fam_link$link),
    "binomial" = ,
    "binary" = binomial(link = fam_link$link),
    "poisson" = poisson(link = fam_link$link),
    "negative_binomial" = ,
    "negbinom" = poisson(link = "log"),
    "gamma" = ,
    # The starting-value GLM for a hurdle gamma is the gamma GLM of its
    # positive part; the zero-mass component has no `stats` counterpart
    # and does not inform the mean-model starting values.
    "hurdle_gamma" = Gamma(link = fam_link$link),
    "beta" = gaussian(link = "identity"),
    gaussian()
  )
}
