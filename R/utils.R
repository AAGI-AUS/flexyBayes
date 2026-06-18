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

# Resolve family and link function
#
# @param family Character: gaussian, binomial, poisson, negative_binomial,
#   gamma, beta
# @param link Character or NULL: override default link
# @return List with family and link
.resolve_family <- function(family, link) {
  defaults <- list(
    gaussian = "identity",
    binomial = "logit",
    binary = "logit",
    poisson = "log",
    negative_binomial = "log",
    negbinom = "log",
    gamma = "log",
    beta = "logit"
  )
  fam <- tolower(family)
  # The generalised extreme value and Dirichlet families are not GLM-link
  # mixed models -- block maxima have no mean-link, and a Dirichlet response
  # is a simplex, not a scalar -- so they are fitted through their own
  # dedicated entry points rather than the formula emit path. Point the user
  # there explicitly instead of refusing them as merely "unsupported".
  if (fam %in% c("gen_extreme_value", "gev", "extreme_value")) {
    stop(
      "Family \"", family, "\" (generalised extreme value) is fitted via ",
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
        "see `?fb_gev` and `?fb_dirichlet`."
      ),
      family_class = "flexybayes_unsupported_family"
    ))
  }
  lnk <- if (!is.null(link)) tolower(link) else defaults[[fam]]
  list(family = fam, link = lnk)
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
    "gamma" = Gamma(link = fam_link$link),
    "beta" = gaussian(link = "identity"),
    gaussian()
  )
}
