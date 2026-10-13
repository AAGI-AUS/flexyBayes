# =============================================================================
# One family gate for every entry.
#
# A family object passed to fb_from_brms(), such as brms::lognormal() or
# brms::student(), goes through .resolve_family(), the same gate as the
# string spelling and the ASReml grammar. It reaches the same engines as
# its string, or refuses with a typed condition.
#
# The families are read from the engines: every family brms carries a
# `.family_<name>()` definition for and brmsfamily() constructs, and the
# `stats` family constructors, with a non-default link, a mixture, a custom
# family and two unevaluated generators beside them. Each object is asked for
# its code on both engines, for its route under `auto`, and by the streamed
# entry, and nothing is sampled. The oracle is the string spelling: an object
# reaches exactly the outcome its name and link reach, and every refusal
# carries a registered reason code.
# =============================================================================

.fog_families <- function() {
  ns <- ls(asNamespace("brms"), all.names = TRUE)
  nm_vec <- sort(sub("^\\.family_", "", grep("^\\.family_", ns, value = TRUE)))
  out <- list()
  for (nm in nm_vec) {
    obj <- tryCatch(brms::brmsfamily(nm), error = function(e) NULL)
    if (!is.null(obj)) {
      out[[paste0("brms::brmsfamily(\"", nm, "\")")]] <- obj
    }
  }
  stats_vec <- c(
    "binomial", "Gamma", "gaussian", "inverse.gaussian", "poisson",
    "quasi", "quasibinomial", "quasipoisson"
  )
  for (nm in stats_vec) {
    out[[paste0("stats::", nm, "()")]] <- get(nm, envir = asNamespace("stats"))()
  }
  out[["stats::binomial(link = \"probit\")"]] <- stats::binomial(link = "probit")
  out[["brms::mixture(gaussian, gaussian)"]] <- suppressMessages(
    brms::mixture(stats::gaussian, stats::gaussian)
  )
  out[["brms::custom_family(\"myfam\")"]] <- brms::custom_family(
    "myfam",
    dpars = "mu",
    links = "identity"
  )
  out[["stats::poisson, unevaluated"]] <- stats::poisson
  out[["brms::student, unevaluated"]] <- brms::student
  out
}

# A 48-row design with a response on the support of the named family.
.fog_data <- function(fam_name) {
  set.seed(1L)
  d <- data.frame(g = factor(rep(seq_len(6L), each = 8L)))
  eta <- 0.3 + stats::rnorm(6L, 0, 0.3)[as.integer(d$g)]
  n <- nrow(d)
  d$y <- if (grepl("bernoulli|binomial", fam_name)) {
    stats::rbinom(n, 1L, stats::plogis(eta))
  } else if (grepl("poisson|negbinomial|geometric|discrete_weibull", fam_name)) {
    stats::rpois(n, exp(eta))
  } else if (grepl("beta", fam_name)) {
    stats::plogis(eta + stats::rnorm(n, 0, 0.3))
  } else if (grepl("hurdle|zero_inflated", fam_name)) {
    ifelse(
      stats::rbinom(n, 1L, 0.7) == 1L,
      exp(eta + stats::rnorm(n, 0, 0.3)),
      0
    )
  } else {
    exp(eta + stats::rnorm(n, 0, 0.3))
  }
  d
}

# The outcome class of one call: admitted, the route a plan chose, the
# registered reason code of a refusal, or the message of an untyped error.
.fog_outcome <- function(expr) {
  out <- tryCatch(suppressMessages(expr), error = function(e) e)
  if (inherits(out, "flexybayes_refusal")) {
    return(paste0("refusal:", out$reason_code))
  }
  if (inherits(out, "error")) {
    return(paste0("untyped:", conditionMessage(out)))
  }
  if (inherits(out, "fb_plan")) {
    return(paste0("plan:", out$backend_chosen))
  }
  "admitted"
}

.fog_calls <- function(family, link, d) {
  list(
    brms = .fog_outcome(flexybayes(
      y ~ 1 + (1 | g), data = d, family = family, link = link,
      backend = "brms", return_code = TRUE, verbose = FALSE
    )),
    inla = .fog_outcome(flexybayes(
      y ~ 1 + (1 | g), data = d, family = family, link = link,
      backend = "inla", return_code = TRUE, verbose = FALSE
    )),
    auto = .fog_outcome(flexybayes(
      y ~ 1 + (1 | g), data = d, family = family, link = link,
      backend = "auto", plan = TRUE, verbose = FALSE
    )),
    stream = .fog_outcome(flexybayes_stream(
      fixed = y ~ 1, random = ~g, source = d, family = family,
      fit = FALSE, verbose = FALSE
    ))
  )
}

test_that("every family object gets the outcome its name gets, on every entry", {
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE
  )
  fam_list <- .fog_families()
  # brms 2.23.0 constructs 49 families. A roster read that found none would
  # leave the loop below with nothing to check.
  expect_gte(length(fam_list), 50L)

  admit_brms <- c(
    "gaussian", "binomial", "poisson", "gamma", "Gamma", "beta",
    "hurdle_gamma", "lognormal"
  )
  admit_inla <- c("gaussian", "binomial", "poisson", "gamma", "Gamma", "beta")
  # INLA takes each family's default link only (test-family-link-gate.R), so
  # stats::Gamma(), whose link is inverse, and binomial(link = "probit") are
  # refused there.
  inla_links <- c(
    gaussian = "identity", binomial = "logit", poisson = "log", gamma = "log",
    beta = "logit"
  )
  # flexybayes_stream() has no `link` argument: it aggregates these three
  # families on the links their names default to.
  stream_links <- c(gaussian = "identity", binomial = "logit", poisson = "log")

  checked <- 0L
  for (i1 in seq_along(fam_list)) {
    key <- names(fam_list)[[i1]]
    fam <- fam_list[[i1]]
    obj <- if (is.function(fam)) fam() else fam
    d <- .fog_data(tolower(obj$family))
    by_object <- .fog_calls(fam, NULL, d)
    by_name <- .fog_calls(obj$family, obj$link, d)

    for (entry in c("brms", "inla", "auto")) {
      expect_identical(
        by_object[[entry]],
        by_name[[entry]],
        label = paste(key, "on", entry, "against its name")
      )
    }
    for (entry in names(by_object)) {
      expect_false(
        startsWith(by_object[[entry]], "untyped:"),
        label = paste(key, "on", entry, "is typed:", by_object[[entry]])
      )
    }
    expect_identical(
      identical(by_object$brms, "admitted"),
      obj$family %in% admit_brms,
      label = paste(key, "admitted on brms")
    )
    expect_identical(
      identical(by_object$inla, "admitted"),
      obj$family %in% admit_inla &&
        identical(unname(inla_links[tolower(obj$family)]), obj$link),
      label = paste(key, "admitted on INLA")
    )
    expect_identical(
      identical(by_object$stream, "admitted"),
      obj$family %in% names(stream_links) &&
        identical(unname(stream_links[obj$family]), obj$link),
      label = paste(key, "admitted by flexybayes_stream()")
    )
    checked <- checked + 1L
  }
  expect_identical(checked, length(fam_list))

  # lognormal is admitted for brms only: INLA refuses it by name on the
  # family row, and `auto` routes it to brms.
  ln <- .fog_calls(brms::lognormal(), NULL, .fog_data("lognormal"))
  expect_identical(ln$brms, "admitted")
  expect_identical(ln$inla, "refusal:inla_gate_refused")
  expect_identical(ln$auto, "plan:brms")
  st <- .fog_calls(brms::student(), NULL, .fog_data("student"))
  expect_identical(st$brms, "refusal:unsupported_family")
})
