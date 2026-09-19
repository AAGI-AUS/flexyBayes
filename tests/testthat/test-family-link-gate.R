# =============================================================================
# Every family and link pair is emitted as that link or refused by name.
#
# The family gate checked the family name and never the link. The brms
# emitter built Beta() and lognormal() without the link, so
# brms::Beta(link = "probit") authored inv_logit code. The INLA emit passes
# no link to INLA at all, so every non-default link was fitted under INLA's
# default while the fit recorded the link the call asked for: stats::Gamma()
# was fitted with a log link, binomial(link = "probit") with a logit link.
# Links a family constructor does not know ended in an untyped error.
#
# The pairs are read from the engines: for each family the gate admits,
# every link brms's own family table lists for the brms family the emitter
# uses, every link INLA lists for its likelihood, the default, and one
# spelling no engine carries. Each pair is asked for its code on both
# engines and for its route under `auto`, and nothing is sampled. On brms the
# oracle is brms's own Stan code for that family and link. On INLA it is the
# link INLA reports fitting when the call names none, read from a 48-row fit
# per likelihood. Every pair matches its oracle or carries a registered
# reason code.
# =============================================================================

# The brms family and the INLA likelihood each flexyBayes family is emitted as.
.flg_brms_name <- c(
  gaussian = "gaussian",
  binomial = "bernoulli",
  binary = "bernoulli",
  poisson = "poisson",
  negative_binomial = "negbinomial",
  negbinom = "negbinomial",
  gamma = "gamma",
  beta = "beta",
  hurdle_gamma = "hurdle_gamma",
  lognormal = "lognormal"
)
.flg_inla_name <- c(
  gaussian = "gaussian",
  binomial = "binomial",
  binary = "binomial",
  poisson = "poisson",
  negative_binomial = "nbinomial",
  negbinom = "nbinomial",
  gamma = "gamma",
  beta = "beta"
)

# The links flexyBayes fits on brms, default first: the pairs whose Stan code
# matched brms's own at 0.12.1 plus U-FAM-GATE-1. Freeze by list, so no link
# is added here without an emitter change that makes it correct.
.flg_brms_links <- list(
  gaussian = c(
    "identity",
    "log",
    "inverse",
    "logit",
    "probit",
    "cloglog",
    "cauchit"
  ),
  binomial = c(
    "logit",
    "probit",
    "probit_approx",
    "cloglog",
    "cauchit",
    "softit",
    "identity",
    "log"
  ),
  binary = c(
    "logit",
    "probit",
    "probit_approx",
    "cloglog",
    "cauchit",
    "softit",
    "identity",
    "log"
  ),
  poisson = c("log", "identity", "sqrt"),
  negative_binomial = c("log", "identity", "sqrt", "softplus", "squareplus"),
  negbinom = c("log", "identity", "sqrt", "softplus", "squareplus"),
  gamma = c("log", "identity", "inverse"),
  beta = "logit",
  hurdle_gamma = c("log", "identity", "inverse", "softplus", "squareplus"),
  lognormal = "identity"
)

.flg_data <- function(fam) {
  set.seed(1L)
  d <- data.frame(g = factor(rep(seq_len(6L), each = 8L)))
  eta <- 0.3 + stats::rnorm(6L, 0, 0.3)[as.integer(d$g)]
  n <- nrow(d)
  d$y <- switch(
    fam,
    binomial = ,
    binary = stats::rbinom(n, 1L, stats::plogis(eta)),
    poisson = ,
    negative_binomial = ,
    negbinom = stats::rpois(n, exp(eta)),
    beta = stats::plogis(eta + stats::rnorm(n, 0, 0.3)),
    hurdle_gamma = ifelse(
      stats::rbinom(n, 1L, 0.7) == 1L,
      exp(eta + stats::rnorm(n, 0, 0.3)),
      0
    ),
    exp(eta + stats::rnorm(n, 0, 0.3))
  )
  d
}

# The likelihood section of a Stan model block: where brms writes the inverse
# link (`mu = Phi(mu);`) or a link-specific density (`bernoulli_logit_lpmf`).
.flg_lik_section <- function(code) {
  lin <- strsplit(
    paste(as.character(code), collapse = "\n"),
    "\n",
    fixed = TRUE
  )[[1L]]
  a <- grep("// likelihood including constants", lin, fixed = TRUE)
  b <- grep("// priors including constants", lin, fixed = TRUE)
  if (length(a) != 1L || length(b) != 1L) {
    return(NA_character_)
  }
  paste(trimws(lin[a:b]), collapse = "\n")
}

# The outcome of one call: the registered reason code of a refusal, the
# message of any other error, or the returned value.
.flg_call <- function(expr) {
  out <- tryCatch(suppressMessages(expr), error = function(e) e)
  if (inherits(out, "flexybayes_refusal")) {
    return(list(class = paste0("refusal:", out$reason_code)))
  }
  if (inherits(out, "error")) {
    return(list(class = paste0("untyped:", conditionMessage(out))))
  }
  list(class = "admitted", value = out)
}

test_that("every family and link pair is emitted as that link or refused by name", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE
  )
  brms_ns <- asNamespace("brms")
  brms_list <- function(fam) {
    get(paste0(".family_", .flg_brms_name[[fam]]), envir = brms_ns)()$links
  }

  # --- the engines' answers ------------------------------------------------ #
  inla_default <- list()
  for (nm in unique(.flg_inla_name)) {
    fam <- names(.flg_inla_name)[match(nm, .flg_inla_name)]
    fit <- INLA::inla(
      y ~ 1 + f(g, model = "iid"),
      family = nm,
      data = .flg_data(fam)
    )
    inla_default[[nm]] <- fit$misc$linkfunctions$names[[1L]]
  }
  ref_cache <- list()
  brms_ref <- function(fam, link) {
    key <- paste(.flg_brms_name[[fam]], link)
    if (is.null(ref_cache[[key]])) {
      ref_cache[[key]] <<- .flg_lik_section(brms::make_stancode(
        y ~ 1 + (1 | g),
        data = .flg_data(fam),
        family = brms::brmsfamily(.flg_brms_name[[fam]], link = link)
      ))
    }
    ref_cache[[key]]
  }

  # The frozen list against the engines: every admitted link is one brms
  # carries, its code is distinguishable from the default link's, and the
  # default is the link INLA fits when the call names none.
  for (fam in names(.flg_brms_links)) {
    lk_vec <- .flg_brms_links[[fam]]
    expect_true(
      all(lk_vec %in% brms_list(fam)),
      label = paste(fam, "links are brms's")
    )
    for (lk in lk_vec[-1L]) {
      expect_false(
        identical(brms_ref(fam, lk), brms_ref(fam, lk_vec[[1L]])),
        label = paste(fam, lk, "authors code distinct from the default link")
      )
    }
    if (fam %in% names(.flg_inla_name)) {
      expect_identical(
        inla_default[[.flg_inla_name[[fam]]]],
        lk_vec[[1L]],
        label = paste(fam, "default link is INLA's")
      )
    }
  }

  # --- the pairs ----------------------------------------------------------- #
  pair_list <- list()
  for (fam in names(.flg_brms_name)) {
    inla_vec <- if (fam %in% names(.flg_inla_name)) {
      setdiff(
        INLA::inla.models()$likelihood[[.flg_inla_name[[fam]]]]$link,
        "default"
      )
    } else {
      character(0)
    }
    for (lk in c(
      NA_character_,
      unique(c(brms_list(fam), inla_vec)),
      "not_a_link"
    )) {
      pair_list[[length(pair_list) + 1L]] <- list(
        key = paste0("\"", fam, "\", link ", lk),
        family = fam,
        link = if (is.na(lk)) NULL else lk,
        name = fam,
        requested = if (is.na(lk)) .flg_brms_links[[fam]][[1L]] else lk
      )
    }
  }
  obj_list <- list(
    "stats::Gamma()" = stats::Gamma(),
    "stats::gaussian(link = \"log\")" = stats::gaussian(link = "log"),
    "stats::poisson(link = \"sqrt\")" = stats::poisson(link = "sqrt"),
    "stats::binomial(link = \"probit\")" = stats::binomial(link = "probit"),
    "brms::Beta(link = \"probit\")" = brms::Beta(link = "probit"),
    "brms::lognormal(link = \"inverse\")" = brms::lognormal(link = "inverse")
  )
  for (nm in names(obj_list)) {
    pair_list[[length(pair_list) + 1L]] <- list(
      key = nm,
      family = obj_list[[nm]],
      link = NULL,
      name = tolower(obj_list[[nm]]$family),
      requested = obj_list[[nm]]$link
    )
  }
  # 121 name and link pairs from brms 2.23.0 and INLA 26.8.7, and six objects.
  # A roster read that found nothing would leave the loop with nothing to check.
  expect_gte(length(pair_list), 100L)

  checked <- 0L
  for (p in pair_list) {
    d <- .flg_data(p$name)
    on_brms <- p$requested %in% .flg_brms_links[[p$name]]
    on_inla <- p$name %in%
      names(.flg_inla_name) &&
      identical(p$requested, inla_default[[.flg_inla_name[[p$name]]]])

    brms_out <- .flg_call(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = p$family,
      link = p$link,
      backend = "brms",
      return_code = TRUE,
      verbose = FALSE
    ))
    # Code authored for a link brms does not carry cannot carry that link, and
    # brms has no code of its own to compare it with.
    brms_class <- if (identical(brms_out$class, "admitted")) {
      if (
        p$requested %in%
          brms_list(p$name) &&
          identical(
            .flg_lik_section(brms_out$value),
            brms_ref(p$name, p$requested)
          )
      ) {
        "emitted"
      } else {
        "changed"
      }
    } else {
      brms_out$class
    }

    inla_out <- .flg_call(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = p$family,
      link = p$link,
      backend = "inla",
      return_code = TRUE,
      verbose = FALSE
    ))
    inla_class <- if (identical(inla_out$class, "admitted")) {
      cf <- inla_out$value$control_family
      used <- cf$link %||%
        cf$control.link$model %||%
        inla_default[[inla_out$value$family]]
      if (identical(used, p$requested)) "emitted" else paste0("changed:", used)
    } else {
      inla_out$class
    }

    auto_out <- .flg_call(flexybayes(
      y ~ 1 + (1 | g),
      data = d,
      family = p$family,
      link = p$link,
      backend = "auto",
      plan = TRUE,
      verbose = FALSE
    ))
    auto_class <- if (identical(auto_out$class, "admitted")) {
      paste0("plan:", auto_out$value$backend_chosen)
    } else {
      auto_out$class
    }

    expect_identical(
      brms_class,
      if (on_brms) "emitted" else "refusal:unsupported_family",
      label = paste(p$key, "on brms")
    )
    expect_identical(
      inla_class,
      if (on_inla) {
        "emitted"
      } else if (on_brms) {
        "refusal:inla_gate_refused"
      } else {
        "refusal:unsupported_family"
      },
      label = paste(p$key, "on INLA")
    )
    expect_identical(
      auto_class,
      if (on_inla) {
        "plan:inla"
      } else if (on_brms) {
        "plan:brms"
      } else {
        "refusal:unsupported_family"
      },
      label = paste(p$key, "under auto")
    )
    checked <- checked + 1L
  }
  expect_identical(checked, length(pair_list))
})
