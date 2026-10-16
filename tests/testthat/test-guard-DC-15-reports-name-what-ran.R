# Class guard DC-15: every report names what ran.
#
# Fits a set of small INLA models that between them reach the per-row and
# the aggregated routes, a Gaussian and a count family, a plain, an
# interaction and a structured random term, and reads every reporting
# surface back against the fitted object and the data:
#
#   * print() names only components the object carries (fca5e567: a $greta
#     slot on every fit);
#   * print() names the fit's own family and no other (3cf50d1d: every
#     aggregated fit printed 'aggregated-gaussian');
#   * ranef() and summary()$varcomp name groups and levels the data has, not
#     an index the engine call generated (8742b252, a15be4d6, 10905299);
#   * summary()$varcomp never calls a prior the engine default where the
#     fit's own prior record sent a package prior (770a19fc, afb569d9);
#   * summary()$engine is the engine the fit ran on.

.g15_data <- function() {
  set.seed(15)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", 1:6), each = 12)),
    trt = factor(rep(c("low", "mid", "high"), 24)),
    x = stats::rnorm(72)
  )
  u_site <- stats::rnorm(6)
  u_st <- stats::rnorm(18, sd = 0.5)
  d$y <- 1 +
    c(0, 0.5, -0.3)[d$trt] +
    0.4 * d$x +
    u_site[d$site] +
    u_st[(as.integer(d$site) - 1L) * 3L + as.integer(d$trt)] +
    stats::rnorm(72, sd = 0.5)
  d$cnt <- stats::rpois(
    72,
    exp(0.5 + 0.3 * as.integer(d$trt) + 0.3 * u_site[d$site])
  )
  d
}

.g15_fits <- function(d) {
  fit_one <- function(...) {
    suppressWarnings(suppressMessages(flexybayes(
      ...,
      data = d,
      backend = "inla",
      verbose = FALSE
    )))
  }
  list(
    gaussian_site = fit_one(y ~ trt + x, random = ~site),
    gaussian_interaction = fit_one(y ~ trt + x, random = ~ site + site:trt),
    gaussian_diag = fit_one(y ~ trt + x, random = ~ diag(trt):site),
    poisson_aggregated = fit_one(
      cnt ~ trt,
      random = ~site,
      family = "poisson",
      aggregate = TRUE
    ),
    poisson_per_row = fit_one(
      cnt ~ trt,
      random = ~site,
      family = "poisson",
      aggregate = FALSE
    )
  )
}

.G15_FAMILIES <- c(
  "gaussian",
  "poisson",
  "binomial",
  "gamma",
  "beta",
  "lognormal",
  "student",
  "negbinomial",
  "nbinomial",
  "weibull",
  "exponential"
)

# Words a group or component name may be built from: the data's column names,
# the levels of its factors, a column name joined to one of its levels as R
# names a coefficient (trthigh), and the names of the formula's special forms.
.g15_vocabulary <- function(d) {
  facs <- d[vapply(d, is.factor, logical(1L))]
  lev <- unlist(lapply(facs, levels), use.names = FALSE)
  coef_names <- unlist(
    lapply(names(facs), function(v) paste0(v, levels(facs[[v]]))),
    use.names = FALSE
  )
  c(
    names(d),
    lev,
    coef_names,
    "diag",
    "us",
    "idh",
    "at",
    "units",
    "sigma",
    "sd",
    "cor",
    "Intercept"
  )
}

.g15_unknown_tokens <- function(label, vocab) {
  tokens <- unlist(strsplit(label, "[^A-Za-z0-9.]+"))
  tokens <- tokens[nzchar(tokens)]
  setdiff(tokens, vocab)
}

test_that("DC-15 guard: print, summary and ranef name what the fit ran", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g15_data()
  fits <- .g15_fits(d)
  vocab <- .g15_vocabulary(d)
  level_vocab <- unique(c(
    unlist(
      lapply(d[vapply(d, is.factor, logical(1L))], levels),
      use.names = FALSE
    ),
    as.vector(outer(levels(d$site), levels(d$trt), paste, sep = ":")),
    as.vector(outer(levels(d$trt), levels(d$site), paste, sep = ":"))
  ))

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    expect_s3_class(fit, "flexybayes")
    printed <- utils::capture.output(print(fit))

    # Components the print names are components the object has.
    named <- unique(sub(
      "^\\$",
      "",
      regmatches(
        printed,
        regexpr("^\\s*\\$[A-Za-z_.][A-Za-z0-9_.]*", printed)
      )
    ))
    named <- trimws(sub("^\\s*\\$", "", named))
    expect_identical(
      setdiff(named, names(fit)),
      character(0),
      label = paste0(nm, ": components print() names that the fit lacks")
    )

    # The print names the fit's family and no other.
    fam <- stats::family(fit)$family
    words <- tolower(unlist(strsplit(
      paste(printed, collapse = " "),
      "[^A-Za-z]+"
    )))
    expect_true(
      fam %in% words,
      label = paste0(nm, ": print() names the family ", fam)
    )
    expect_identical(
      intersect(setdiff(.G15_FAMILIES, fam), words),
      character(0),
      label = paste0(nm, ": other family names in print()")
    )

    # ranef() groups and levels come from the data.
    r <- ranef(fit)
    expect_gt(length(r), 0L)
    for (g in names(r)) {
      expect_identical(
        .g15_unknown_tokens(g, vocab),
        character(0),
        label = paste0(
          nm,
          ": ranef() group '",
          g,
          "' names something the data lacks"
        )
      )
      lv <- as.character(r[[g]]$level)
      expect_identical(
        setdiff(lv, level_vocab),
        character(0),
        label = paste0(nm, ": ranef() levels of '", g, "' the data lacks")
      )
    }

    # summary()$varcomp names data terms, and a package prior is not called
    # the engine default.
    s <- summary(fit)
    vc <- s$varcomp
    for (comp in vc$component) {
      expect_identical(
        .g15_unknown_tokens(comp, vocab),
        character(0),
        label = paste0(
          nm,
          ": varcomp component '",
          comp,
          "' names something the data lacks"
        )
      )
    }
    slots <- prior_summary(fit)$slots
    if (is.data.frame(slots) && nrow(slots) && "prior" %in% names(vc)) {
      sent <- slots$target_key[
        slots$section == "random" & slots$source == "package"
      ]
      for (k in stats::na.omit(sent)) {
        row <- vc[vc$component == paste0("sd_", k), , drop = FALSE]
        if (nrow(row)) {
          expect_false(
            grepl("engine", row$prior[1], ignore.case = TRUE),
            label = paste0(
              nm,
              ": varcomp calls the package prior on sd_",
              k,
              " an engine default"
            )
          )
        }
      }
    }

    # summary() names the engine the fit ran on.
    expect_identical(s$engine, "inla", label = paste0(nm, ": summary()$engine"))
  }
})
