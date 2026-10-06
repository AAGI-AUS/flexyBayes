# DEF-197: the flexybayes_stream() help states the prior a fit with
# `prior = NULL` runs under, and a streamed fit runs under that prior.

# The text of one argument on the flexybayes_stream() help page, from the
# source tree when the tests run there and from the installed help otherwise.
.sdp_arg_text <- function(arg) {
  root <- testthat::test_path("..", "..")
  rd_file <- file.path(root, "man", "flexybayes_stream.Rd")
  rd <- if (file.exists(rd_file)) {
    tools::parse_Rd(rd_file)
  } else {
    tools::Rd_db("flexyBayes")[["flexybayes_stream.Rd"]]
  }
  tags <- vapply(rd, attr, character(1L), which = "Rd_tag")
  args <- rd[[which(tags == "\\arguments")]]
  for (item in args) {
    if (identical(attr(item, "Rd_tag"), "\\item") &&
      identical(paste(unlist(item[[1L]]), collapse = ""), arg)) {
      return(gsub("\\s+", " ", trimws(paste(unlist(item[[2L]]), collapse = ""))))
    }
  }
  NA_character_
}

test_that("DEF-197: the help names the package default, not the backend's", {
  txt <- .sdp_arg_text("prior")
  expect_false(is.na(txt))
  expect_false(grepl("backend default", txt, fixed = TRUE))
  expect_match(txt, "half_normal(scale = 3)", fixed = TRUE)
  expect_match(txt, "half_normal(scale = 5)", fixed = TRUE)
  expect_match(txt, "fixed effects keep INLA's default priors", fixed = TRUE)
})

test_that("DEF-197: streamed fits with no prior run under the priors the help states", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  set.seed(4)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", 1:8), each = 12)),
    trt = factor(rep(c("a", "b", "c"), 32))
  )
  d$cnt <- stats::rpois(96, exp(0.5 + 0.2 * as.integer(d$trt) + stats::rnorm(8, sd = 0.5)[d$site]))
  fit <- suppressWarnings(flexybayes_stream(
    cnt ~ trt, random = ~site, source = d, family = "poisson", verbose = FALSE
  ))
  d$k <- 10L
  d$s <- pmin(d$cnt, 10L)
  fit_binomial <- suppressWarnings(flexybayes_stream(
    s ~ trt, random = ~site, source = d, family = "binomial", trials = "k", verbose = FALSE
  ))
  for (case in list(list(fit = fit, scale = "s=3;"), list(fit = fit_binomial, scale = "s=5;"))) {
    slots <- prior_summary(case$fit)$slots
    sd_slot <- slots[slots$section == "random", , drop = FALSE]
    expect_identical(sd_slot$source, "package")
    expect_match(sd_slot$readback, case$scale, fixed = TRUE)
    expect_true(all(slots$source[slots$section == "fixed"] == "engine_default"))
  }
})
