# Tier-3 stress test for the Stage 2 MVP preflight layer
# (ADR 0021 / v0.3.0). Asserts the lazy-O(terms) property: a 1e8-row
# metadata-only descriptor can be preflighted in <= 60 s wall time
# with < 200 MB RSS delta -- because preflight reads only
# fb_dataset metadata (n_rows, dictionaries, col_types) and never
# materialises a 1e8-row data.table.
#
# Gated by `FLEXYBAYES_RUN_STRESS=true` because it allocates the
# 1e8-element factor-dictionary character vector (a ~1.6 GB R object
# on most platforms; safe to construct but a no-go on small CI
# runners). Run locally via:
#   FLEXYBAYES_RUN_STRESS=true R -e 'devtools::test(filter = "preflight-stress")'

skip_if_not_stress <- function() {
  if (!identical(Sys.getenv("FLEXYBAYES_RUN_STRESS"), "true")) {
    testthat::skip("FLEXYBAYES_RUN_STRESS != \"true\"; skip Tier-3 stress.")
  }
}

# Process-resident-set-size probe. Uses `ps -o rss=` on Unix (returns
# KB); returns NA on platforms where the probe is unavailable so the
# RSS assertion can be skipped cleanly.
.process_rss_mb <- function() {
  sys <- Sys.info()[["sysname"]]
  if (sys %in% c("Darwin", "Linux")) {
    out <- tryCatch(
      system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE),
      error = function(e) character(0),
      warning = function(w) character(0)
    )
    if (length(out) >= 1L) {
      val <- suppressWarnings(as.numeric(trimws(out)))
      if (!is.na(val) && val > 0) return(val / 1024) # KB -> MB
    }
  }
  NA_real_
}

test_that("1e8-row metadata-only preflight: lazy O(terms) wall time + RSS", {
  skip_if_not_stress()

  # Constructable: the dictionary is 1e5 group levels (1.6 MB-ish);
  # the n_rows = 1e8 is just an integer scalar -- no row-wise allocation.
  t0_dict <- Sys.time()
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e8,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1e5)))
  )
  t1_dict <- Sys.time()
  expect_true(.fb_dataset_is_metadata(ds))
  expect_null(ds$data)
  expect_identical(ds$n_rows, 100000000L)
  # Dataset construction itself is fast on metadata-only path
  expect_lte(as.numeric(difftime(t1_dict, t0_dict, units = "secs")), 10)

  # IR representative of the worked workflow: random intercept on the
  # 1e5-level grouping factor + a fixed numeric covariate. Avoids
  # smooth so the test exercises the lazy path through every estimator
  # branch without dense_baseline triggering by itself.
  fb <- structure(
    list(
      response = "y",
      family = "gaussian",
      link = "identity",
      intercept = TRUE,
      fixed_terms = list(list(type = "numeric", var = "x")),
      random_terms = list(list(type = "simple", var = "g", var_n = 1e5)),
      rcov_terms = list(),
      addition_terms = list(),
      priors = list(),
      data_summary = list(n = 1e8),
      capabilities = character(),
      source = "stress_test"
    ),
    class = c("fb_terms", "list")
  )

  rss_before <- .process_rss_mb()
  t0_pf <- Sys.time()
  pf <- .fb_preflight(fb, ds, memory_ceiling_gb = 32)
  t1_pf <- Sys.time()
  rss_after <- .process_rss_mb()

  # Wall time gate: the lazy estimator runs in O(terms) so the 1e8
  # n_rows does not enter the time complexity -- only the per-term
  # arithmetic does. ADR target: <= 60 s.
  pf_wall_s <- as.numeric(difftime(t1_pf, t0_pf, units = "secs"))
  message(sprintf("preflight wall time at n_rows = 1e8: %.3f s", pf_wall_s))
  expect_lte(pf_wall_s, 60)

  # RSS delta gate: < 200 MB. Skip the assertion if the probe is
  # unavailable on this platform (Windows / minimal containers).
  if (!is.na(rss_before) && !is.na(rss_after)) {
    rss_delta <- rss_after - rss_before
    message(sprintf("preflight RSS delta at n_rows = 1e8: %.1f MB", rss_delta))
    expect_lte(rss_delta, 200)
  } else {
    message(
      "RSS probe unavailable on this platform; ",
      "skipped the RSS assertion."
    )
  }

  # Sanity: the per-term entries carry the correct n_rows-driven
  # scaling and the aggregation flag matches the in-scope envelope.
  expect_identical(pf$n_rows, 1e8)
  expect_identical(length(pf$per_term_estimate), 2L)
  expect_true(pf$per_term_estimate$x$aggregated_likelihood_candidate)
  expect_true(pf$per_term_estimate[["(1 | g)"]]$aggregated_likelihood_candidate)
  # Random-intercept estimate is dominated by the integer index column:
  # 4 * 1e8 = 4e8 bytes. Allow a small overhead slack.
  expect_gte(pf$per_term_estimate[["(1 | g)"]]$design_memory_bytes, 4e8 - 1024)
})

test_that("1e8-row preflight: skipped cleanly when FLEXYBAYES_RUN_STRESS != true", {
  # Default-off: a routine devtools::test() invocation skips the
  # heavy path without spuriously reporting it as a failure.
  if (identical(Sys.getenv("FLEXYBAYES_RUN_STRESS"), "true")) {
    succeed() # if running stress, the prior test_that covers this
    return(invisible())
  }
  expect_silent({
    # If the gate wraps correctly, calling skip_if_not_stress() inside
    # a try-catch yields nothing observable in capture.output() or in
    # the message stream.
    res <- tryCatch(skip_if_not_stress(), skip = function(s) "skipped")
    expect_identical(res, "skipped")
  })
})
