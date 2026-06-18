# Parser-derived preflight tests (audit P1.2, 2026-05-23). The
# hand-built IR tests in test-fb-preflight.R exercise the byte-formula
# math; these tests exercise the END-TO-END path: real formula -> live
# parser (parse_formula.R + .enrich()) -> live preflight. They catch
# the specific class of bug the audit named: parser-emitted type
# strings that the preflight estimator did not match, silently routing
# to the fallthrough and under-estimating by a factor of `k` or `df`.
#
# Uses fb_from_asreml() rather than flexybayes() to skip the backend
# dispatch/fit path -- we only need the parsed IR.

test_that("preflight on parsed s(x, k = 6) matches object.size() of the basis", {
  skip_if_not_installed("mgcv")

  set.seed(2026)
  n <- 5000L
  dx <- data.frame(x = sort(stats::runif(n, 0, 10)))
  dx$y <- sin(dx$x) + stats::rnorm(n, 0, 0.2)

  fb <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ s(x, k = 6L),
    data = dx
  )

  # Confirm the live IR shape -- audit's central observation.
  smooth_idx <- which(vapply(
    fb$random_terms,
    function(t) identical(t$type, "smooth_mgcv"),
    logical(1L)
  ))
  expect_length(smooth_idx, 1L)
  smooth_term <- fb$random_terms[[smooth_idx]]

  # .enrich() populates term$X (the N x k basis) and term$k after
  # mgcv::smoothCon() absorbs the identifiability constraint, which
  # reduces user-supplied k = 6 to effective k = 5.
  expect_true(!is.null(smooth_term$X))
  expect_equal(smooth_term$k, 5L)

  # Reference: object.size() of the basis matrix that .enrich() built.
  actual_basis_bytes <- as.numeric(utils::object.size(smooth_term$X))

  ds <- .fb_dataset(dx)
  pf <- .fb_preflight(fb, ds)

  smooth_entries <- Filter(
    function(e) identical(e$term_kind, "random_smooth"),
    pf$per_term_estimate
  )
  expect_length(smooth_entries, 1L)

  # Tolerance per ADR 0021: per-term estimate accuracy <= 5% vs the
  # actually-allocated representation. The preflight estimator
  # uses 8 * N * k + 8 * k + 256; object.size() includes a few extra
  # bytes for SEXP header + dimnames + class attribute.
  expect_equal(
    smooth_entries[[1L]]$design_memory_bytes,
    actual_basis_bytes,
    tolerance = 0.05
  )

  # Sanity: the new estimate is much larger than the old fallthrough
  # (12N + 16k + 256 = 60,336 at N=5000, k=5) which was the audit-
  # confirmed under-estimate.
  expect_gt(smooth_entries[[1L]]$design_memory_bytes, 100000)
})

test_that("preflight on parsed spl(x) matches splines::bs() basis size", {
  set.seed(2026)
  n <- 5000L
  dx <- data.frame(x = stats::runif(n, 0, 10))
  dx$y <- dx$x + stats::rnorm(n, 0, 0.2)

  fb <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ spl(x),
    data = dx
  )

  spline_idx <- which(vapply(
    fb$random_terms,
    function(t) identical(t$type, "spline"),
    logical(1L)
  ))
  expect_length(spline_idx, 1L)

  # Reference: the basis that codegen would build at fit time.
  # codegen hardcodes splines::bs(x_std, df = 8, degree = 3,
  # intercept = FALSE).
  x_std <- as.numeric(scale(dx$x))
  B_ref <- splines::bs(x_std, df = 8, degree = 3, intercept = FALSE)
  actual_basis_bytes <- as.numeric(utils::object.size(B_ref))

  ds <- .fb_dataset(dx)
  pf <- .fb_preflight(fb, ds)

  spline_entries <- Filter(
    function(e) identical(e$term_kind, "random_spline"),
    pf$per_term_estimate
  )
  expect_length(spline_entries, 1L)

  expect_equal(
    spline_entries[[1L]]$design_memory_bytes,
    actual_basis_bytes,
    tolerance = 0.05
  )

  # 8 * 5000 * 8 + 8 * 8 + 256 = 320,320 -- the audit-confirmed
  # expected number, vs the v0.3.0 fallthrough 60,336.
  expect_equal(
    spline_entries[[1L]]$design_memory_bytes,
    8 * 5000 * 8 + 8 * 8 + 256
  )
})

test_that("preflight on parsed y ~ f1 * f2 sizes the factor:factor block", {
  set.seed(2026)
  n <- 1000L
  df <- data.frame(
    y = rnorm(n),
    f1 = factor(rep(letters[1:10], n / 10L)),
    f2 = factor(rep(LETTERS[1:5], each = n / 5L))
  )

  fb <- fb_from_asreml(
    fixed = y ~ f1 * f2,
    data = df
  )

  # The parser produces three fixed-effect terms beyond the intercept:
  # f1, f2, f1:f2 (factor_interaction).
  fi_idx <- which(vapply(
    fb$fixed_terms,
    function(t) identical(t$type, "factor_interaction"),
    logical(1L)
  ))
  expect_length(fi_idx, 1L)

  # Reference: the dense treatment-coded interaction block sized
  # as a clean numeric matrix + coefficient vector (matches the
  # idiom used by the other estimate-accuracy tests in
  # test-fb-preflight.R; excludes dimnames overhead from
  # model.matrix() which is not part of the design-representation
  # contract).
  L1 <- nlevels(df$f1)
  L2 <- nlevels(df$f2)
  n_dummies <- (L1 - 1L) * (L2 - 1L)
  actual_inter_bytes <-
    as.numeric(utils::object.size(matrix(0, nrow(df), n_dummies))) +
    as.numeric(utils::object.size(numeric(n_dummies)))

  ds <- .fb_dataset(df)
  pf <- .fb_preflight(fb, ds)

  fi_entries <- Filter(
    function(e) identical(e$term_kind, "fixed_factor_interaction"),
    pf$per_term_estimate
  )
  expect_length(fi_entries, 1L)

  expect_equal(
    fi_entries[[1L]]$design_memory_bytes,
    actual_inter_bytes,
    tolerance = 0.05
  )

  # Sanity: vastly larger than the v0.3.0 fallthrough (8N + 8 + 256
  # = 8,264 bytes).
  expect_gt(fi_entries[[1L]]$design_memory_bytes, 100000)
})
