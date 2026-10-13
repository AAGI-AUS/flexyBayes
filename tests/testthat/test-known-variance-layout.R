# How a known variance or covariance is carried with the rows of the data:
# a column of variances, a row index into a covariance matrix, and a zero on
# a row the design grid adds.

fbkl_trials <- function(seed = 7L, n = 40L, n_cluster = 12L) {
  .fb_with_seed(seed, {
    cluster <- factor(sprintf("c%02d", rep(seq_len(n_cluster),
                                             length.out = n)))
    x <- matrix(stats::rnorm(n * 6L), n, 6L,
                dimnames = list(NULL, paste0("x", 1:6)))
    a <- matrix(stats::rnorm(n * n, 0, 0.05), n, n)
    k <- diag(stats::runif(n, 0.05, 0.3)) + tcrossprod(a)
    l <- stats::rnorm(n_cluster, 0, 0.6)[as.integer(cluster)]
    m <- 2 + 0.8 * x[, 1L] + l + stats::rnorm(n, 0, 0.4) +
      drop(t(chol(k)) %*% stats::rnorm(n))
    list(d = data.frame(m = m, cluster = cluster, x), k = k)
  })
}

test_that("a matrix is carried as its diagonal and a row index, and follows its rows", {
  tr <- fbkl_trials()
  fb <- fb_from_asreml(m ~ 1, random = ~cluster, data = tr$d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  attached <- .fb_known_variance_attach(fb, tr$d, tr$k, aggregate = "auto")
  expect_identical(attached$data$flexybayes_known_variance, diag(tr$k))
  expect_identical(attached$data$flexybayes_known_row, seq_len(nrow(tr$d)))
  expect_identical(
    .fb_known_covariance_values(attached$fb, attached$data),
    tr$k
  )
  # Dropping and reordering rows takes the matching rows and columns. A row
  # the design grid adds has no index and no known part.
  rows <- c(9L, 2L, 30L)
  kept <- attached$data[rows, ]
  expect_identical(
    .fb_known_covariance_values(attached$fb, kept),
    tr$k[rows, rows]
  )
  grown <- rbind(kept, kept[1L, ])
  grown <- .fb_known_variance_fill(attached$fb, grown, 1L)
  got <- .fb_known_covariance_values(attached$fb, grown)
  expect_identical(got[1:3, 1:3], tr$k[rows, rows])
  expect_true(all(got[4L, ] == 0) && all(got[, 4L] == 0))
  expect_identical(
    .fb_known_variance_reissue(attached$fb, kept),
    tr$k[rows, rows]
  )
  # A vector re-issues as the vector.
  vec <- .fb_known_variance_attach(fb, tr$d, diag(tr$k), aggregate = "auto")
  expect_identical(.fb_known_variance_reissue(vec$fb, vec$data), diag(tr$k))
})

test_that("a plot the design grid adds carries no known variance", {
  # A plot absent from a field's grid is added with no response. The data
  # the engine would receive are read where flexybayes() hands them over.
  g <- .fb_with_seed(5L, {
    grid <- expand.grid(col = 1:5, row = 1:6)
    out <- data.frame(row = factor(grid$row), col = factor(grid$col))
    out$y <- stats::rnorm(nrow(out))
    out$v <- stats::runif(nrow(out), 0.1, 0.5)
    out[-7L, ]
  })
  local_mocked_bindings(.dispatch_backend = function(fb, data, ...) data)
  handed <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = g,
                       known_variance = g$v, backend = "inla",
                       verbose = FALSE)
  added <- is.na(handed$y)
  expect_identical(sum(added), 1L)
  expect_identical(handed$flexybayes_known_variance[added], 0)
  key <- paste(handed$row, handed$col)[!added]
  expect_identical(
    handed$flexybayes_known_variance[!added],
    g$v[match(key, paste(g$row, g$col))]
  )
})

