# Blocks-format vm() and ped() on INLA refuse by name.
#
# A blocks carrier is one block-diagonal covariance, bdiag(V_1, ..., V_K), and
# a vm() or ped() term over it has one standard deviation. INLA's emitter fits
# each block as its own generic0 field with its own precision, which for two
# or more blocks is a different model, and the one-block fit, identical to the
# precision carrier of the same matrix, did not meet its brms oracle within
# the tolerance fixed for it. The gate refuses the carrier on
# INLA and names the precision carrier of the block-diagonal inverse, which
# fits the covariance with one standard deviation.
#
# Fixture: two unrelated families of 15 under one additive relationship
# matrix, three records per individual (90 rows, 30 levels).

.vc1_fixture <- function() {
  # Numerator relationship matrix by the tabular method. Parents precede
  # their offspring and 0 marks an unknown parent.
  a_matrix <- function(sire, dam) {
    n <- length(sire)
    A <- matrix(0, n, n)
    for (i in seq_len(n)) {
      s <- sire[[i]]
      d <- dam[[i]]
      A[i, i] <- 1 + if (s > 0L && d > 0L) 0.5 * A[s, d] else 0
      for (j in seq_len(i - 1L)) {
        v <- 0.5 * ((if (s > 0L) A[j, s] else 0) + (if (d > 0L) A[j, d] else 0))
        A[i, j] <- v
        A[j, i] <- v
      }
    }
    A
  }
  sire <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 3L, 3L, 3L, 3L, 5L, 5L, 9L)
  dam <- c(0L, 0L, 0L, 0L, 2L, 2L, 2L, 2L, 4L, 4L, 4L, 4L, 10L, 11L, 6L)
  family <- a_matrix(sire, dam)
  ids <- sprintf("g%02d", 1:30)
  A1 <- family
  A2 <- family
  dimnames(A1) <- list(ids[1:15], ids[1:15])
  dimnames(A2) <- list(ids[16:30], ids[16:30])
  A <- matrix(0, 30L, 30L, dimnames = list(ids, ids))
  A[1:15, 1:15] <- A1
  A[16:30, 16:30] <- A2

  data <- data.frame(geno = factor(rep(ids, each = 3L), levels = ids))
  withr::with_seed(20260916L, {
    u <- as.vector(t(chol(A)) %*% stats::rnorm(30L))
    data$y <- 10 + u[as.integer(data$geno)] + stats::rnorm(nrow(data), sd = 0.7)
  })
  list(
    data = data,
    A = A,
    Q = solve(A),
    blocks1 = list(A),
    blocks2 = list(A1, A2)
  )
}

test_that("a blocks carrier refuses on INLA by name, for vm() and ped(), with one block or two", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fx <- .vc1_fixture()
  terms <- list(
    vm = ~ vm(geno, cov = fb_cov(B, type = "blocks")),
    ped = ~ ped(geno, cov = fb_cov(B, type = "blocks"))
  )
  for (nm in names(terms)) {
    for (blocks in list(fx$blocks1, fx$blocks2)) {
      for (return_code in c(FALSE, TRUE)) {
        err <- tryCatch(
          flexybayes(
            fixed = y ~ 1,
            random = terms[[nm]],
            data = fx$data,
            known_matrices = list(B = blocks),
            backend = "inla",
            return_code = return_code,
            verbose = FALSE
          ),
          flexybayes_refusal = identity
        )
        label <- paste(
          nm,
          length(blocks),
          "block(s), return_code =",
          return_code
        )
        expect_s3_class(err, "flexybayes_lgm_random_term_type_inla")
        expect_identical(err$reason_code, "inla_gate_refused", label = label)
        expect_match(
          conditionMessage(err),
          "blocks carrier is refused on INLA",
          fixed = TRUE,
          label = label
        )
        expect_match(
          conditionMessage(err),
          "type = \"precision\"",
          fixed = TRUE,
          label = label
        )
      }
    }
  }
})

test_that("the precision carrier the refusal names fits the covariance with one standard deviation", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("Matrix")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fx <- .vc1_fixture()
  Q_named <- Matrix::bdiag(lapply(fx$blocks2, solve))
  terms <- list(
    vm = ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    ped = ~ ped(geno, cov = fb_cov(Q, type = "precision"))
  )
  for (nm in names(terms)) {
    fits <- lapply(list(named = Q_named, full = fx$Q), function(Q) {
      flexybayes(
        fixed = y ~ 1,
        random = terms[[nm]],
        data = fx$data,
        known_matrices = list(Q = Q),
        backend = "inla",
        verbose = FALSE
      )
    })
    vc <- fits$named$extras$variance_comps
    expect_identical(vc$component, fits$full$extras$variance_comps$component)
    expect_equal(sum(vc$component != "sigma"), 1L, label = nm)
    expect_true(all(is.finite(c(vc$estimate, vc$q2.5, vc$q97.5))), label = nm)
    expect_true(all(fits$named$extras$support_audit$inside), label = nm)
    # Q_named and solve(A) are the same matrix, and INLA repeats on one
    # thread to within 3e-6 of the value.
    full <- fits$full$extras$variance_comps
    expect_lte(
      max(
        abs(c(vc$q2.5, vc$q97.5) - c(full$q2.5, full$q97.5)) /
          abs(c(full$q2.5, full$q97.5))
      ),
      1e-5,
      label = nm
    )
  }
})
