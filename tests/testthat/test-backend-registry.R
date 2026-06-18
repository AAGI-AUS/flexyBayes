# Backend registry (ADR 0031; v0.5.0 backend-axis recovery, Phase 1).
#
# The registry is additive and behaviour-neutral in Phase 1. These tests
# pin two things: (a) the registry's own shape + lock invariants, mirrored
# from the other four registries' test posture; and (b) the CONSISTENCY
# guard that makes the registry the single source of truth -- every
# backend name the current hard-coded dispatch references (the match.arg
# vocabularies + .routing_policy_table()) must be a registered name, and
# auto's default candidate set must match the registry's default_in_auto
# flags. This guard is what lets Phase 2 rewire dispatch onto the registry
# without silent drift.

test_that("backend registry is locked after .onLoad()", {
  expect_true(environmentIsLocked(flexyBayes:::.backend_registry))
  # A post-load registration must raise, not mutate.
  expect_error(
    flexyBayes:::.register_backend(
      name = "jags",
      status = "active",
      engine = "emit_jags",
      grammars = "brms",
      paradigm = "mcmc_gibbs",
      available_pkg = "rjags",
      default_in_auto = FALSE
    ),
    "locked"
  )
})

test_that("every backend entry carries the full field schema", {
  required <- c(
    "name",
    "status",
    "engine",
    "grammars",
    "paradigm",
    "available_pkg",
    "default_in_auto",
    "capability_predicate",
    "rename_to",
    "registered_in_adr"
  )
  for (nm in flexyBayes:::.registered_backend_names()) {
    e <- flexyBayes:::.lookup_backend(nm)
    expect_true(all(required %in% names(e)), info = nm)
    expect_identical(e$name, nm)
    expect_true(
      e$status %in% flexyBayes:::.BACKEND_STATUS_VOCABULARY,
      info = nm
    )
    expect_true(
      all(
        e$grammars %in%
          flexyBayes:::.BACKEND_GRAMMAR_VOCABULARY
      ),
      info = nm
    )
    expect_true(
      is.logical(e$default_in_auto) &&
        !is.na(e$default_in_auto),
      info = nm
    )
  }
})

test_that("the v0.5.0 population registers the expected backends", {
  expect_setequal(
    flexyBayes:::.registered_backend_names(),
    c("greta", "inla", "brms", "gretaR")
  )
  expect_identical(flexyBayes:::.lookup_backend("greta")$status, "active")
  expect_identical(flexyBayes:::.lookup_backend("inla")$status, "active")
  expect_identical(flexyBayes:::.lookup_backend("brms")$status, "active")
  expect_identical(flexyBayes:::.lookup_backend("gretaR")$status, "active")
  # (koine, the dormant 4th-backend slot, moved to flexyBayesOrchestra in the
  # lean-core split, 2026-06-06; it is no longer registered in the core.)
  # brms is retained as the engine label (the brms -> stan rename was
  # reversed 2026-05-31); no backend carries a rename target.
  expect_true(is.na(flexyBayes:::.lookup_backend("brms")$rename_to))
})

test_that("CONSISTENCY: routing-policy-table engines are all registered", {
  tbl <- flexyBayes:::.routing_policy_table()
  reg <- flexyBayes:::.registered_backend_names()
  # user_request values: backend names plus the meta-value "auto".
  ureq <- setdiff(unique(tbl$user_request), "auto")
  expect_true(
    all(ureq %in% reg),
    info = paste(
      "unregistered user_request:",
      paste(setdiff(ureq, reg), collapse = ", ")
    )
  )
  # chosen_backend values (drop the NA explicit-inla-refusal row).
  chosen <- stats::na.omit(unique(tbl$chosen_backend))
  expect_true(
    all(chosen %in% reg),
    info = paste(
      "unregistered chosen_backend:",
      paste(setdiff(chosen, reg), collapse = ", ")
    )
  )
})

test_that("CONSISTENCY: verb match.arg vocabularies are all registered", {
  reg <- flexyBayes:::.registered_backend_names()
  fx <- setdiff(eval(formals(flexyBayes::flexybayes)$backend), "auto")
  fb <- setdiff(eval(formals(flexyBayes::fb_brms)$backend), "auto")
  expect_true(
    all(fx %in% reg),
    info = paste("flexybayes:", paste(fx, collapse = ", "))
  )
  expect_true(
    all(fb %in% reg),
    info = paste("fb_brms:", paste(fb, collapse = ", "))
  )
})

test_that("CONSISTENCY: auto candidate set matches default_in_auto flags", {
  # The current dispatch candidate list lives in .resolve_routing();
  # the registry's default_in_auto-TRUE names must equal it.
  expect_setequal(
    flexyBayes:::.auto_default_backend_names(),
    c("greta", "inla", "gretaR")
  )
  # brms is deliberately NOT auto-default (opt-in only; Stan compile
  # latency would break the auto fast-path promise -- ADR 0024).
  expect_false(flexyBayes:::.lookup_backend("brms")$default_in_auto)
})

test_that(".available_backend_names() are active and installed", {
  avail <- flexyBayes:::.available_backend_names()
  for (nm in avail) {
    e <- flexyBayes:::.lookup_backend(nm)
    expect_identical(e$status, "active", info = nm)
    if (!is.na(e$available_pkg)) {
      expect_true(requireNamespace(e$available_pkg, quietly = TRUE), info = nm)
    }
  }
  # gretaR is activated; it is "available" iff the gretaR package is installed
  # (the version-floor probe is enforced downstream, in emit_gretaR()).
  expect_identical("gretaR" %in% avail, nzchar(system.file(package = "gretaR")))
})

test_that("capability predicates: greta universal, brms refuses structured-cov", {
  d <- data.frame(
    y = rnorm(40),
    env = factor(rep(letters[1:4], each = 10L)),
    geno = factor(rep(seq_len(10L), times = 4L))
  )
  fb_ri <- flexyBayes:::fb_from_asreml(y ~ env, random = ~geno, data = d)

  # A plain Gaussian random-intercept model: every active backend capable.
  expect_true(flexyBayes:::.capability_greta(fb_ri))
  expect_true(flexyBayes:::.capability_brms(fb_ri))
  expect_true(flexyBayes:::.capability_inla(fb_ri))
  expect_true(flexyBayes:::.backend_can_fit("greta", fb_ri)$ok)
  expect_true(flexyBayes:::.backend_can_fit("brms", fb_ri)$ok)
  expect_true(flexyBayes:::.backend_can_fit("inla", fb_ri)$ok)

  # Inject an asreml structured-covariance term: brms/Stan cannot
  # represent it; greta still can.
  fb_struct <- fb_ri
  fb_struct$random_terms <- c(
    fb_struct$random_terms,
    list(list(type = "fa", var = "env"))
  )
  expect_identical(
    flexyBayes:::.capability_brms(fb_struct),
    "stan_cannot_represent_structured_cov"
  )
  expect_true(flexyBayes:::.capability_greta(fb_struct))
  cf <- flexyBayes:::.backend_can_fit("brms", fb_struct)
  expect_false(cf$ok)
  expect_identical(cf$reason_code, "stan_cannot_represent_structured_cov")

  # Inject a low_rank approximation: .collect_approx() keys off
  # `approx_spec`, so a term carrying one is what brms/Stan must refuse.
  fb_lr <- fb_ri
  fb_lr$random_terms <- list(list(
    type = "spl",
    var = "x",
    approx_spec = list(scheme = "low_rank")
  ))
  expect_length(flexyBayes:::.collect_approx(fb_lr$random_terms), 1L)
  expect_identical(
    flexyBayes:::.capability_brms(fb_lr),
    "stan_cannot_represent_low_rank_approx"
  )
  expect_identical(
    flexyBayes:::.backend_can_fit("brms", fb_lr)$reason_code,
    "stan_cannot_represent_low_rank_approx"
  )
  expect_true(flexyBayes:::.capability_greta(fb_lr))
})

test_that("flexybayes(backend='brms') allows dense vm/ped but refuses non-dense-able structured carriers", {
  skip_if_not_installed("brms")
  # Genomics expansion (G1): vm() / ped() with an exact dense-able
  # carrier now reach brms via its native known-covariance group term
  # (1 | gr(var, cov = K)) -- brms Cholesky-factors K internally. The
  # remaining structured-covariance carriers (block-diagonal, low-rank,
  # and the fa / us / ar1 terms) still have no lossless Stan translation
  # and are refused structurally rather than failing mid-emit.
  d <- data.frame(
    y = rnorm(40),
    env = factor(rep(letters[1:4], each = 10L)),
    geno = factor(rep(seq_len(10L), times = 4L))
  )
  G <- diag(10L)
  dimnames(G) <- list(as.character(seq_len(10L)), as.character(seq_len(10L)))

  # Dense GBLUP now generates Stan code (the gr(cov = K) route).
  code <- flexybayes(
    y ~ env,
    random = ~ vm(geno, G),
    data = d,
    known_matrices = list(G = G),
    backend = "brms",
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl("Lcov", code, fixed = TRUE))

  # A block-diagonal vm() carrier is still refused at the capability layer.
  Bs <- list(diag(5L), diag(5L))
  expect_error(
    flexybayes(
      y ~ env,
      random = ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
      data = d,
      known_matrices = list(Bs = Bs),
      backend = "brms",
      verbose = FALSE,
      mcmc_verbose = FALSE
    ),
    regexp = "structured-covariance term"
  )
  # brms is reachable from flexybayes() now (vocabulary widened).
  expect_true("brms" %in% eval(formals(flexyBayes::flexybayes)$backend))
})

test_that(".register_backend() refuses unknown status / grammar / duplicate", {
  # These run against the locked registry, so the vocabulary checks
  # (validated before the lock check) are what must fire first.
  expect_error(
    flexyBayes:::.register_backend(
      name = "x",
      status = "experimental",
      engine = NA_character_,
      grammars = "brms",
      paradigm = "p",
      available_pkg = NA_character_,
      default_in_auto = FALSE
    ),
    "unknown status"
  )
  expect_error(
    flexyBayes:::.register_backend(
      name = "x",
      status = "active",
      engine = "e",
      grammars = "pymc",
      paradigm = "p",
      available_pkg = NA_character_,
      default_in_auto = FALSE
    ),
    "unknown grammar"
  )
})
