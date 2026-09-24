# Every random-effect index the INLA emitter generates has a specification.
#
# emit_inla() indexes a structured random term by an integer column it builds,
# and INLA names the summary.random entry after that column and numbers its
# rows 1..n. ranef() can only print a group and levels a reader recognises if
# something maps each generated index back, and that mapping used to grow one
# route at a time: the level map covered interactions and the correlated slope
# and nothing else, so diag(), fa(), vm(), ped() and the AR1 fields reported
# levels "1", "2", ... (review finding N5).
#
# The tests below take the set of generated indices from the EMITTER -- the
# term types .build_inla_formula() branches on, and the f() terms it writes
# for a fixture of each type -- and require a specification for every one.
# Reading the set off the specification registry instead would only prove the
# registry agrees with itself.

# The random-term types the formula builder lowers, read off its own source:
# every `term$type` literal it compares against inside the loop over
# `fb$random_terms`. A new branch there is a new index type.
.ris_emitter_types <- function(fun) {
  is_type_ref <- function(x) {
    is.call(x) && identical(x[[1L]], as.name("$")) &&
      identical(x[[2L]], as.name("term")) &&
      identical(x[[3L]], as.name("type"))
  }
  found <- character(0)
  collect <- function(x) {
    head <- x[[1L]]
    if (identical(head, as.name("identical")) && length(x) == 3L &&
          is_type_ref(x[[2L]]) && is.character(x[[3L]])) {
      found <<- c(found, x[[3L]])
    }
    if (identical(head, as.name("%in%")) && is_type_ref(x[[2L]])) {
      found <<- c(found, eval(x[[3L]], baseenv()))
    }
    if (identical(head, as.name("switch")) && is_type_ref(x[[2L]])) {
      nms <- names(x)[-c(1L, 2L)]
      found <<- c(found, nms[nzchar(nms)])
    }
    for (i in seq_along(x)[-1L]) {
      if (is.call(x[[i]])) collect(x[[i]])
    }
  }
  find_loop <- function(x) {
    if (identical(x[[1L]], as.name("for")) &&
          identical(x[[3L]], quote(fb$random_terms))) {
      collect(x[[4L]])
      return(invisible(NULL))
    }
    for (i in seq_along(x)[-1L]) {
      if (is.call(x[[i]])) find_loop(x[[i]])
    }
  }
  find_loop(body(fun))
  unique(found)
}

# The indices an emitted formula carries: the first argument of every f()
# term, which is what INLA names a summary.random entry after, and the
# group = / replicate = indices, which the same f() term reads without a
# table of their own.
.ris_emitted_indices <- function(form) {
  out <- list(tables = character(0), carried = character(0))
  walk <- function(x) {
    if (identical(x[[1L]], as.name("f"))) {
      args <- as.list(x)[-1L]
      nms <- names(args)
      if (is.null(nms)) nms <- rep("", length(args))
      first <- args[[which(!nzchar(nms))[1L]]]
      out$tables <<- c(out$tables, as.character(first))
      for (k in intersect(c("group", "replicate"), nms)) {
        out$carried <<- c(out$carried, as.character(args[[k]]))
      }
      return(invisible(NULL))
    }
    for (i in seq_along(x)[-1L]) {
      if (is.call(x[[i]])) walk(x[[i]])
    }
  }
  walk(form[[3L]])
  out
}

# One fixture per random-term type, built through the public ingest adapters.
.ris_fixtures <- function() {
  set.seed(7L)
  met <- expand.grid(
    gen = factor(paste0("G", 1:6)),
    env = factor(c("E1", "E2", "E3")),
    rep = factor(c("R1", "R2"))
  )
  met$y <- stats::rnorm(nrow(met))
  met$x <- stats::rnorm(nrow(met))
  met$z <- stats::runif(nrow(met), 0, 10)

  levs <- paste0("g", 1:5)
  kin <- diag(5) + 0.2
  dimnames(kin) <- list(levs, levs)
  gen_tab <- data.frame(geno = factor(rep(levs, 4)), y = stats::rnorm(20))
  blocks <- list(diag(2) + 0.1, diag(3) + 0.2)

  grid <- expand.grid(
    col = factor(paste0("C", 1:3)),
    row = factor(paste0("W", 1:4)),
    trial = factor(c("TA", "TB"))
  )
  grid$y <- stats::rnorm(nrow(grid))
  series <- expand.grid(
    t = factor(paste0("T", 1:5)),
    plot = factor(c("P1", "P2"))
  )
  series$y <- stats::rnorm(nrow(series))

  asr <- function(random, data, known_matrices = list()) {
    list(
      fb = fb_from_asreml(y ~ 1, random = random, data = data,
                          known_matrices = known_matrices),
      data = data,
      known_matrices = known_matrices
    )
  }
  # The grammar parses both id() and ide() to the IR type "ide", and the
  # formula builder still branches on "id", so that term is built directly.
  id_term <- asr(~ ide(gen), met)
  id_term$fb$random_terms[[1L]]$type <- "id"

  list(
    simple = asr(~gen, met),
    ide = asr(~ ide(gen), met),
    id = id_term,
    spline = asr(~ spl(z), met),
    nested = asr(~ rep:gen, met),
    combo = asr(~ rep:env:gen, met),
    uncor_slope = list(fb = fb_from_brms(y ~ x + (x || gen), data = met),
                       data = met, known_matrices = list()),
    cor_slope = list(fb = fb_from_brms(y ~ x + (x | gen), data = met),
                     data = met, known_matrices = list()),
    at_simple = asr(~ diag(env):gen, met),
    us_gxe = asr(~ us(env):gen, met),
    fa_gxe = asr(~ fa(env, 1):gen, met),
    vm = asr(~ vm(geno, cov = fb_cov(Qk, type = "precision")), gen_tab,
             list(Qk = solve(kin))),
    ped = asr(~ ped(geno, Ak, use_sparse_precision = TRUE), gen_tab,
              list(Ak = solve(kin))),
    vm_blocks = asr(~ vm(geno, cov = fb_cov(Bk, type = "blocks")), gen_tab,
                    list(Bk = blocks)),
    ar1 = asr(~ ar1(t), series[series$plot == "P1", ]),
    ar1_grouped = asr(~ ar1(t):plot, series),
    ar1_spatial = asr(~ ar1(row):ar1(col), grid[grid$trial == "TA", ]),
    ar1_replicated = asr(~ id(trial):ar1(row):ar1(col), grid)
  )
}

test_that("every term type the INLA formula builder lowers has a specification", {
  types <- .ris_emitter_types(.build_inla_formula)
  # A walk that found nothing would pass the comparison below vacuously.
  expect_gt(length(types), 5L)
  registry <- .fb_inla_ranef_registry()
  expect_setequal(names(registry), types)
  for (tt in names(registry)) {
    entry <- registry[[tt]]
    ok <- is.function(entry) ||
      (is.character(entry) && length(entry) == 1L && nzchar(entry))
    expect_true(ok, info = tt)
  }
})

test_that("the fixture set covers every term type the builder lowers", {
  fixtures <- .ris_fixtures()
  covered <- unique(unlist(lapply(fixtures, function(fx) {
    vapply(fx$fb$random_terms, function(t) t$type, character(1L))
  })))
  expect_true(
    all(.ris_emitter_types(.build_inla_formula) %in% covered),
    info = paste(
      setdiff(.ris_emitter_types(.build_inla_formula), covered),
      collapse = ", "
    )
  )
})

test_that("every index the emitter writes has exactly one specification", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  for (nm in names(.ris_fixtures())) {
    fx <- .ris_fixtures()[[nm]]
    code <- suppressMessages(emit_inla(
      fx$fb,
      data = fx$data,
      known_matrices = fx$known_matrices,
      return_code = TRUE,
      verbose = FALSE
    ))
    emitted <- .ris_emitted_indices(code$formula)
    specs <- .fb_inla_index_specs(
      fx$fb,
      data = fx$data,
      known_matrices = fx$known_matrices
    )
    expect_false(anyDuplicated(names(specs)) > 0L, info = nm)

    for (idx in emitted$tables) {
      if (idx %in% names(specs)) {
        s <- specs[[idx]]
        expect_true(nzchar(s$label), info = paste(nm, idx))
        expect_false(is.null(s$levels), info = paste(nm, idx))
        expect_identical(length(s$levels), length(s$rows),
                         info = paste(nm, idx))
        expect_true(max(s$rows) <= s$n_field, info = paste(nm, idx))
        expect_false(any(grepl("^[0-9]+$", s$levels)), info = paste(nm, idx))
      } else {
        # No specification is correct only where INLA was handed a column
        # the user supplied, whose values are already the levels.
        expect_true(idx %in% names(fx$data), info = paste(nm, idx))
      }
    }
    # And no specification for an index the emitter does not write.
    expect_true(all(names(specs) %in% emitted$tables),
                info = paste(nm, paste(setdiff(names(specs), emitted$tables),
                                       collapse = ", ")))
    # A group = or replicate = index has no table; the specification of the
    # table that carries it must say so.
    carried <- unlist(lapply(specs, function(s) s$carries))
    expect_setequal(carried %||% character(0), emitted$carried)
  }
})

test_that("the levels a specification reports are the data's, in index order", {
  # The expected levels are composed here from the data columns, in the
  # order each emitter index is built in: interaction levels are sorted
  # with the first variable slowest, and an unstructured or grid field
  # enumerates its first index fastest.
  fixtures <- .ris_fixtures()
  met <- fixtures$us_gxe$data
  gen <- levels(met$gen)
  env <- levels(met$env)

  s <- .fb_inla_index_specs(fixtures$us_gxe$fb, data = met)
  expect_identical(s$gen_by_env_id$label, "gen:env")
  expect_identical(
    s$gen_by_env_id$levels,
    paste0(rep(gen, times = 3L), ":", rep(env, each = 6L))
  )

  s <- .fb_inla_index_specs(fixtures$at_simple$fb, data = met)
  expect_identical(s$gen_at_env_E2$label, "gen (env = E2)")
  expect_identical(s$gen_at_env_E2$levels, gen)

  s <- .fb_inla_index_specs(
    fixtures$vm_blocks$fb,
    data = fixtures$vm_blocks$data,
    known_matrices = fixtures$vm_blocks$known_matrices
  )
  expect_identical(s$geno_id_block_1$levels, c("g1", "g2"))
  expect_identical(s$geno_id_block_2$levels, c("g3", "g4", "g5"))

  grid <- fixtures$ar1_replicated$data
  s <- .fb_inla_index_specs(fixtures$ar1_replicated$fb, data = grid)
  cells <- expand.grid(
    row = levels(grid$row), col = levels(grid$col), trial = levels(grid$trial),
    stringsAsFactors = FALSE
  )
  expect_identical(s$row_id$label, "id(trial):ar1(row):ar1(col)")
  expect_identical(
    s$row_id$levels,
    paste(cells$trial, cells$row, cells$col, sep = ":")
  )
  expect_setequal(s$row_id$carries, c("col_id", "trial_id"))
})


# --- the accessor, on engine tables of a known shape --------------------- #
#
# A table shaped as INLA returns it -- `ID` 1..n and one row per index value
# -- so the accessor's use of a specification is checked without a fit. The
# row count is taken from the data, not from the specification.

.ris_mock_fit <- function(fb, data, tables, known_matrices = list(),
                          attach_terms = TRUE) {
  structure(
    list(
      inla = list(summary.random = tables),
      fb = fb,
      data = data,
      level_labels = list(),
      extras = list(
        fb_terms = if (attach_terms) fb else NULL,
        call_info = list(known_matrices = known_matrices)
      )
    ),
    class = c("flexybayes_inla", "flexybayes")
  )
}

.ris_table <- function(n) {
  data.frame(
    ID = seq_len(n),
    mean = seq_len(n) / 10,
    sd = rep(0.5, n),
    `0.025quant` = seq_len(n) / 10 - 1,
    `0.975quant` = seq_len(n) / 10 + 1,
    check.names = FALSE
  )
}

test_that("ranef() puts the data's levels on a generated index without moving a number", {
  fx <- .ris_fixtures()$at_simple
  n_gen <- nlevels(fx$data$gen)
  tables <- list(gen_at_env_E1 = .ris_table(n_gen))
  r <- .fb_random_inla(.ris_mock_fit(fx$fb, fx$data, tables))
  expect_named(r, "gen (env = E1)")
  expect_identical(r[[1L]]$level, levels(fx$data$gen))
  expect_identical(r[[1L]]$estimate, tables[[1L]]$mean)
  expect_identical(r[[1L]]$std.error, tables[[1L]]$sd)
  expect_identical(r[[1L]]$conf.low, tables[[1L]]$`0.025quant`)
  expect_identical(r[[1L]]$conf.high, tables[[1L]]$`0.975quant`)
})

test_that("a table whose length disagrees with its specification is reported", {
  # A table of another length is not the field the term indexed. Its rows
  # are neither relabelled nor cut to fit, and the disagreement is said out
  # loud rather than leaving integer levels to pass as a result.
  fx <- .ris_fixtures()$at_simple
  n_gen <- nlevels(fx$data$gen)
  tables <- list(gen_at_env_E1 = .ris_table(n_gen + 1L))
  expect_warning(
    r <- .fb_random_inla(.ris_mock_fit(fx$fb, fx$data, tables)),
    "gen (env = E1)",
    fixed = TRUE
  )
  expect_identical(nrow(r[[1L]]), n_gen + 1L)
  expect_identical(r[[1L]]$estimate, tables[[1L]]$mean)
})

test_that("a fit made by emit_inla() directly is relabelled from its own IR", {
  # emit_inla() stores the IR as `fb`; only flexybayes() also attaches it as
  # `extras$fb_terms`. us(f):g and (x | g) are refused at the gate, so a
  # direct emit is the one route by which they reach ranef() on INLA.
  fx <- .ris_fixtures()$us_gxe
  n <- nlevels(fx$data$gen) * nlevels(fx$data$env)
  tables <- list(gen_by_env_id = .ris_table(n))
  r <- .fb_random_inla(
    .ris_mock_fit(fx$fb, fx$data, tables, attach_terms = FALSE)
  )
  expect_named(r, "gen:env")
  # Every genotype-environment cell of the data, genotype fastest.
  cells <- unique(fx$data[c("gen", "env")])
  cells <- cells[order(as.integer(cells$env), as.integer(cells$gen)), ]
  expect_identical(r[[1L]]$level, paste0(cells$gen, ":", cells$env))
})


# --- the table is matched to its specification by identity -------------- #
#
# `ID` is the index value INLA placed on each row of a table. A table is the
# field a term indexed only when those values are, element by element, the
# ones the specification says its rows carry. A count of rows agrees with a
# table in another order, and relabelling that table by position attaches
# every estimate to another level (review finding P2.1).

test_that("a table in another order at the same length is not relabelled by position", {
  fx <- .ris_fixtures()$at_simple
  n_gen <- nlevels(fx$data$gen)
  tab <- .ris_table(n_gen)
  tab <- tab[rev(seq_len(n_gen)), ]
  rownames(tab) <- NULL
  # The row count is the field's, so a check on it sees nothing.
  expect_identical(nrow(tab), n_gen)
  expect_warning(
    r <- .fb_random_inla(
      .ris_mock_fit(fx$fb, fx$data, list(gen_at_env_E1 = tab))
    ),
    "gen (env = E1)",
    fixed = TRUE
  )
  # INLA's own index values are kept, each beside its own estimate.
  expect_identical(r[[1L]]$level, as.character(rev(seq_len(n_gen))))
  expect_identical(r[[1L]]$estimate, tab$mean)
  expect_false(any(r[[1L]]$level %in% levels(fx$data$gen)))
})

test_that("every specification relabels a table carrying its index values and no permutation of them", {
  fixtures <- .ris_fixtures()
  checked <- 0L
  for (nm in names(fixtures)) {
    fx <- fixtures[[nm]]
    specs <- .fb_inla_index_specs(
      fx$fb,
      data = fx$data,
      known_matrices = fx$known_matrices
    )
    for (idx in names(specs)) {
      s <- specs[[idx]]
      if (is.null(s$levels)) {
        next
      }
      info <- paste(nm, idx)
      expect_identical(length(s$ids), as.integer(s$n_field), info = info)

      tab <- .ris_table(s$n_field)
      tab$ID <- s$ids
      expect_no_warning(
        r <- .fb_random_inla(.ris_mock_fit(
          fx$fb, fx$data, stats::setNames(list(tab), idx), fx$known_matrices
        ))
      )
      expect_identical(r[[s$label]]$level, s$levels, info = info)
      expect_identical(r[[s$label]]$estimate, tab$mean[s$rows], info = info)

      permuted <- tab
      permuted$ID <- rev(s$ids)
      expect_false(identical(permuted$ID, tab$ID), info = info)
      expect_warning(
        r <- .fb_random_inla(.ris_mock_fit(
          fx$fb, fx$data, stats::setNames(list(permuted), idx),
          fx$known_matrices
        )),
        s$label,
        fixed = TRUE
      )
      expect_identical(
        r[[s$label]]$level,
        as.character(permuted$ID),
        info = info
      )
      checked <- checked + 1L
    }
  }
  # Every generated index of every fixture with levels, none skipped.
  expect_identical(checked, 23L)
})

test_that("every specification's index values are the IDs a fitted INLA table carries", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    # This is an index-construction contract. Posterior spatial-collapse
    # diagnostics are exercised in test-spatial-field-collapse.R.
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  fixtures <- .ris_fixtures()
  checked <- 0L
  for (nm in names(fixtures)) {
    fx <- fixtures[[nm]]
    fit <- suppressMessages(emit_inla(
      fx$fb,
      data = fx$data,
      known_matrices = fx$known_matrices,
      verbose = FALSE
    ))
    specs <- .fb_inla_index_specs(
      fx$fb,
      data = fit$data,
      level_labels = fit$level_labels,
      known_matrices = fx$known_matrices
    )
    for (idx in intersect(names(specs), names(fit$inla$summary.random))) {
      s <- specs[[idx]]
      if (is.null(s$levels)) {
        next
      }
      expect_identical(
        as.numeric(fit$inla$summary.random[[idx]]$ID),
        as.numeric(s$ids),
        info = paste(nm, idx)
      )
      checked <- checked + 1L
    }
    expect_no_warning(.fb_random_inla(fit))
  }
  expect_identical(checked, 23L)
})
