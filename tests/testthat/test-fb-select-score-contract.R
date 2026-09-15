# =============================================================================
# fb_select(): an out-of-fold comparison exists only under one score contract,
# and at one model stage.
#
# The defects this closes. The side-by-side stanza compared the SET of fold
# labels with setequal(), so two different assignments of units to folds
# labelled 1..K passed; it received no unit-to-fold map, no row identity, no
# prediction target and no design signature, yet certified one partition and
# one target; and it placed the Bayesian reference beside the optimiser's
# SELECTED model, two different model stages, in one table.
#
# A score contract now carries the canonical unit-to-fold map, the scored
# rows, the prediction target, the response, the metric, the scoring unit,
# the candidate design, the preprocessing boundary and the model stage, and
# an arm enters a stage table only when its contract agrees with the primary
# arm's field for field. No sampler is involved.
#
# The contract identifies the values it scored, not only their names. The
# response and every design column enter through a canonical byte encoding
# with explicit types, lengths and byte order, digested with MD5, so a
# rescaled covariate or a transformed response under an unchanged name gives
# a different contract, and one input gives the same hash under every R
# version. A contract whose hash cannot be verified refuses by name.
# =============================================================================

fbsc_data <- function() {
  set.seed(31L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(12L)), each = 4L)))
  d$p1 <- stats::rnorm(nrow(d))
  d$x1 <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d))
  d
}

# The design both scorers record: the protected and candidate columns over
# all rows, before fold-local preprocessing.
fbsc_design <- function(data) {
  stats::model.matrix(~ p1 + x1, data)
}

fbsc_contract <- function(folds,
                          data = fbsc_data(),
                          stage = "reference",
                          prediction = "new_cluster",
                          metric = "nlpd",
                          design = fbsc_design(data),
                          learned_on = "analysis_rows",
                          scoring_unit = NULL,
                          response = "y",
                          offsets = character(0)) {
  .fb_select_score_contract(
    folds = folds,
    data = data,
    response = response,
    prediction = prediction,
    metric = metric,
    design = design,
    learned_on = learned_on,
    stage = stage,
    scoring_unit = scoring_unit,
    offsets = offsets
  )
}

# A per-unit score table of the shape both scorers write, one row per unit.
fbsc_clusters <- function(folds, model = "full", shift = 0) {
  key <- attr(folds, "unit_key")
  do.call(rbind, lapply(folds, function(s) {
    units <- unique(key[s$assess_rows])
    data.frame(
      model = model, fold = s$label, repeat_id = s$repeat_id,
      cluster = units, n = 4L,
      mse = 1 + shift, rmse = 1 + shift,
      elpd = -(seq_along(units) + shift),
      nlpd = seq_along(units) / 10 + shift,
      stringsAsFactors = FALSE
    )
  }))
}

fbsc_folds <- function(seed = 3L, data = fbsc_data()) {
  .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), data, seed = seed)
}

# The fixture behind the recorded hash: no random-number generator, no
# transcendental function and explicit folds, with the values the encoding
# treats specially (NA, NaN, -0, Inf, an offset, a unit label outside ASCII).
fbsc_golden <- function() {
  i <- seq_len(24L)
  d <- data.frame(
    site = factor(rep(c("A01", "B02", "Süd", "Z04"), each = 6L)),
    p1 = (i %% 5L - 2L) / 4,
    x1 = i / 3,
    o = i / 8,
    y = ((i * 37L) %% 101L) / 10
  )
  d$y[3L] <- NA_real_
  d$y[4L] <- NaN
  d$x1[6L] <- -0
  d$p1[7L] <- Inf
  d
}

fbsc_golden_contract <- function() {
  d <- fbsc_golden()
  folds <- .fb_sel_folds(
    fb_resamples(unit = ~site, fold_id = rep(c(1L, 2L, 1L, 2L), each = 6L)),
    d,
    seed = NULL
  )
  fbsc_contract(folds, data = d, offsets = "offset(o)")
}


test_that("one computation gives one hash, and the hash is label-invariant", {
  folds <- fbsc_folds()
  a <- fbsc_contract(folds)
  b <- fbsc_contract(folds)
  expect_s3_class(a, "fb_select_score_contract")
  expect_identical(a$hash, b$hash)
  expect_match(a$hash, "^[0-9a-f]{32}$")
  expect_identical(a$version, 2L)

  # The same partition under different fold labels and in a different order
  # is the same partition: the map is canonical, not a copy of the labels.
  relabelled <- folds
  for (i in seq_along(relabelled)) {
    relabelled[[i]]$label <- paste0("F", 4L - i)
    relabelled[[i]]$fold_id <- paste0("F", 4L - i)
  }
  relabelled <- structure(
    rev(unclass(relabelled)),
    class = class(folds),
    unit = attr(folds, "unit"), unit_key = attr(folds, "unit_key"),
    n_rows = attr(folds, "n_rows")
  )
  expect_identical(fbsc_contract(relabelled)$hash, a$hash)
})


test_that("two assignments with the same fold labels have different contracts", {
  # The defect, restated as data. Both partitions deal twelve sites into
  # folds labelled 1..3, so a comparison of label sets passes them as one.
  a <- fbsc_folds(seed = 3L)
  b <- fbsc_folds(seed = 4L)
  labels <- function(folds) vapply(folds, function(s) s$label, character(1L))
  expect_true(setequal(labels(a), labels(b)))
  expect_false(identical(
    lapply(a, function(s) s$assess_rows),
    lapply(b, function(s) s$assess_rows)
  ))

  ca <- fbsc_contract(a)
  cb <- fbsc_contract(b)
  expect_false(identical(ca$hash, cb$hash))
  expect_true("unit_to_fold" %in% .fb_select_contract_differences(ca, cb))

  # And the stage table refuses to place them side by side, naming the field.
  table <- .fb_select_stage_table(
    "reference",
    list(
      list(arm = "bayes", model = "full", clusters = fbsc_clusters(a),
           contract = ca),
      list(arm = "gls", model = "full", clusters = fbsc_clusters(b),
           contract = cb)
    ),
    metric = "nlpd"
  )
  expect_false(table$available)
  expect_identical(table$excluded$arm, "gls")
  expect_match(table$reason, "unit_to_fold", fixed = TRUE)
})


test_that("every contract field moves the hash", {
  d <- fbsc_data()
  folds <- fbsc_folds(data = d)
  base <- fbsc_contract(folds, data = d)
  shifted <- d
  shifted$y <- shifted$y + 1
  rescaled <- fbsc_design(d)
  rescaled[, "x1"] <- 1000 * rescaled[, "x1"] + 17

  variants <- list(
    prediction = fbsc_contract(folds, data = d, prediction = "known_cluster"),
    response = fbsc_contract(folds, data = shifted),
    metric = fbsc_contract(folds, data = d, metric = "rmse"),
    design_columns = fbsc_contract(
      folds, data = d, design = stats::model.matrix(~p1, d)
    ),
    design_values = fbsc_contract(folds, data = d, design = rescaled),
    learned_on = fbsc_contract(folds, data = d, learned_on = "all_rows"),
    stage = fbsc_contract(folds, data = d, stage = "decision"),
    scoring_unit = fbsc_contract(folds, data = d, scoring_unit = "fold")
  )
  # Dropping a column changes the names and so the values; every other
  # variant changes one field.
  expected <- list(design_columns = c("design_columns", "design_values"))
  for (field in names(variants)) {
    expect_false(identical(variants[[field]]$hash, base$hash), info = field)
    expect_identical(
      .fb_select_contract_differences(base, variants[[field]]),
      expected[[field]] %||% field,
      info = field
    )
  }
})


test_that("a transformed scoring frame with unchanged names has a different contract", {
  d <- fbsc_data()
  folds <- fbsc_folds(data = d)
  base <- fbsc_contract(folds, data = d)
  # A fixed permutation, so the check does not depend on a random draw.
  permutation <- c(seq.int(2L, nrow(d)), 1L)

  # A design column transformed under its own name.
  frames <- list(
    rescaled = within(d, x1 <- 1000 * x1 + 17),
    shifted = within(d, x1 <- x1 + 1),
    one_value = within(d, x1[5L] <- x1[5L] + 1e-8),
    permuted_within_column = within(d, x1 <- x1[permutation])
  )
  for (nm in names(frames)) {
    variant <- fbsc_contract(folds, data = frames[[nm]])
    expect_identical(
      colnames(fbsc_design(frames[[nm]])), colnames(fbsc_design(d)),
      info = nm
    )
    expect_false(identical(variant$hash, base$hash), info = nm)
    expect_identical(
      .fb_select_contract_differences(base, variant), "design_values",
      info = nm
    )
  }

  # The response transformed under its own name.
  responses <- list(
    rescaled = within(d, y <- 10 * y),
    one_value = within(d, y[5L] <- y[5L] + 1e-8),
    permuted_within_column = within(d, y <- y[permutation])
  )
  for (nm in names(responses)) {
    variant <- fbsc_contract(folds, data = responses[[nm]])
    expect_false(identical(variant$hash, base$hash), info = nm)
    expect_identical(
      .fb_select_contract_differences(base, variant), "response", info = nm
    )
  }

  # A response written as an expression is read through the fit's resolver,
  # so its values enter: log(y) over y and over 2 * y are different.
  positive <- within(d, y <- exp(y))
  doubled <- within(positive, y <- 2 * y)
  log_a <- fbsc_contract(folds, data = positive, response = "log(y)")
  log_b <- fbsc_contract(folds, data = doubled, response = "log(y)")
  expect_false(identical(log_a$hash, log_b$hash))
  expect_identical(.fb_select_contract_differences(log_a, log_b), "response")

  # An offset is a transformation of the scored response: its values enter.
  with_offset <- within(d, o <- p1 / 2)
  off_a <- fbsc_contract(folds, data = with_offset, offsets = "offset(o)")
  off_b <- fbsc_contract(
    folds, data = within(with_offset, o <- o + 1), offsets = "offset(o)"
  )
  expect_false(identical(off_a$hash, off_b$hash))
  expect_identical(.fb_select_contract_differences(off_a, off_b), "response")
  expect_identical(
    .fb_select_contract_differences(base, off_a), "response"
  )

  # A response that neither names a column nor evaluates refuses by name.
  expect_error(
    fbsc_contract(folds, data = d, response = "log(absent)"),
    class = "flexybayes_refusal_selection_fold_variable_not_in_data"
  )
})


test_that("identical and rebuilt frames have one contract", {
  d <- fbsc_data()
  folds <- fbsc_folds(data = d)
  a <- fbsc_contract(folds, data = d)
  expect_identical(fbsc_contract(folds, data = d)$hash, a$hash)

  # A fresh object from the same values: other column order, other row names,
  # the factor rebuilt from its labels.
  rebuilt <- data.frame(
    y = as.numeric(unlist(d$y)),
    x1 = c(d$x1),
    p1 = d$p1 + 0,
    site = factor(as.character(d$site), levels = levels(d$site)),
    row.names = paste0("r", seq_len(nrow(d)))
  )
  expect_false(identical(names(rebuilt), names(d)))
  b <- fbsc_contract(folds, data = rebuilt)
  expect_identical(b$hash, a$hash)
  expect_identical(.fb_select_contract_differences(a, b), character(0))

  # A design with its columns in another order is the same design.
  X <- fbsc_design(d)
  expect_identical(
    fbsc_contract(folds, data = d, design = X[, c("x1", "(Intercept)", "p1")])$hash,
    a$hash
  )
})


test_that("a stage table excludes an arm scored on a transformed frame with unchanged names", {
  d <- fbsc_data()
  folds <- fbsc_folds(data = d)
  # Shaped as the out-of-fold objects both scorers return.
  arm <- function(X) {
    list(
      folds = folds, clusters = fbsc_clusters(folds),
      prediction = "new_cluster",
      preprocessing = list(learned_on = "analysis_rows"),
      substrates = list(F1 = list(
        x_protected = X[, "(Intercept)", drop = FALSE],
        x_candidate = X[, c("p1", "x1"), drop = FALSE]
      )),
      design = X
    )
  }
  X <- fbsc_design(d)
  X_transformed <- X
  X_transformed[, "x1"] <- 1000 * X_transformed[, "x1"] + 17
  tables <- function(gls_design) {
    .fb_select_stage_tables(
      base = list(data = d, fb = list(response = "y")),
      axes = list(workflow = "bayes", prior = "none", decision = "full"),
      cv = arm(X),
      benchmark = list(arms = list(gls = list()), cv = arm(gls_design)),
      metric = "nlpd"
    )
  }

  same <- tables(X)$reference
  expect_true(same$available)
  expect_setequal(same$arms$arm, c("bayes", "gls"))

  transformed <- tables(X_transformed)$reference
  expect_false(transformed$available)
  expect_identical(transformed$excluded$arm, "gls")
  expect_match(transformed$excluded$reason, "design_values", fixed = TRUE)
})


test_that("a stage table is written only when every arm shares the primary's contract", {
  folds <- fbsc_folds()
  contract <- fbsc_contract(folds)
  entries <- list(
    list(arm = "bayes", model = "full", clusters = fbsc_clusters(folds),
         contract = contract),
    list(arm = "gls", model = "full",
         clusters = fbsc_clusters(folds, shift = 0.2), contract = contract)
  )
  table <- .fb_select_stage_table("reference", entries, metric = "nlpd")
  expect_true(table$available)
  expect_setequal(table$arms$arm, c("bayes", "gls"))
  expect_true(all(table$arms$stage == "reference"))
  expect_identical(unique(table$arms$contract_hash), contract$hash)
  expect_true(all(table$arms$n_units == 12L))
  # The scores are re-aggregated from the per-unit tables by one function.
  expect_equal(
    table$arms$score[table$arms$arm == "gls"],
    .fb_sel_aggregate_bayes_scores(
      entries[[2L]]$clusters, metric = "nlpd", unit = "cluster"
    )$score
  )

  # A stage mismatch is a contract mismatch: the reference posterior is not
  # placed beside a selected model.
  selected <- entries
  selected[[2L]]$model <- "benchmark_7"
  selected[[2L]]$clusters <- fbsc_clusters(folds, model = "benchmark_7")
  selected[[2L]]$contract <- fbsc_contract(folds, stage = "decision")
  mixed <- .fb_select_stage_table("reference", selected, metric = "nlpd")
  expect_false(mixed$available)
  expect_match(mixed$reason, "stage", fixed = TRUE)

  # Equal unit COUNTS are not the same units. One arm scores a unit its
  # contract does not name in place of one it does.
  swapped <- entries
  swapped[[2L]]$clusters$cluster[[1L]] <- "S99"
  expect_identical(
    length(unique(swapped[[2L]]$clusters$cluster)), 12L
  )
  swapped_table <- .fb_select_stage_table("reference", swapped, metric = "nlpd")
  expect_false(swapped_table$available)
  expect_match(swapped_table$excluded$reason, "units", fixed = TRUE)

  # A table with no primary score, or with one arm, says why and holds no
  # rows.
  alone <- .fb_select_stage_table("decision", entries[1L], metric = "nlpd")
  expect_false(alone$available)
  expect_null(alone$arms)
  named <- .fb_select_stage_table(
    "decision", entries, metric = "nlpd",
    unavailable = "a posterior threshold defines no predictive distribution"
  )
  expect_false(named$available)
  expect_match(named$reason, "posterior threshold", fixed = TRUE)
})


test_that("a contract whose hash cannot be verified refuses by name", {
  folds <- fbsc_folds()
  contract <- fbsc_contract(folds)
  entries <- function(primary, second) {
    list(
      list(arm = "bayes", model = "full", clusters = fbsc_clusters(folds),
           contract = primary),
      list(arm = "gls", model = "full",
           clusters = fbsc_clusters(folds, shift = 0.2), contract = second)
    )
  }
  refusal <- function(primary, second) {
    tryCatch(
      .fb_select_stage_table(
        "reference", entries(primary, second), metric = "nlpd"
      ),
      flexybayes_refusal_selection_score_contract_unverifiable = function(e) e
    )
  }
  expect_true(.fb_select_stage_table(
    "reference", entries(contract, contract), metric = "nlpd"
  )$available)

  # A contract written before this release: no version field, and the
  # polynomial checksum over names in place of the MD5 hash.
  legacy <- contract
  legacy$version <- NULL
  legacy$hash <- "620415bc-7d39e25e"
  err <- refusal(contract, legacy)
  expect_s3_class(
    err, "flexybayes_refusal_selection_score_contract_unverifiable"
  )
  expect_match(conditionMessage(err), "`gls` arm", fixed = TRUE)
  expect_match(conditionMessage(err), "no version field", fixed = TRUE)
  # The primary arm is verified too.
  err <- refusal(legacy, contract)
  expect_s3_class(
    err, "flexybayes_refusal_selection_score_contract_unverifiable"
  )
  expect_match(conditionMessage(err), "`bayes` arm", fixed = TRUE)

  # A version field does not make a polynomial checksum verifiable.
  polynomial <- contract
  polynomial$hash <- "620415bc-7d39e25e"
  expect_match(
    conditionMessage(refusal(contract, polynomial)), "not an MD5", fixed = TRUE
  )

  # Another version is not read as this one.
  later <- contract
  later$version <- 3L
  expect_match(
    conditionMessage(refusal(contract, later)), "version 3", fixed = TRUE
  )

  # A field changed after the hash was computed.
  tampered <- contract
  tampered$design_values <- contract$response$values
  expect_match(
    conditionMessage(refusal(contract, tampered)), "does not match its fields",
    fixed = TRUE
  )

  # An arm whose source recorded no design carries no contract. Without a
  # score table it is excluded for that reason; with one it refuses.
  no_table <- entries(contract, NULL)
  no_table[[2L]]$clusters <- NULL
  excluded <- .fb_select_stage_table("reference", no_table, metric = "nlpd")
  expect_false(excluded$available)
  expect_match(
    excluded$excluded$reason, "no out-of-fold score table", fixed = TRUE
  )
  expect_match(
    conditionMessage(refusal(contract, NULL)), "not a score contract",
    fixed = TRUE
  )
})


test_that("the canonical encoding writes explicit types, lengths and little-endian bytes", {
  con <- rawConnection(raw(0), open = "wb")
  on.exit(close(con))
  .fb_select_write_field(con, "x", c(1, NA, NaN, -0))
  .fb_select_write_field(con, "s", c("é", NA))
  .fb_select_write_field(con, "f", factor(c("b", NA), levels = c("a", "b")))
  .fb_select_write_field(con, "l", c(TRUE, NA))
  bytes <- paste(as.character(rawConnectionValue(con)), collapse = "")

  # Written from IEEE 754 and the layout, not from the encoder.
  expected <- paste0(
    # the name "x": S, one element of one byte, "x"
    "53", "01000000", "01000000", "78",
    # c(1, NA, NaN, -0): D, four values, NA at 2, NaN at 3, then 1.0 and
    # three +0 (NA, NaN and -0 are not written as bit patterns)
    "44", "04000000",
    "49", "01000000", "02000000",
    "49", "01000000", "03000000",
    "000000000000f03f", "0000000000000000", "0000000000000000",
    "0000000000000000",
    # the name "s", then c("e acute", NA): lengths 2 and -1, UTF-8 c3 a9
    "53", "01000000", "01000000", "73",
    "53", "02000000", "02000000", "ffffffff", "c3a9",
    # the name "f", then the factor: levels "a", "b"; codes 2 and NA
    "53", "01000000", "01000000", "66",
    "46", "53", "02000000", "01000000", "01000000", "6162",
    "49", "02000000", "02000000", "00000080",
    # the name "l", then c(TRUE, NA)
    "53", "01000000", "01000000", "6c",
    "4c", "02000000", "01000000", "00000080"
  )
  expect_identical(bytes, expected)
})


test_that("one fixture's contract is the same bytes under R 4.5.2 and R 4.6.1", {
  contract <- fbsc_golden_contract()
  # Recorded under R 4.5.2 (aarch64-apple-darwin20) and R 4.6.1
  # (aarch64-apple-darwin23), identical under both. A change in the encoding,
  # or a byte that depends on the R version, moves these strings.
  expect_identical(contract$response$values, "53517c7013eaf7027dce3fce3ec066ba")
  expect_identical(contract$design_values, "3ce2587bbc9a7c5d99dd8e672fecccf6")
  expect_identical(contract$hash, "25f49f50bdd7fc2e5bd26ecba37d89bd")
})


test_that("a digest is an MD5 of the encoding in a temporary file that is removed", {
  digest_of <- function(value) {
    .fb_select_digest("t", function(con) .fb_select_write_field(con, "v", value))
  }
  before <- list.files(tempdir(), pattern = "^fb_score_contract_")
  one <- digest_of(1)
  expect_match(one, "^[0-9a-f]{32}$")
  expect_identical(digest_of(1), one)
  expect_false(identical(digest_of(1 + 2^-52), one))
  expect_false(identical(digest_of(NA_real_), digest_of(NaN)))
  expect_identical(digest_of(-0), digest_of(0))
  expect_identical(list.files(tempdir(), pattern = "^fb_score_contract_"), before)

  # The file route computes MD5: the RFC 1321 test vector for "abc".
  path <- tempfile()
  on.exit(unlink(path))
  writeBin(charToRaw("abc"), path)
  expect_identical(
    unname(tools::md5sum(path)), "900150983cd24fb0d6963f7d28e17f72"
  )
})
