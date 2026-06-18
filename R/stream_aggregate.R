# stream_aggregate.R -- Out-of-core streaming exact aggregation.
#
# Fits a cell-aggregatable mixed model on a dataset that is too large to
# hold in memory at once. The data is read in row-range chunks from an
# on-disk source (an `.fst` file, a delimited file, or a user-supplied
# chunk generator); per design-cell sufficient statistics are
# accumulated additively across chunks and never expanded back to N
# rows. The aggregated likelihood is algebraically identical to the
# per-row likelihood (gaussian: sum of squares; binomial: sum of
# successes / trials; poisson: sum of counts / exposure), so the fitted
# posterior matches a full-data fit -- this is exact compression, not an
# approximation.
#
# Peak memory is bounded by one chunk plus the K-cell accumulator, not
# by N. The binding constraint is the cell count K (distinct combinations
# of the factor design), not the row count: factorial designs with heavy
# replication compress to a tiny K and fit in seconds regardless of N;
# designs keyed on a continuous covariate do not compress and are
# refused before any data is read.
#
# References. The sufficient-statistic aggregation contract is the
# gaussian aggregation contract extended to the exponential-family
# count cases; the in-memory counterpart is `.fb_aggregate_gaussian()`.

# ---------------------------------------------------------------- #
# User-facing entry                                                 #
# ---------------------------------------------------------------- #

#' Fit a mixed model to an out-of-core dataset by streaming aggregation
#'
#' Fit a cell-aggregatable mixed model to a dataset held on disk, reading
#' it in row-range chunks and accumulating exact per-cell sufficient
#' statistics rather than materialising the full table in memory.
#'
#' @description
#' `flexybayes_stream()` is the large-data entry point: it reads the data
#' a chunk at a time from `source`, accumulates additive sufficient
#' statistics per design cell, and fits the resulting aggregated model
#' through the same greta or INLA emit path as [flexybayes()]. The
#' aggregated likelihood is algebraically identical to the per-row
#' likelihood, so the posterior is the full-data posterior, not an
#' approximation.
#'
#' @details
#' The method scales by cell count `K` (the number of distinct factor
#' design combinations), not by row count `N`. A replicated factorial
#' design collapses billions of rows to a few thousand cells and fits in
#' seconds; a design that keys on a continuous covariate does not
#' compress and is refused before any chunk is read. Continuous fixed
#' effects, random slopes, structured covariance, and smooth terms break
#' the cell-constant linear-predictor property and are likewise refused
#' -- pass the data to [flexybayes()] per row for those models.
#'
#' Supported families are `"gaussian"` (identity link), `"binomial"`,
#' and `"poisson"`. For binomial data supply `trials` to name the column
#' of trial counts (a 0/1 response is treated as Bernoulli with one trial
#' per row). For poisson data supply `exposure` to name an exposure /
#' offset column.
#'
#' @param fixed A two-sided formula `response ~ fixed_terms`. The fixed
#'   terms must be factors or factor interactions; a continuous term
#'   forces one cell per row and is refused.
#' @param random A one-sided formula of random-intercept grouping factors,
#'   ASReml style (for example `~ geno`), or `NULL` for no random terms.
#' @param source The data source. One of: a length-1 character path to an
#'   `.fst` file; a length-1 character path to a delimited file readable
#'   by [data.table::fread()]; an in-memory `data.frame` / `data.table`
#'   (chunked internally, mainly for testing); or a function of one
#'   argument `i` returning the `i`-th chunk as a `data.frame` and `NULL`
#'   once the chunks are exhausted.
#' @param family A length-1 character family: `"gaussian"`, `"binomial"`,
#'   or `"poisson"`.
#' @param trials For `family = "binomial"`, the name of the trials-count
#'   column, or `NULL` for Bernoulli (one trial per row).
#' @param exposure For `family = "poisson"`, the name of the exposure /
#'   offset column, or `NULL` for unit exposure.
#' @param backend The estimation backend, `"inla"` (default) or
#'   `"greta"`.
#' @param chunk_rows The number of rows to read per chunk. Larger chunks
#'   read faster but use more peak memory; the default 5e6 keeps peak
#'   memory modest while amortising read overhead.
#' @param prior An optional [fb_prior()] object. When `NULL` the backend
#'   default priors are used.
#' @param fit When `TRUE` (default) the aggregated model is fitted and a
#'   `<flexybayes>` object is returned. When `FALSE` the function returns
#'   the `<fb_aggregated>` sufficient-statistics object without fitting --
#'   useful for inspecting the compression a design will achieve.
#' @param verbose When `TRUE`, report streaming progress and the
#'   compression achieved.
#' @param ... Further arguments passed to the aggregated emit (for
#'   example `n_samples`, `warmup`, `chains` for the greta backend).
#'
#' @returns When `fit = TRUE`, a `<flexybayes>` fit object carrying
#'   `extras$aggregation_meta$streamed == TRUE`. When `fit = FALSE`, an
#'   `<fb_aggregated>` object.
#'
#' @seealso [flexybayes()] for the in-memory entry point;
#'   [fb_prior()] for prior specification.
#'
#' @examplesIf interactive() && requireNamespace("fst", quietly = TRUE)
#' set.seed(1L)
#' n <- 1e6
#' df <- data.frame(
#'   env = factor(sample(letters[1:6], n, replace = TRUE)),
#'   geno = factor(sample(1:50, n, replace = TRUE)),
#'   y = rnorm(n)
#' )
#' path <- tempfile(fileext = ".fst")
#' fst::write_fst(df, path)
#' fit <- flexybayes_stream(y ~ env, random = ~ geno, source = path,
#'                          backend = "inla")
#' summary(fit)
#'
#' @author Max Moldovan, \email{max.moldovan@@adelaide.edu.au}
#' @export
flexybayes_stream <- function(
  fixed,
  random = NULL,
  source,
  family = "gaussian",
  trials = NULL,
  exposure = NULL,
  backend = c("inla", "greta"),
  chunk_rows = 5e6,
  prior = NULL,
  fit = TRUE,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  the_call <- match.call()
  .check_stream_family(family)
  .check_stream_chunk_rows(chunk_rows)

  # IR is built from a head sample so term classification and factor
  # contrasts match the in-memory path exactly. An `.fst` source carries
  # the full factor-level dictionary in its metadata, so any row range
  # reads back with the complete level set and cell keys stay aligned
  # across chunks.
  src <- .fb_stream_source(source, chunk_rows = as.integer(chunk_rows))
  head <- .fb_stream_head(src, fixed, random, trials, exposure)
  fb <- fb_from_asreml(
    fixed = fixed,
    random = random,
    data = head,
    family = family,
    prior = prior
  )

  agg <- .fb_stream_aggregate(
    fb = fb,
    src = src,
    trials = trials,
    exposure = exposure,
    verbose = verbose
  )

  if (!isTRUE(fit)) {
    return(agg)
  }

  .fb_stream_emit(
    fb = fb,
    agg = agg,
    head = head,
    backend = backend,
    fixed = fixed,
    random = random,
    family = family,
    the_call = the_call,
    verbose = verbose,
    ...
  )
}


# ---------------------------------------------------------------- #
# Validation helpers                                                #
# ---------------------------------------------------------------- #

#' @noRd
#' @keywords internal
.check_stream_family <- function(family) {
  ok <- c("gaussian", "binomial", "poisson")
  if (!is.character(family) || length(family) != 1L || !family %in% ok) {
    stop(
      "`family` must be one of ",
      paste(dQuote(ok), collapse = ", "),
      "; got: ",
      deparse(family),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
#' @keywords internal
.check_stream_chunk_rows <- function(chunk_rows) {
  if (
    !is.numeric(chunk_rows) ||
      length(chunk_rows) != 1L ||
      !is.finite(chunk_rows) ||
      chunk_rows < 1
  ) {
    stop(
      "`chunk_rows` must be a positive numeric scalar; got: ",
      deparse(chunk_rows),
      call. = FALSE
    )
  }
  if (chunk_rows >= 2^31) {
    stop(
      "`chunk_rows` must be below 2^31 (a single chunk is read into ",
      "memory at once).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------- #
# Source abstraction                                                #
# ---------------------------------------------------------------- #
# Normalises the four accepted source shapes to a single reader
# protocol. Returns a list with:
#   $kind       "fst" | "delim" | "memory" | "generator"
#   $n_rows     integer total rows, or NA_integer_ for a generator
#   $n_chunks   integer chunk count, or NA_integer_ for a generator
#   $columns    available column names, or NULL for a generator
#   $read(i, cols)  read chunk i (1-based) as a data.table, selecting
#                   `cols` when the source supports column projection;
#                   returns NULL once chunks are exhausted.
#   $chunk_rows the configured chunk size.

#' @noRd
#' @keywords internal
.fb_stream_source <- function(source, chunk_rows) {
  if (is.function(source)) {
    return(.fb_stream_source_generator(source, chunk_rows))
  }
  if (inherits(source, "data.frame")) {
    return(.fb_stream_source_memory(source, chunk_rows))
  }
  if (is.character(source)) {
    miss <- source[!file.exists(source)]
    if (length(miss)) {
      stop(
        "`source` file(s) do not exist: ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    is_fst <- grepl("\\.fst$", source, ignore.case = TRUE)
    if (length(source) > 1L) {
      if (!all(is_fst)) {
        stop(
          "A multi-file `source` must be a vector of .fst shard ",
          "paths (partitioned dataset).",
          call. = FALSE
        )
      }
      return(.fb_stream_source_fst_multi(source, chunk_rows))
    }
    if (is_fst) {
      return(.fb_stream_source_fst(source, chunk_rows))
    }
    return(.fb_stream_source_delim(source, chunk_rows))
  }
  stop(
    "`source` must be an .fst path (or vector of .fst shard paths), ",
    "a delimited-file path, a data.frame / data.table, or a ",
    "chunk-generator function.",
    call. = FALSE
  )
}

#' @noRd
#' @keywords internal
.fb_stream_source_fst <- function(path, chunk_rows) {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop(
      "Reading an .fst source requires the 'fst' package. Install ",
      "with install.packages('fst').",
      call. = FALSE
    )
  }
  proxy <- fst::fst(path)
  n_rows <- nrow(proxy)
  cols <- names(proxy)
  n_chunks <- as.integer(ceiling(n_rows / chunk_rows))
  read <- function(i, cols_sel = NULL) {
    from <- (i - 1L) * chunk_rows + 1L
    if (from > n_rows) {
      return(NULL)
    }
    to <- min(i * chunk_rows, n_rows)
    sel <- if (is.null(cols_sel)) cols else intersect(cols_sel, cols)
    fst::read_fst(path, columns = sel, from = from, to = to)
  }
  head <- function(n) {
    to <- as.integer(min(n, n_rows))
    fst::read_fst(path, from = 1L, to = to)
  }
  list(
    kind = "fst",
    n_rows = as.integer(n_rows),
    n_chunks = n_chunks,
    columns = cols,
    read = read,
    head = head,
    chunk_rows = chunk_rows
  )
}

#' @noRd
#' @keywords internal
.fb_stream_source_fst_multi <- function(paths, chunk_rows) {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Reading .fst shards requires the 'fst' package.", call. = FALSE)
  }
  # A partitioned dataset: each shard is read in row-range chunks, and a
  # single chunk plan walks (shard, from, to) triples across all shards
  # in order. Peak memory stays bounded by one chunk regardless of the
  # total row count across shards.
  per_file <- vapply(paths, function(p) nrow(fst::fst(p)), numeric(1L))
  cols <- names(fst::fst(paths[[1L]]))
  plan <- list()
  for (f in seq_along(paths)) {
    nf <- per_file[[f]]
    if (nf == 0) {
      next
    }
    starts <- seq.int(1L, nf, by = chunk_rows)
    for (s in starts) {
      plan[[length(plan) + 1L]] <-
        list(path = paths[[f]], from = s, to = min(s + chunk_rows - 1L, nf))
    }
  }
  read <- function(i, cols_sel = NULL) {
    if (i > length(plan)) {
      return(NULL)
    }
    pl <- plan[[i]]
    sel <- if (is.null(cols_sel)) cols else intersect(cols_sel, cols)
    fst::read_fst(pl$path, columns = sel, from = pl$from, to = pl$to)
  }
  head <- function(n) {
    fst::read_fst(
      paths[[1L]],
      from = 1L,
      to = as.integer(min(n, per_file[[1L]]))
    )
  }
  # Total rows are carried as a double, not an integer: a partitioned
  # dataset can exceed the 2^31 integer ceiling (a 5-billion-row dataset
  # is 20 shards of 250M). Doubles are exact for integer counts to 2^53
  # (~9e15 rows), and the per-cell accumulators are doubles for the same
  # reason, so the streaming total never overflows.
  list(
    kind = "fst_multi",
    n_rows = sum(per_file),
    n_chunks = length(plan),
    columns = cols,
    read = read,
    head = head,
    chunk_rows = chunk_rows
  )
}

#' @noRd
#' @keywords internal
.fb_stream_source_delim <- function(path, chunk_rows) {
  # Row count via a header-aware line count; data.table::fread reads a
  # row range with `skip` (data rows) + `nrows`. The header is re-read
  # per chunk via `fread(nrows = 0)` to recover column names.
  header <- names(data.table::fread(path, nrows = 0L))
  n_lines <- .fb_stream_count_lines(path)
  n_rows <- max(n_lines - 1L, 0L)
  n_chunks <- as.integer(ceiling(n_rows / chunk_rows))
  read <- function(i, cols_sel = NULL) {
    skip <- (i - 1L) * chunk_rows
    if (skip >= n_rows) {
      return(NULL)
    }
    nr <- min(chunk_rows, n_rows - skip)
    dt <- data.table::fread(
      path,
      skip = skip + 1L,
      nrows = nr,
      header = FALSE,
      col.names = header,
      select = if (is.null(cols_sel)) NULL else intersect(cols_sel, header)
    )
    as.data.frame(dt)
  }
  head <- function(n) {
    as.data.frame(data.table::fread(path, nrows = as.integer(n)))
  }
  list(
    kind = "delim",
    n_rows = as.integer(n_rows),
    n_chunks = n_chunks,
    columns = header,
    read = read,
    head = head,
    chunk_rows = chunk_rows
  )
}

#' @noRd
#' @keywords internal
.fb_stream_source_memory <- function(data, chunk_rows) {
  df <- as.data.frame(data)
  n_rows <- nrow(df)
  n_chunks <- as.integer(max(ceiling(n_rows / chunk_rows), 1L))
  read <- function(i, cols_sel = NULL) {
    from <- (i - 1L) * chunk_rows + 1L
    if (from > n_rows) {
      return(NULL)
    }
    to <- min(i * chunk_rows, n_rows)
    sel <- if (is.null(cols_sel)) names(df) else intersect(cols_sel, names(df))
    df[from:to, sel, drop = FALSE]
  }
  head <- function(n) df[seq_len(as.integer(min(n, n_rows))), , drop = FALSE]
  list(
    kind = "memory",
    n_rows = as.integer(n_rows),
    n_chunks = n_chunks,
    columns = names(df),
    read = read,
    head = head,
    chunk_rows = chunk_rows
  )
}

#' @noRd
#' @keywords internal
.fb_stream_source_generator <- function(gen, chunk_rows) {
  read <- function(i, cols_sel = NULL) {
    chunk <- gen(i)
    if (is.null(chunk)) {
      return(NULL)
    }
    df <- as.data.frame(chunk)
    if (!is.null(cols_sel)) {
      df <- df[, intersect(cols_sel, names(df)), drop = FALSE]
    }
    df
  }
  head <- function(n) {
    chunk <- gen(1L)
    if (is.null(chunk)) {
      return(NULL)
    }
    df <- as.data.frame(chunk)
    df[seq_len(min(nrow(df), as.integer(n))), , drop = FALSE]
  }
  list(
    kind = "generator",
    n_rows = NA_integer_,
    n_chunks = NA_integer_,
    columns = NULL,
    read = read,
    head = head,
    chunk_rows = chunk_rows
  )
}

#' @noRd
#' @keywords internal
.fb_stream_count_lines <- function(path) {
  # Stream the file in binary blocks counting newlines; avoids reading
  # the whole file into memory just to size it.
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  block <- 1048576L
  n <- 0L
  repeat {
    buf <- readBin(con, what = "raw", n = block)
    if (length(buf) == 0L) {
      break
    }
    n <- n + sum(buf == as.raw(10L))
  }
  # A file with no trailing newline still has a final line.
  n + 1L
}


# ---------------------------------------------------------------- #
# Head sample (IR construction + level dictionary)                  #
# ---------------------------------------------------------------- #

#' @noRd
#' @keywords internal
.fb_stream_head <- function(
  src,
  fixed,
  random,
  trials,
  exposure,
  head_n = 1e5L
) {
  # The head sample drives IR construction (term classification, factor
  # contrasts). It is read directly by the source so column types and
  # .fst factor-level dictionaries come from the data itself.
  head <- src$head(head_n)
  if (is.null(head) || nrow(head) == 0L) {
    stop(
      "The source returned no rows for the head sample; no data to ",
      "build the model from.",
      call. = FALSE
    )
  }
  model_cols <- .fb_stream_model_cols(fixed, random, trials, exposure)
  missing <- setdiff(model_cols, names(head))
  if (length(missing)) {
    stop(
      "Model columns absent from the data source: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  head[]
}

#' @noRd
#' @keywords internal
.fb_stream_model_cols <- function(fixed, random, trials, exposure) {
  cols <- all.vars(fixed)
  if (!is.null(random)) {
    cols <- c(cols, all.vars(random))
  }
  if (!is.null(trials)) {
    cols <- c(cols, trials)
  }
  if (!is.null(exposure)) {
    cols <- c(cols, exposure)
  }
  unique(cols)
}


# ---------------------------------------------------------------- #
# Streaming accumulation                                            #
# ---------------------------------------------------------------- #
# Reads every chunk, computes per-cell partial sufficient statistics,
# and folds them into a running K-row accumulator. The fold is a plain
# additive sum keyed by the cell columns -- exact and order-independent
# up to floating-point summation. Counts and sums are carried as doubles
# so totals stay exact past the 2^31 integer ceiling (doubles are exact
# integers to 2^53, roughly 9e15).

#' @noRd
#' @keywords internal
.fb_stream_aggregate <- function(fb, src, trials, exposure, verbose) {
  .assert_aggregate_in_scope_stream(fb)

  key_cols <- .fb_stream_key_cols(fb)
  if (!length(key_cols)) {
    stop(
      "Streaming aggregation needs at least one factor cell-key term ",
      "(a fixed factor or a random grouping factor). An intercept-only ",
      "model does not compress.",
      call. = FALSE
    )
  }

  response <- fb$response
  family <- fb$family
  read_cols <- unique(c(response, key_cols, trials, exposure))

  acc <- NULL
  n_total <- 0
  i <- 1L
  repeat {
    chunk <- src$read(i, read_cols)
    if (is.null(chunk)) {
      break
    }
    if (nrow(chunk) == 0L) {
      i <- i + 1L
      next
    }

    partial <- .fb_stream_chunk_partial(
      chunk,
      key_cols,
      response,
      family,
      trials,
      exposure
    )
    acc <- .fb_stream_fold(acc, partial, key_cols)
    n_total <- n_total + nrow(chunk)

    if (isTRUE(verbose)) {
      .fb_stream_report_progress(i, src$n_chunks, n_total, nrow(acc))
    }
    i <- i + 1L
  }
  if (is.null(acc) || nrow(acc) == 0L) {
    stop(
      "No rows were read from the source; nothing to aggregate.",
      call. = FALSE
    )
  }

  .fb_stream_finalise(fb, acc, key_cols, n_total, family)
}

#' @noRd
#' @keywords internal
.fb_stream_key_cols <- function(fb) {
  fixed_factor <- character(0L)
  for (t in fb$fixed_terms) {
    ttype <- t$type %||% "expression"
    if (ttype %in% c("factor", "categorical")) {
      fixed_factor <- c(fixed_factor, as.character(t$var))
    } else if (identical(ttype, "factor_interaction")) {
      fixed_factor <- c(fixed_factor, as.character(t$vars))
    }
  }
  random_grp <- character(0L)
  for (t in fb$random_terms) {
    if (identical(t$type %||% "simple", "simple")) {
      random_grp <- c(random_grp, as.character(t$var))
    }
  }
  unique(c(fixed_factor, random_grp))
}

#' @noRd
#' @keywords internal
.fb_stream_chunk_partial <- function(
  chunk,
  key_cols,
  response,
  family,
  trials,
  exposure
) {
  # The value columns are named `.y` / `.m` / `.e` -- symbols that are
  # NOT local variables -- so the substitute() idiom below inlines the
  # data.table value (`kdt`) and the `by` columns (`key_cols`) while
  # leaving `.N`, `.y`, `.m`, `.e` symbolic for data.table to resolve as
  # columns. Counts and sums are doubles for overflow safety past 2^31.
  kdt <- data.table::as.data.table(chunk[key_cols])
  data.table::set(kdt, j = ".y", value = as.numeric(chunk[[response]]))

  if (identical(family, "gaussian")) {
    expr <- substitute(
      kdt[,
        list(n_k = as.double(.N), S1_k = sum(.y), S2_k = sum(.y * .y)),
        by = key_cols
      ],
      list(kdt = kdt, key_cols = key_cols)
    )
  } else if (identical(family, "binomial")) {
    mv <- if (is.null(trials)) {
      rep(1, nrow(kdt))
    } else {
      as.numeric(chunk[[trials]])
    }
    data.table::set(kdt, j = ".m", value = mv)
    expr <- substitute(
      kdt[,
        list(n_k = as.double(.N), succ_k = sum(.y), trials_k = sum(.m)),
        by = key_cols
      ],
      list(kdt = kdt, key_cols = key_cols)
    )
  } else {
    ev <- if (is.null(exposure)) {
      rep(1, nrow(kdt))
    } else {
      as.numeric(chunk[[exposure]])
    }
    data.table::set(kdt, j = ".e", value = ev)
    expr <- substitute(
      kdt[,
        list(n_k = as.double(.N), count_k = sum(.y), expo_k = sum(.e)),
        by = key_cols
      ],
      list(kdt = kdt, key_cols = key_cols)
    )
  }
  eval(expr, envir = asNamespace("data.table"))
}

#' @noRd
#' @keywords internal
.fb_stream_fold <- function(acc, partial, key_cols) {
  if (is.null(acc)) {
    return(partial)
  }
  stat_cols <- setdiff(names(partial), key_cols)
  combined <- data.table::rbindlist(list(acc, partial))

  # Sum each statistic column by the cell key. The j-list call references
  # the stat columns by their (non-local) symbol names, so substitute()
  # inlines only `combined` and `key_cols`.
  j_list <- as.call(c(
    quote(list),
    stats::setNames(
      lapply(stat_cols, function(cn) bquote(sum(.(as.name(cn))))),
      stat_cols
    )
  ))
  expr <- substitute(
    combined[, JLIST, by = key_cols],
    list(combined = combined, key_cols = key_cols, JLIST = j_list)
  )
  eval(expr, envir = asNamespace("data.table"))
}

#' @noRd
#' @keywords internal
.fb_stream_report_progress <- function(i, n_chunks, n_total, k_now) {
  cli::cli_inform(
    paste0(
      "chunk {i}/{if (is.na(n_chunks)) '?' else as.character(n_chunks)}: ",
      "{format(n_total, big.mark = ' ')} rows read, K = {k_now} cells"
    )
  )
}


# ---------------------------------------------------------------- #
# Finalisation into an <fb_aggregated>                              #
# ---------------------------------------------------------------- #
# Builds the same object shape the in-memory `.fb_aggregate_gaussian()`
# returns, so the existing emit consumes a streamed aggregation
# unchanged. The fixed-effect contrast columns are reconstructed once,
# on the K unique cells, via model.matrix() -- the only place a design
# matrix is built, and it is K x p, never N x p.

#' @noRd
#' @keywords internal
.fb_stream_finalise <- function(fb, acc, key_cols, n_total, family) {
  fixed_form <- .fb_aggregate_fixed_formula(fb)

  random_cols <- character(0L)
  for (t in fb$random_terms) {
    if (identical(t$type %||% "simple", "simple")) {
      random_cols <- c(random_cols, as.character(t$var))
    }
  }

  cell_df <- as.data.frame(acc)[, key_cols, drop = FALSE]
  if (is.null(fixed_form)) {
    X <- matrix(
      1,
      nrow = nrow(acc),
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )
  } else {
    X <- stats::model.matrix(fixed_form, data = cell_df)
  }
  fixed_cols <- colnames(X)

  ss <- data.table::as.data.table(X)
  data.table::setnames(ss, fixed_cols)
  for (g in random_cols) {
    data.table::set(ss, j = g, value = acc[[g]])
  }

  stat_cols <- setdiff(names(acc), key_cols)
  for (cn in stat_cols) {
    data.table::set(ss, j = cn, value = acc[[cn]])
  }

  K <- nrow(ss)
  data.table::set(ss, j = "cell_id", value = seq_len(K))
  data.table::setcolorder(
    ss,
    c("cell_id", fixed_cols, random_cols, stat_cols)
  )

  cell_design <- as.matrix(as.data.frame(ss)[, fixed_cols, drop = FALSE])
  colnames(cell_design) <- fixed_cols

  structure(
    list(
      cell_design = cell_design,
      sufficient_stats = ss,
      cell_key_cols = c(fixed_cols, random_cols),
      fixed_cols = fixed_cols,
      random_cols = random_cols,
      K = as.integer(K),
      N = n_total,
      compression = as.numeric(K) / as.numeric(n_total),
      response = fb$response,
      family = family,
      streamed = TRUE
    ),
    class = c("fb_aggregated", "list")
  )
}


# ---------------------------------------------------------------- #
# Scope gate (streaming variant)                                    #
# ---------------------------------------------------------------- #
# The streaming path admits the three exponential-family cases gaussian
# / binomial / poisson; the term-shape refusals (continuous fixed
# effect, smooth, random slope, structured covariance) match
# `.assert_aggregate_in_scope()` so a model refused in memory is refused
# the same way when streamed.

#' @noRd
#' @keywords internal
.assert_aggregate_in_scope_stream <- function(fb) {
  ok_family <- c("gaussian", "binomial", "poisson")
  if (!fb$family %in% ok_family) {
    .stop_aggregate_out_of_scope(
      reason_code = "non_aggregatable_family",
      detail = paste0(
        "family = '",
        fb$family,
        "'; streaming aggregation supports ",
        paste(ok_family, collapse = ", "),
        "."
      )
    )
  }
  if (
    identical(fb$family, "gaussian") &&
      !(is.null(fb$link) || identical(fb$link, "identity"))
  ) {
    .stop_aggregate_out_of_scope(
      reason_code = "non_identity_link",
      detail = paste0(
        "link = '",
        fb$link,
        "'; gaussian aggregation requires identity."
      )
    )
  }

  for (t in fb$fixed_terms) {
    ttype <- t$type %||% "expression"
    if (ttype %in% c("smooth", "s", "t2", "smooth_mgcv", "spline")) {
      .stop_aggregate_out_of_scope(
        reason_code = "smooth_term_not_aggregatable",
        detail = paste0(
          "smooth fixed-effect term '",
          t$var %||% "?",
          "' breaks the cell-constant linear predictor."
        )
      )
    }
    if (
      ttype %in%
        c(
          "numeric",
          "continuous",
          "I",
          "expression",
          "polynomial",
          "factor_numeric_interaction"
        )
    ) {
      .stop_aggregate_out_of_scope(
        reason_code = "continuous_cell_key_data_dependent",
        detail = paste0(
          "fixed term '",
          t$label %||% t$var %||% "?",
          "' is continuous; it forces one cell per ",
          "row, so aggregation does not compress. ",
          "Use flexybayes() per row."
        )
      )
    }
  }

  for (t in fb$random_terms) {
    rtype <- t$type %||% "simple"
    if (identical(rtype, "simple")) {
      next
    }
    reason <- if (rtype %in% c("smooth", "s", "t2", "smooth_mgcv", "spline")) {
      "smooth_term_not_aggregatable"
    } else if (rtype %in% c("simple_slope_uncor", "slope", "random_slope")) {
      "random_slope_not_aggregatable"
    } else {
      "structured_random_not_aggregatable"
    }
    .stop_aggregate_out_of_scope(
      reason_code = reason,
      detail = paste0(
        "random-effect term type = '",
        rtype,
        "' breaks the cell-constant linear predictor."
      )
    )
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------- #
# Emit dispatch                                                     #
# ---------------------------------------------------------------- #

#' @noRd
#' @keywords internal
.fb_stream_emit <- function(
  fb,
  agg,
  head,
  backend,
  fixed,
  random,
  family,
  the_call,
  verbose,
  ...
) {
  data_name <- "streamed"

  if (isTRUE(verbose)) {
    cli::cli_inform(
      paste0(
        "fitting {agg$K} cells ",
        "(compression {sprintf('%.0f:1', 1 / agg$compression)}) ",
        "on backend '{backend}'"
      )
    )
  }

  fit <- if (identical(family, "gaussian")) {
    emit_gaussian_aggregated(
      fb = fb,
      fb_aggregated = agg,
      data = head,
      backend = backend,
      verbose = verbose,
      the_call = the_call,
      fixed = fixed,
      random = random,
      family = family,
      data_name = data_name,
      ...
    )
  } else {
    emit_count_aggregated(
      fb = fb,
      fb_aggregated = agg,
      data = head,
      backend = backend,
      verbose = verbose,
      the_call = the_call,
      fixed = fixed,
      random = random,
      family = family,
      data_name = data_name,
      ...
    )
  }

  if (!is.null(fit$extras) && !is.null(fit$extras$aggregation_meta)) {
    fit$extras$aggregation_meta$streamed <- TRUE
  }
  fit
}
