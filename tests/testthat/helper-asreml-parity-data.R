# helper-asreml-parity-data.R -- the data for the stored ASReml parity fits.
#
# data-raw/asreml_parity.R sources this file to build the data it hands to
# ASReml, and test-asreml-parity.R builds the same frames here, so both sides
# fit byte-identical rows. Each frame's MD5 is stored with its fit and
# checked before any comparison.

# One frame per shape, from agridat 1.26.
.ap_data <- function(shape) {
  e <- new.env()
  if (shape %in% c("01_iid_intercept", "02_nested_stratum", "15_weights")) {
    utils::data("john.alpha", package = "agridat", envir = e)
    d <- e$john.alpha[, c("plot", "rep", "block", "gen", "yield")]
    d$rep <- factor(d$rep)
    d$block <- factor(d$block)
    d$gen <- factor(d$gen)
    if (identical(shape, "15_weights")) {
      d$wt <- rep(c(1, 2, 4, 2), length.out = nrow(d))
    }
  } else if (identical(shape, "04_dsum_residual")) {
    utils::data("yan.winterwheat", package = "agridat", envir = e)
    d <- e$yan.winterwheat[, c("gen", "env", "yield")]
    d$gen <- factor(d$gen)
    d$env <- factor(d$env)
  } else {
    utils::data("besag.met", package = "agridat", envir = e)
    d <- e$besag.met[!is.na(e$besag.met$yield), ]
    d <- data.frame(
      env = factor(d$county),
      rep = factor(d$rep),
      gen = factor(d$gen),
      yield = d$yield
    )
    n_env <- switch(shape, "06_us_dim2" = 2L, "07_us_dim3" = 3L, 6L)
    d <- d[d$env %in% levels(d$env)[seq_len(n_env)], ]
    d <- droplevels(d)
  }
  ord <- if ("env" %in% names(d)) order(d$env, d$gen) else seq_len(nrow(d))
  d <- d[ord, ]
  rownames(d) <- NULL
  d
}

# MD5 of the frame as text, which does not depend on R's serialisation
# version.
.ap_data_md5 <- function(d) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(format(d, digits = 15L), path, row.names = FALSE)
  unname(tools::md5sum(path))
}

