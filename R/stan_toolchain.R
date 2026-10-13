# stan_toolchain.R, telling a broken compiler apart from a broken model.
#
# The brms path compiles a Stan program, so it needs a working C++
# toolchain. When that toolchain is missing or mismatched, a macOS
# Command Line Tools / SDK skew is the common one, the compiler's own
# message ("fatal error: 'cmath' file not found") goes to the console and
# what reaches R is a downstream symptom: `invalid connection`. Reported
# verbatim, that sends a user looking for a fault in their model, and there
# is none. The R Journal review of 0.10.0 arrived carrying exactly this.
#
# The fix is not to pattern-match the message, which would mean guessing at
# strings from compilers and platforms this machine cannot produce. It is
# to ASK the toolchain, once, on the failure path: compile a three-line
# program that includes the header the failure implicates. The answer is
# evidence about this machine instead of a reading of an error string.

# Compile a minimal C++ translation unit through R's own build machinery.
# TRUE when the toolchain works, FALSE when it does not. Runs only after a
# fit has already failed, so its cost never falls on a working install.
.fb_cxx_toolchain_works <- function() {
  dir <- file.path(tempdir(), paste0("fb_cxx_", as.integer(runif(1L, 1, 1e8))))
  if (!dir.create(dir, showWarnings = FALSE, recursive = TRUE)) {
    return(NA)
  }
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  src <- file.path(dir, "fb_probe.cpp")
  writeLines(
    c(
      "#include <cmath>",
      "extern \"C\" double fb_probe(double x) { return std::sqrt(x); }"
    ),
    src
  )
  status <- tryCatch(
    suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", shQuote(src)),
      stdout = FALSE,
      stderr = FALSE
    )),
    error = function(e) NA_integer_
  )
  if (is.na(status)) NA else identical(as.integer(status), 0L)
}

# Raise the failure a brms fit ended on. When the toolchain probe says
# this machine cannot compile at all, the failure gets a name and a
# repair, when it compiles fine, the fault is somewhere else and the
# engine's own words are what the user needs. Naming a toolchain fault on
# a fit that failed for another reason would send the reader to fix
# something that is not broken, which is worse than saying less.
.fb_stop_brms_failure <- function(cond_message) {
  if (isFALSE(.fb_cxx_toolchain_works())) {
    stop(.fb_refusal_condition(
      reason_code = "brms_stan_toolchain_compile_failure",
      message = paste0(
        "brms::brm() failed, and this machine cannot compile a three-line ",
        "C++ file either, so the fault is in the C++ toolchain, not in ",
        "the model.\n",
        "  Check it with pkgbuild::check_build_tools(debug = TRUE). On ",
        "macOS this is usually a Command Line Tools or SDK mismatch, ",
        "repaired by xcode-select --install.\n",
        "  backend = \"inla\" compiles nothing and needs no toolchain.\n",
        "  The engine reported: ",
        cond_message
      )
    ))
  }
  stop("brms::brm() failed: ", cond_message, call. = FALSE)
}
