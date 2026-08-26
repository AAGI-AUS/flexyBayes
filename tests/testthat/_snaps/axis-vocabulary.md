# axis-label rendering is stable across the registered pairs

    Code
      for (k in sort(ls(reg))) {
        entry <- get(k, envir = reg, inherits = FALSE)
        cat(paste(entry$pair, collapse = "--"), ": ", flexyBayes:::.format_independence_axes(
          entry$axes), "\n", sep = "")
      }
    Output
      brms--inla: algorithmic + implementation

