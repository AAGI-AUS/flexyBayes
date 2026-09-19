# the registered pairs and their axes are stable

    Code
      for (k in sort(ls(reg))) {
        entry <- get(k, envir = reg, inherits = FALSE)
        cat(paste(entry$pair, collapse = "--"), ": ", paste(entry$axes, collapse = " + "),
        "\n", sep = "")
      }
    Output
      brms--inla: algorithmic + implementation

