# The refusal API for plug-in packages: fb_register_refusals() adds a
# package's codes to its own partition of the one refusal registry,
# fb_refuse() raises any registered code, and fb_refusals() lists every
# partition with a `package` column.
#
# Fake plug-ins are registered under made-up package names and removed when
# each test ends, so the live registry other tests read is unchanged.

.fake_refusals <- function(codes, description = "A fake plug-in refusal.") {
  data.frame(
    reason_code = codes,
    description = rep(description, length(codes)),
    since_version = rep("0.0.1", length(codes)),
    stringsAsFactors = FALSE
  )
}

.drop_partition_on_exit <- function(package, env = parent.frame()) {
  withr::defer(
    flexyBayes:::.unregister_refusal_partition(package),
    envir = env
  )
}

# An environment R treats as the namespace of `name`, so a call evaluated in
# it sees `name` as its calling package.
.fake_namespace <- function(name) {
  ns <- new.env(parent = globalenv())
  info <- new.env(parent = baseenv())
  assign("spec", c(name = name, version = "0.0.1"), envir = info)
  assign(".__NAMESPACE__.", info, envir = ns)
  ns
}

test_that("a plug-in partition is listed by fb_refusals() under its package", {
  .drop_partition_on_exit("fakePluginA")
  fb_register_refusals(
    .fake_refusals(c("fakea_widget_missing", "fakea_gadget_unfit")),
    package = "fakePluginA"
  )

  tab <- fb_refusals()
  expect_true("package" %in% names(tab))
  expect_identical(
    tab$package[tab$reason_code %in% c("fakea_widget_missing", "fakea_gadget_unfit")],
    c("fakePluginA", "fakePluginA")
  )
  expect_identical(
    tab$package[tab$reason_code == "unsupported_family"],
    "flexyBayes"
  )

  one <- fb_refusals(reason_code = "fakea_widget_missing")
  expect_identical(nrow(one), 1L)
  expect_identical(one$description, "A fake plug-in refusal.")
  expect_identical(one$since_version, "0.0.1")
  expect_identical(one$plan_field, NA_character_)
})

test_that("fb_register_refusals() takes the calling package by default", {
  .drop_partition_on_exit("fakePluginDefault")
  ns <- .fake_namespace("fakePluginDefault")
  assign("tab", .fake_refusals("fakedefault_thing_absent"), envir = ns)
  eval(quote(flexyBayes::fb_register_refusals(tab)), envir = ns)

  tab <- fb_refusals(reason_code = "fakedefault_thing_absent")
  expect_identical(tab$package, "fakePluginDefault")
})

test_that("fb_register_refusals() outside a package asks for `package`", {
  outside <- new.env(parent = globalenv())
  assign("tab", .fake_refusals("fake_orphan"), envir = outside)
  expect_error(
    eval(quote(flexyBayes::fb_register_refusals(tab)), envir = outside),
    "outside a package"
  )
  expect_false("fake_orphan" %in% fb_refusals()$reason_code)
})

test_that("fb_refuse() raises a flexybayes_refusal carrying the code and package", {
  .drop_partition_on_exit("fakePluginA")
  fb_register_refusals(
    .fake_refusals("fakea_widget_missing", "The widget is missing."),
    package = "fakePluginA"
  )

  cnd <- tryCatch(
    fb_refuse(
      "fakea_widget_missing",
      "No widget in `data`.",
      widget = "w1",
      class = "fakepluginA_error"
    ),
    error = function(e) e
  )
  expect_s3_class(cnd, "flexybayes_refusal")
  expect_s3_class(cnd, "flexybayes_refusal_fakea_widget_missing")
  expect_s3_class(cnd, "fakepluginA_error")
  expect_identical(
    class(cnd),
    c(
      "flexybayes_refusal_fakea_widget_missing",
      "fakepluginA_error",
      "flexybayes_refusal",
      "error",
      "condition"
    )
  )
  expect_identical(cnd$reason_code, "fakea_widget_missing")
  expect_identical(cnd$package, "fakePluginA")
  expect_identical(cnd$widget, "w1")
  expect_identical(conditionMessage(cnd), "No widget in `data`.")
})

test_that("fb_refuse() with defaults uses the registered description", {
  .drop_partition_on_exit("fakePluginA")
  fb_register_refusals(
    .fake_refusals("fakea_widget_missing", "The widget is missing."),
    package = "fakePluginA"
  )
  expect_error(
    fb_refuse("fakea_widget_missing"),
    "The widget is missing.",
    fixed = TRUE,
    class = "flexybayes_refusal_fakea_widget_missing"
  )
})

test_that("fb_refuse() raises core codes, which carry package flexyBayes", {
  cnd <- tryCatch(fb_refuse("unsupported_family"), error = function(e) e)
  expect_s3_class(cnd, "flexybayes_refusal_unsupported_family")
  expect_identical(cnd$package, "flexyBayes")
})

test_that("fb_refuse() refuses an unregistered code", {
  expect_error(
    fb_refuse("fake_never_registered", "x"),
    "not registered"
  )
  err <- tryCatch(fb_refuse("fake_never_registered"), error = function(e) e)
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "fb_register_refusals()", fixed = TRUE)
})

test_that("fb_refuse() refuses a field that would replace a standard one", {
  expect_error(
    fb_refuse("unsupported_family", "x", package = "plotly"),
    "package"
  )
})

test_that("a code already in the core partition is refused", {
  .drop_partition_on_exit("fakePluginA")
  expect_error(
    fb_register_refusals(
      .fake_refusals(c("fakea_fresh_code", "unsupported_family")),
      package = "fakePluginA"
    ),
    "unsupported_family"
  )
  # Nothing from a refused table is registered.
  expect_false("fakea_fresh_code" %in% fb_refusals()$reason_code)
  expect_identical(
    fb_refusals(reason_code = "unsupported_family")$package,
    "flexyBayes"
  )
})

test_that("a code already in another plug-in partition is refused", {
  .drop_partition_on_exit("fakePluginA")
  .drop_partition_on_exit("fakePluginB")
  fb_register_refusals(.fake_refusals("fake_shared_code"), package = "fakePluginA")
  expect_error(
    fb_register_refusals(.fake_refusals("fake_shared_code"), package = "fakePluginB"),
    "fakePluginA"
  )
  expect_identical(
    fb_refusals(reason_code = "fake_shared_code")$package,
    "fakePluginA"
  )
})

test_that("a package re-registering its own code is accepted only unchanged", {
  .drop_partition_on_exit("fakePluginA")
  fb_register_refusals(.fake_refusals("fakea_widget_missing"), package = "fakePluginA")
  expect_no_error(
    fb_register_refusals(.fake_refusals("fakea_widget_missing"), package = "fakePluginA")
  )
  expect_error(
    fb_register_refusals(
      .fake_refusals("fakea_widget_missing", "A different description."),
      package = "fakePluginA"
    ),
    "already registered"
  )
  expect_identical(
    fb_refusals(reason_code = "fakea_widget_missing")$description,
    "A fake plug-in refusal."
  )
})

test_that("a code given twice in one table is refused", {
  .drop_partition_on_exit("fakePluginA")
  expect_error(
    fb_register_refusals(
      .fake_refusals(c("fakea_twice", "fakea_twice")),
      package = "fakePluginA"
    ),
    "fakea_twice"
  )
})

test_that("the core partition cannot be modified", {
  n_core <- length(ls(flexyBayes:::.refusal_registry))
  expect_error(
    fb_register_refusals(.fake_refusals("fake_into_core"), package = "flexyBayes"),
    "flexyBayes"
  )
  expect_error(
    assign(
      "fake_into_core",
      list(),
      envir = flexyBayes:::.refusal_registry
    )
  )
  expect_true(environmentIsLocked(flexyBayes:::.refusal_registry))
  expect_identical(length(ls(flexyBayes:::.refusal_registry)), n_core)
  expect_false("fake_into_core" %in% fb_refusals()$reason_code)
})

test_that("a plug-in entry cannot be overwritten in place", {
  .drop_partition_on_exit("fakePluginA")
  fb_register_refusals(.fake_refusals("fakea_widget_missing"), package = "fakePluginA")
  partition <- get(
    "fakePluginA",
    envir = flexyBayes:::.refusal_plugin_registry
  )
  expect_error(assign("fakea_widget_missing", list(), envir = partition))
})

test_that("fb_register_refusals() validates the table", {
  expect_error(
    fb_register_refusals(list(reason_code = "x"), package = "fakePluginA"),
    "data frame"
  )
  expect_error(
    fb_register_refusals(
      data.frame(reason_code = "fakea_x", stringsAsFactors = FALSE),
      package = "fakePluginA"
    ),
    "description"
  )
  bad <- .fake_refusals("fakea_x")
  bad$description <- NA_character_
  expect_error(fb_register_refusals(bad, package = "fakePluginA"), "description")
  expect_error(
    fb_register_refusals(.fake_refusals("Not A Code"), package = "fakePluginA"),
    "reason_code"
  )
  expect_error(
    fb_register_refusals(.fake_refusals("fakea_x"), package = "not a package"),
    "package"
  )
  expect_false("fakea_x" %in% fb_refusals()$reason_code)
})

test_that("the selection layer registers in the flexySelect partition", {
  tab <- fb_refusals()
  sel <- tab$reason_code[tab$package == "flexySelect"]
  expect_true("selection_seed_required" %in% sel)
  expect_true("penalty_block_overlap" %in% sel)
  expect_true("loss_surface_requires_plotly" %in% sel)
  expect_true("selection_loss_surface_undefined" %in% sel)
  expect_false("loss_surface_requires_plotly" %in% ls(flexyBayes:::.refusal_registry))

  cnd <- tryCatch(
    stop(flexyBayes:::.fb_refusal_condition(
      reason_code = "selection_seed_required",
      message = "probe"
    )),
    error = function(e) e
  )
  expect_s3_class(cnd, "flexybayes_refusal_selection_seed_required")
  expect_identical(cnd$package, "flexySelect")
})
