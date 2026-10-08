# Entry-argument refusals --- the checks every fitting and planning route
# runs on its arguments and on the parsed model, and the registry entries
# for them and for the aggregated routes' refusals (DEF-212 to DEF-218,
# DEF-221, DEF-223).

# .fb_refuse_aggregate_argument() --- `aggregate` that is not TRUE, FALSE or
# "auto".
.fb_refuse_aggregate_argument <- function(x) {
  stop(.fb_refusal_condition(
    reason_code = "aggregate_argument_invalid",
    message = paste0(
      "`aggregate` must be TRUE, FALSE, or \"auto\"; got: ",
      paste(deparse(x), collapse = " ")
    ),
    supplied = x
  ))
}

# .fb_check_weights() --- observation weights are one finite number per row.
.fb_check_weights <- function(weights, data) {
  n <- nrow(data)
  problem <- if (!is.numeric(weights)) {
    paste0("a ", class(weights)[[1L]], " vector, not a numeric one")
  } else if (length(weights) != n) {
    paste0(length(weights), " values for ", n, " rows")
  }
  if (is.null(problem)) {
    return(invisible(TRUE))
  }
  stop(.fb_refusal_condition(
    reason_code = "weights_invalid",
    message = paste0(
      "`weights` must be a numeric vector of length nrow(data); it is ",
      problem, "."
    )
  ))
}

# .fb_check_prior_argument() --- `prior` is an fb_prior() object or NULL.
.fb_check_prior_argument <- function(prior, arg = "prior") {
  if (is.null(prior) || inherits(prior, "fb_prior")) {
    return(invisible(TRUE))
  }
  stop(.fb_refusal_condition(
    reason_code = "prior_not_fb_prior",
    message = paste0(
      "`", arg, "` must be an fb_prior() object or NULL; it is a ",
      class(prior)[[1L]], ". Write the prior as ",
      "fb_prior(sd(group = \"g\") ~ half_normal(scale = 1)), for example."
    ),
    supplied_class = class(prior)
  ))
}

# .fb_refuse_absent_model_variables() --- a variable the parsed model reads
# that is not a column of the data. Without this the engines met it as
# "object 'x' not found", INLA after running.
.fb_refuse_absent_model_variables <- function(fb, data) {
  absent <- setdiff(.fb_model_vars(fb), names(data))
  if (!length(absent)) {
    return(invisible(TRUE))
  }
  stop(.fb_refusal_condition(
    reason_code = "model_variable_not_in_data",
    message = paste0(
      "The model names ",
      paste0("`", absent, "`", collapse = ", "),
      if (length(absent) == 1L) {
        ", which is not a column"
      } else {
        ", which are not columns"
      },
      " of the data. Correct the name in the formula, or add the column."
    ),
    variables = absent
  ))
}

# .fb_refuse_non_numeric_response() --- a factor or character response.
# No family flexyBayes fits takes one; on INLA a factor was fitted on its
# level codes. A binomial response keeps its own refusal.
.fb_refuse_non_numeric_response <- function(fb, data) {
  yv <- .fb_response_values(fb, data)
  if (!is.factor(yv) && !is.character(yv)) {
    return(invisible(TRUE))
  }
  if (identical(fb$family, "binomial")) {
    .refuse_non_binary_binomial_response(fb, data)
  }
  stop(.fb_refusal_condition(
    reason_code = "response_not_numeric",
    message = paste0(
      "The response `", fb$response, "` is a ", class(yv)[[1L]], ", and ",
      "family = \"", fb$family, "\" takes a numeric response. Recode it ",
      "to the numbers the family models."
    ),
    column = fb$response
  ))
}

# .fb_refuse_method_argument() --- an argument a sibling method takes that
# this method does not implement.
.fb_refuse_method_argument <- function(method, argument, why) {
  stop(.fb_refusal_condition(
    reason_code = "method_argument_not_supported",
    message = paste0(method, " does not take `", argument, "`: ", why),
    method = method,
    argument = argument
  ))
}


.populate_refusal_registry_entry_arguments <- function() {
  reg <- function(code, description) {
    .register_refusal(
      reason_code = code,
      description = description,
      message_template = "%s",
      registered_in_adr = "PLAN 0.16.0 DEF-212 to DEF-224",
      plan_field = NA_character_,
      since_version = "0.16.0"
    )
  }
  reg(
    "model_variable_not_in_data",
    "A fixed or random term names a variable that is not a column of the data."
  )
  reg(
    "aggregate_argument_invalid",
    "`aggregate` is not TRUE, FALSE or \"auto\"."
  )
  reg(
    "weights_invalid",
    "`weights` is not a numeric vector with one value per row of the data."
  )
  reg(
    "prior_not_fb_prior",
    "`prior` is neither an fb_prior() object nor NULL."
  )
  reg(
    "response_not_numeric",
    paste0(
      "The response is a factor or character vector, which no family ",
      "flexyBayes fits takes."
    )
  )
  reg(
    "stream_source_invalid",
    paste0(
      "flexybayes_stream() was given a `source` that is not an existing ",
      "file, a vector of .fst shards, a data frame or a chunk generator."
    )
  )
  reg(
    "method_argument_not_supported",
    paste0(
      "A method was given an argument that the same method on another fit ",
      "class takes and this one does not implement."
    )
  )
  reg(
    "aggregation_plan_ineligible",
    paste0(
      "`aggregate = TRUE` on a model the aggregation plan cannot compress; ",
      "the plan's reasons are carried in `reason_codes`."
    )
  )
  reg(
    "aggregation_out_of_scope",
    paste0(
      "The aggregated Gaussian route or flexybayes_stream() was given a model ",
      "outside the aggregation scope: a family or link it does not carry, a ",
      "smooth, a random slope, a structured random term or a continuous ",
      "fixed term. The gate's reason is carried in `scope_reason`."
    )
  )
}
