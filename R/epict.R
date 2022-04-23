#' Model Cycle Thresholds
#'
#' Provides a user friendly interface around package functionality
#' to model cycle thresholds, and symptom onsets from observed preprocessed
#' data, a piecewise linear model, and a linear adjustment model.
#'
#' @param convert_to_list Convert observations and model settings into
#' a list for use in stan. Defaults to using [epict_convert_to_list()].
#'
#' @param inits A function that combined with a list of data returns 
#' sample initial conditions for use during model fitting based
#' on model priors. Defaults to [epict_inits()]
#'
#' @param ... Additional arguments passed to [cmdstanr].
#'
#' @return A [cmdstanr] model fit containing posterior samples
#'
#' @inheritParams epict_convert_to_list
#' @inheritParams epict_inits
#' @family epict
#' @export
epict <- function(obs,
                  model = epict_model(),
                  piecewise_formula = epict::piecewise_formula(~1, obs),
                  adjustment_formula = epict::adjustment_formula(~1, obs),
                  priors = epict::epict_priors(),
                  model_opts = epict::epict_model_opts()
                  inference_opts = epict::epict_inference_opts()
                  convert_to_list = epict_convert_to_list,
                  inits = epict_inits,
                  ...) {
  model_list <- convert_to_list(
    obs,  piecewise_formula = piecewise_formula, 
    adjustment_formula = adjustment_formula, priors = priors,
    model_opts = model_opts, inference_opts = inference_opts
  )

  fit <- model$sample(
    data = model_list,
    init = inits(model_list),
    ...
  )
  class(fit) <- c("epict", class(fit))
  return(fit)
}
