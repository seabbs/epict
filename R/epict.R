#' @title Model Cycle Thresholds
#'
#' @description Provides a user friendly interface around package functionality
#' to model cycle tresholds, and symptom onsets from observed preprocessed
#' data, a piecewise linear model, and a global adjustment model.
#'
#' @param as_data_list PARAM_DESCRIPTION
#'
#' @param inits PARAM DESCRIPTION
#'
#'
#' @param ... Additional arguments passed to [cmdstanr].
#'
#' @return A [cmdstanr] model fit containing posterior samples
#'
#' @inheritParams epict_to_stan
#' @inheritParams epict_inits
#' @family epict
#' @export
epict <- function(obs,
                  model =  epict_model(),
                  as_data_list = epict_to_stan,
                  inits = epict_inits,
                  ct_model = subject_design(~1, obs),
                  adjustment_model = test_design(~1, obs),
                  individual_variation = 0.2, individual_correlation = 1,
                  censoring_threshold = 40, switch = TRUE,
                  onsets = TRUE, incubation_period = get_inc_period(),
                  likelihood = TRUE, output_loglik = FALSE, ...) {
  stan_data <- as_data_list(
    obs,
    ct_model = ct_model,
    adjustment_model = adjustment_model,
    individual_variation = individual_variation,
    individual_correlation = individual_correlation,
    censoring_threshold = censoring_threshold,
    switch = switch,
    onsets = onsets, incubation_period = incubation_period,
    likelihood = likelihood, output_loglik = output_loglik
  )

  fit <- model$sample(
    data = stan_data,
    init = inits(stan_data),
    ...
  )
  class(fit) <- c("epict", class(fit))
  return(fit)
}
