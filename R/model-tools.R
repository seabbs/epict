#' Format observed symptom onset data for use with stan
#'
#' @param onsets Logical, defaults to `TRUE`. Should symptom onsets 
#' observations be included in the model if available.
#' [epict_obs_as_data_list()].
#' @inheritParams epict_check_obs
#' @return A list as required by stan.
#' @family modeltools
#' @importFrom data.table copy
#' @export
epict_onset_obs_as_data_list <- function(obs, onsets = TRUE) {

  obs <- data.table::copy(obs)
  P <- length(unique(obs$id))

  if (is.null(obs$onset_time) | !onsets) {
    data <- list(
      any_onsets = 0,
      onset_avail = rep(0, P),
      onset_time = rep(0, P),
      onset_window = rep(0, P)
    )
  } else {
    onset_obs <- suppressWarnings(
      obs[,
        .(onset_time = min(onset_t_rel_uncensored, na.rm = TRUE), id),
        by = "id"
      ][
        is.infinite(onset_time), onset_time := NA
      ]
    )
    data <- list(
      any_onsets = 1,
      onset_avail = as.numeric(!is.na(onset_obs$onset_time)),
      nonsets = sum(as.numeric(!is.na(onset_obs$onset_time))),
      ids_with_onsets = onset_obs[!is.na(onset_time), id],
      onset_time = ifelse(is.na(onset_obs$onset_time),
                            0, onset_obs$onset_time),
      onset_window = rep(1, P)
    )
  }
  return(data)
}

#' Format observed data for use with stan
#'
#' 
#' @inheritParams epict_check_obs
#' @return A list as required by stan.
#' @family modeltools
#' @importFrom data.table copy
#' @export
epict_obs_as_data_list <- function(obs) {

  obs <- data.tabler::copy(obs)[order(id)]
  obs[, obs_id := 1:.N]
  tests_per_id <- obs[, .(n = .N), by = "id"]$n

  # Format indexing and observed data
  # See stan code for docs on what all of these are
  data <- list(
    N = obs[, .N],
    P = length(unique(obs$id)),
    id = obs[, id],
    tests_per_id = tests_per_id,
    cum_tests_per_id = cumsum(tests_per_id),
    day_rel = obs$t_rel_uncensored,
    ct_value = obs$ct_value,
    ncensored = length(obs[censored == 1]$obs),
    censored = obs[censored == 1]$obs,
    nuncensored = length(obs[censored == 0]$obs),
    uncensored = obs[censored == 0]$obs,
    uncensored_by_test = abs(obs$censored - 1),
    c_lod = min(obs[censored == 1]$ct_value)
  )
  return(data)
}

params_avail_to_adjust <- function(params = "all") {
  choices <- c("t_p", "t_s", "t_clear", "c_p", "c_s", "inc_mean", "inc_sd")
  params <- match.arg(params, c(choices, "all"), several.ok = TRUE)
  if (any(params %in% "all")) {
    params <- choices
  }
  params_list <- as.list(choices)
  names(params_list) <- choices
  params_list <- purrr::map(params_list, ~ as.numeric(any(params %in% .)))
  return(params_list)
}

test_design <- function(formula = ~1, data, preds_sd = 1) {
  design <- model.matrix(formula, data = data)

  out <- list(
    design = design, preds_sd = preds_sd
  )
  return(out)
}

subject_design <- function(formula = ~1, data, preds_sd = 0.1,
                                 params = "all") {
  params <- params_avail_to_adjust(params)

  subjects <- extract_subjects(data)
  design <- model.matrix(formula, data = subjects)

  out <- list(
    design = design, subjects = subjects, params = params,
    preds_sd = preds_sd
  )
  return(out)
}

#' Format formula data for use with stan
#' 
#' @param switch Logical, default to `FALSE`. Should a secondary breakpoint be
#' included in the piecewise linear Cycle threshold model.
#' 
#' @param latent_infections Logical, defaults to `TRUE`. Should latent time of
#' infection be modelled and used to adjust time from first positive test.
#' 
#' @param variation A character string indicating the type of individual level
#' variation to include. Defaults to "correlated" (a random effect with
#' modelled correlation structure). Other options include "uncorrelated" (
#' a random effect with no modelled correlation structure), and "none" (for no
#' individual level variation).
#' 
#' @return A list as required by stan.
#' @inheritParams epict_opts_as_data_list
#' @family modeltools
#' @export
epict_model_opts <- function(onsets = TRUE, switch = FALSE,
                             latent_infections = TRUE,
                             variation = "correlated") {
  variation <- match.arg(
    variation, 
    choices = c("none", "correlated", "uncorrelated")
  )

  data <- list(
    switch = as.numeric(switch),
    onsets = as.numeric(onsets),
    latent_inf = as.numeric(latent_infections),
    ind_var_m = as.numeric(variation != "none"),
    ind_corr = as.numeric(variation == "correlated")
  )
  return(opts)
}

#' Format formula data for use with stan
#' @param model_opts A list of Cycle threshold model options as returned
#' by [epict_ct_model_opts()].
#' 
#' @param piecewise_formula A list describing the piecewise linear cycle threshold
#' formula as described by [subject_design()].
#'
#' @param adjustment_model A list describing the cycle threshold linear
#' adjustment formula (shift and scale) as described by [test_design()].
#'
#' @return A list as required by stan.
#' @inheritParams epict_opts_as_data_list
#' @family modeltools
#' @export
epict_formula_as_data_list <- function(piecewise_formula,
                                       adjustment_formula,
                                       model_opts = epict::model_opts()) {
  data <- list(
      preds = ncol(piecewise_formula$design) - 1,
      preds_sd = piecewise_formula$preds_sd,
      design = piecewise_formula$design,
      adj_t_p = piecewise_formula$params[["t_p"]],
      adj_t_s = min(
        piecewise_formula$params[["t_s"]], model_opts$switch
      ),
      adj_t_clear = piecewise_formula$params[["t_clear"]],
      adj_c_p = piecewise_formula$params[["c_p"]],
      adj_c_s = min(
        piecewise_formula$params[["c_s"]], model_opts$switch
      ),
      adj_inc_mean = piecewise_formula$params[["inc_mean"]],
      adj_inc_sd = piecewise_formula$params[["inc_sd"]]
  )

  # Add adjustment formula data
  data <- c(fdata, list(
    adj_ct = as.numeric(
      (ncol(adjustment_formula$design) - 1) > 0
    ),
    ct_preds = ncol(adjustment_formula$design) - 1,
    ct_preds_sd = adjustment_formula$preds_sd,
    ct_design = adjustment_formula$design
  ))
  return(data)
}

#' Format population priors for use in stan
#'
#'
#' @param priors A data.table of population priors as produced
#' by [epict_priors()].
#'
#' @return A list of named priors with the mean and standard deivation
#' for each.
#' 
#' @family modeltools
#' @importFrom data.table copy
#' @importFrom purrr map
#' @export
#' @examples
#' priors <- epict_priors()
#' epict_population_priors_as_data_list(priors)
epict_population_priors_as_data_list <- function(priors) {
  priors <- data.table::copy(priors)
  priors[, variable := paste0(variable, "_p")]
  priors <- priors[, .(variable, intercept_mean, intercept_sd)]
  priors <- split(priors, by = "variable", keep.by = FALSE)
  priors <- purrr::map(priors, ~ as.vector(t(.)))
  return(priors)
}

#' Format 
#'
#'
#' @return RETURN_DESCRIPTION
#' @family modeltools
#' @importFrom data.table copy
#' @importFrom purrr map
#' @export
#' @examples
#' priors <- epict_priors()
#' epict_individual_priors_as_data_list(priors)
epict_individual_priors_as_data_list <- function(priors) {
  priors <- data.table::copy(priors)
  priors[, variable := paste0(variable, "_p")]
  priors <- priors[!is.na(individual_variation_mean)]
  priors <- priors[, 
    .(variable, individual_variation_mean, individual_variation_sd)
  ]
  priors <- list(ind_var_mean = priors$individual_variation_mean,
                 ind_var_sd = priors$individual_variation_sd)
  return(priors)
}

#' Format inference options for use with stan
#' 
#' @param pp Logical, defaults to `FALSE`. Should posterior predictions be made
#' for observed data. Useful for evaluating the performance of the model.
#'
#' @param likelihood Logical, defaults to `TRUE`. Should the likelihood be
#' included in the model
#'
#' @param output_loglik Logical, defaults to `FALSE`. Should the
#' log-likelihood be output. Disabling this will speed up fitting
#' if evaluating the model fit is not required.
#'
#' @param debug Logical, defaults to `FALSE`. Should within model debug
#' information be returned.
#' 
#' @return A list as required by stan.
#' @importFrom data.table fcase
#' @family modeltools
#' @export
epict_inference_opts <- function(pp = FALSE, likelihood = TRUE,
                                    debug = FALSE, output_loglik = FALSE) {

  data <- list(
    debug = as.numeric(debug),
    likelihood = as.numeric(likelihood),
    pp = as.numeric(pp),
    output_loglik = as.numeric(output_loglik)
  )
  return(data)
}

#' Extract summarised posteriors to use as priors
#'
#' Extract and summarise posteriors from a [epict()] to use as
#' priors (based on the distributions given in [epict_priors()])
#' assuming they are independent normal
#'
#' @param fit A [cmdstanr] fit as returned by [epict()]. 
#'
#' @param priors A data frame of priors to update as defined in
#' [epict_priors()].
#'
#' @param variables A character vector of variables both in the
#' posterior and in the default priors.
#'
#' @param sub A character string indicating the part of variable names in the
#' posterior to remove. This can be used to link variable names in the model
#' object to those expected as input.
#' 
#' @param scale Numeric, defaults to 5. Amount to scale posterior standard
#' deviations by.
#'
#' @return A data frame of priors
#' @family modeltools
#' @importFrom data.table setDT
#' @export
epict_posterior_as_prior <- function(fit, priors = epict::epict_priors(),
                                     variables = c(), sub = "_int",  scale = 5) {
  posteriors <- fit$summary(variables)
  posteriors <- setDT(posteriors)[, sd := sd * scale]
  posteriors <- posteriors[, .(variable, mean, sd)]
  posteriors <- posteriors[, variable := gsub(sub, "", variable)]
  priors <- priors[!(variable %in% variables)]
  priors <- rbind(priors, posteriors, fill = TRUE)
  return(priors[])
}