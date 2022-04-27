#' Format observed symptom onset data for use with stan
#'
#' The package stan code (accessed using [epict_model()]) has
#' more extensive documentation.
#' 
#' @param onsets Logical, defaults to `TRUE`. Should symptom onsets 
#' observations be included in the model if available.
#' 
#' @return A list as required by stan.
#' @inheritParams epict_check_obs
#' @family modeltools
#' @importFrom data.table copy
#' @author Sam Abbott
#' @export
epict_onset_obs_as_list <- function(obs, onsets = TRUE) {

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
#' Formats observations for use in the stan model. Note that the
#' censoring limit is set internally in this function based on the minimum
#' cycle threshold value present for a censored test. The package stan code
#' (accessed using [epict_model()]) has more extensive documentation.
#' 
#' @return A list as required by stan.
#' @inheritParams epict_check_obs
#' @family modeltools
#' @importFrom data.table copy
#' @author Sam Abbott
#' @export
epict_obs_as_list <- function(obs) {

  obs <- data.table::copy(obs)[order(id)]
  obs[, obs_id := 1:.N]
  tests_per_id <- obs[, .(n = .N), by = "id"]$n

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

#' Select the parameters to adjust in the piecewise model
#'
#' @param params A character string indicating the paramters
#' to adjust in the piecewise model. Defaults to "all". Options 
#' are: the time at peak ("t_p"), time at switch ("t_p"),
#' time at clearance ("t_clear"), Cycle threshold (Ct) at peak ("c_p"),
#' Ct at switch ("c_s"), incubation period mean ("inc_mean"),
#' and incubation period standard deviation ("inc_sd"). 
#'
#' @return A named list of parameters to adjust in the piecewise 
#' model.
#' @importFrom purrr map
#' @family modeltools
#' @author Sam Abbott
#' @export
#' @examples
#' select_piecewise_parameters()
select_piecewise_parameters <- function(params = "all") {
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

#' Specify the cycle threshold adjustment model formula 
#'
#' @param formula A model formula defaulting to `~1`.
#' 
#' @param beta_default A vector of length two containing the 
#' default mean and standard deviation to use as the prior 
#' for all covariate effect sizes.
#' 
#' @return A named list including the design matrix ("design") 
#' and a `data.table` of priors for covariate effects ("beta").
#' @inheritParams epict_check_obs
#' @inheritParams select_piecewise_parameters
#' @importFrom data.table CJ
#' @family modeltools
#' @export
#' @author Sam Abbott
piecewise_formula <- function(formula = ~1, obs, beta_default = c(0, 0.1),
                              params = "all") {
  params <- select_piecewise_parameters(params)

  subjects <- extract_subjects(obs)
  design <- model.matrix(formula, data = subjects)

  beta <- data.table::CJ(
    parameter = params,
    effect = colnames(design)[-1]
  )

  beta[, `:=`(
    mean = beta_default[1],
    sd = beta_default[2]
  )]

  out <- list(
    design = design, subjects = subjects, params = params,
    beta = beta[]
  )
  return(out)
}

#' Specify the cycle threshold adjustment model formula 
#' 
#' @return A named list including the design matrix ("design") 
#' and a `data.table` of priors for covariate effects ("beta").
#' @inheritParams piecewise_formula
#' @export
#' @importFrom data.table CJ
#' @family modeltools
#' @author Sam Abbott
#' @examples 
#' obs <- data.frame(
#'  age = c(1, 2, 3), cats = c(1, 2, 1), status = c("h", "l", "h")
#' )
#' adjustment_formula(~ cats + age + status, obs)
adjustment_formula <- function(formula = ~1, obs, beta_default = c(0, 0.1)) {
  design <- model.matrix(formula, data = obs)
  
  beta <- data.table::CJ(
    parameter = c("ct_shift", "ct_scale"),
    effect = colnames(design)[-1]
  )

  beta[, `:=`(
    mean = beta_default[1],
    sd = beta_default[2]
  )]

  out <- list(
    design = design, beta = beta[]
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
#' @inheritParams epict_onset_obs_as_list
#' @family modeltools
#' @author Sam Abbott
#' @export
#' @examples
#' # Default options
#' epict_model_opts()
#' 
#' # No variation and use the piecewise switch
#' epict_model_opts(switch = TRUE, variation = "none")
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
    ind_corr = as.numeric(variation == "correlated"),
    K = ifelse(switch, 5, 3)
  )
  return(data)
}

#' Format formula data for use with stan
#' 
#' The package stan code (accessed using [epict_model()]) has
#' extensive more extensive documentation.
#' 
#' @param piecewise_formula A list describing the piecewise linear cycle threshold
#' formula as described by [piecewise_formula()].
#'
#' @param adjustment_formula A list describing the cycle threshold linear
#' adjustment formula (shift and scale) as described by [adjustment_formula()].
#'
#' @param model_opts A list of Cycle threshold model options as returned
#' by [epict_ct_model_opts()].
#' 
#' @return A list as required by stan.
#' @importFrom data.table as.data.table
#' @family modeltools
#' @author Sam Abbott
#' @export
epict_formula_as_list <- function(piecewise_formula, adjustment_formula,
                                  model_opts = epict::model_opts()) {
  piecewise_formula$beta <- data.table::as.data.table(piecewise_formula$beta)
  adjustment_formula$beta <- data.table::as.data.table(adjustment_formula$beta)

  data <- list(
      preds = ncol(piecewise_formula$design) - 1,
      design = piecewise_formula$design,
      adj_t_p = piecewise_formula$params[["t_p"]],
      beta_t_p_m = piecewise_formula$beta[parameter == "t_p"]$mean,
      beta_t_p_sd = piecewise_formula$beta[parameter == "t_p"]$sd,
      adj_t_s = min(
        piecewise_formula$params[["t_s"]], model_opts$switch
      ),
      beta_t_s_m = piecewise_formula$beta[parameter == "t_s"]$mean,
      beta_t_s_sd = piecewise_formula$beta[parameter == "t_s"]$sd,
      adj_t_clear = piecewise_formula$params[["t_clear"]],
      beta_t_clear_m = piecewise_formula$beta[parameter == "t_clear"]$mean,
      beta_t_clear_sd = piecewise_formula$beta[parameter == "t_clear"]$sd,
      adj_c_p = piecewise_formula$params[["c_p"]],
      beta_c_p_m = piecewise_formula$beta[parameter == "c_p"]$mean,
      beta_c_p_sd = piecewise_formula$beta[parameter == "c_p"]$sd,
      adj_c_s = min(
        piecewise_formula$params[["c_s"]], model_opts$switch
      ),
      beta_c_s_m = piecewise_formula$beta[parameter == "c_s"]$mean,
      beta_c_s_sd = piecewise_formula$beta[parameter == "c_s"]$sd,
      adj_inc_mean = piecewise_formula$params[["inc_mean"]],
      beta_inc_mean_m = piecewise_formula$beta[parameter == "inc_mean"]$mean,
      beta_inc_mean_sd = piecewise_formula$beta[parameter == "inc_mean"]$sd,
      adj_inc_sd = piecewise_formula$params[["inc_sd"]],
      beta_inc_sd_m = piecewise_formula$beta[parameter == "inc_sd"]$mean,
      beta_inc_sd_sd = piecewise_formula$beta[parameter == "inc_sd"]$sd
  )

  data <- c(fdata, list(
    adj_ct = as.numeric(
      (ncol(adjustment_formula$design) - 1) > 0
    ),
    ct_preds = ncol(adjustment_formula$design) - 1,
    ct_design = adjustment_formula$design,
    beta_ct_shift_m = adjustment_formula$beta[parameter == "ct_shift"]$mean,
    beta_ct_shift_sd = adjustment_formula$beta[parameter == "ct_shift"]$sd,
    bt_ct_scale_m = adjustment_formula$beta[parameter == "ct_scale"]$mean,
    beta_scale_sd = adjustment_formula$beta[parameter == "ct_scale"]$sd
  ))
  return(data)
}

#' Format population-level priors for use in stan
#'
#' @param priors A data.table of population-level priors as produced
#' by [epict_priors()].
#'
#' @return A list of named priors with the mean and standard deivation
#' for each.
#' 
#' @family modeltools
#' @importFrom data.table copy
#' @importFrom purrr map
#' @export
#' @author Sam Abbott
#' @examples
#' priors <- epict_priors()
#' epict_population_priors_as_list(priors)
epict_population_priors_as_list <- function(priors) {
  priors <- data.table::copy(priors)
  priors[, variable := paste0(variable, "_p")]
  priors <- priors[, .(variable, intercept_mean, intercept_sd)]
  priors <- split(priors, by = "variable", keep.by = FALSE)
  priors <- purrr::map(priors, ~ as.vector(t(.)))
  priors <- purrr::map(priors, ~ .[!is.na(.)])
  return(priors)
}

#' Format individual-level priors for use in stan
#'
#' @param priors A data.table of individual-level priors as produced
#' by [epict_priors()].
#'
#' @return A list of named priors with the mean and standard deivation
#' for each.
#' 
#' @family modeltools
#' @importFrom data.table copy
#' @importFrom purrr map
#' @export
#' @author Sam Abbott
#' @examples
#' priors <- epict_priors()
#' epict_individual_priors_as_list(priors)
epict_individual_priors_as_list <- function(priors) {
  priors <- data.table::copy(priors)
  default_priors <- epict::epict_priors()
  default_priors[, .(variable)]
  priors <- default_priors[priors, on = "variable"]
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
#' @author Sam Abbott
#' @export
#' @examples 
#' epict_inference_opts()
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
#' @author Sam Abbott
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
