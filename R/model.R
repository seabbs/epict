#' Default model priors
#'
#' @param individual_variation A vector of length two specifying the 
#' default mean and standard deviation of individual-level variation.
#' 
#' @param include_descriptive Logical, defaults `FALSE`. Should an extended
#' list of priors be returned including transformed parameters. This is useful
#' when transforming model output and for user exploration.
#' 
#' @return A data frame summarising the model priors.
#'
#' @importFrom data.table data.table
#' @family model
#' @export
#' @author Sam Abbott
#' @examples
#' epict_priors()
#' 
#' # Also include descriptive parameters
#' epict_priors(include_descriptive = TRUE)
epict_priors <- function(individual_variation = c(0, 0.05),
                         include_descriptive = FALSE) {
  priors <- data.table::data.table(
    variable = c(
      "t_inf",
      "c_thres",
      "c_p",
      "t_p",
      "c_s",
      "t_s",
      "t_clear",
      "inc_mean",
      "inc_sd",
      "sigma",
      "lkj",
      "ct_shift",
      "ct_scale",
      "nat_inc_mean",
      "nat_inc_sd"
    ),
    name = c(
      "Time between infection and first positive test",
      "Ct value at clearance of infection",
      "Ct value at peak",
      "Time at peak Ct",
      "Ct value at switch",
      "Time at switch Ct",
      "Time at clearance of infection",
      "Incubation period (log) mean",
      "Incubation period (log) standard deviation",
      "Observation standard deviation",
      "Lewandowski-Kurowicka-Joe distribution",
      "Ct intercept adjustment",
      "Ct multiplicative adjustment",
      "Incubation period (mean)",
      "Incubation period (standard deviation)"
    ),
    detail = c(
      "Offset by onset if available",
      "Offset by user specified limit of Ct detection",
      "Logit scale offset by latent limit of detection/Ct at switch",
      "Log scale",
      "Logit scale offset by latent limit of detection",
      "Log scale relative to the time at peak Ct",
      "Relative to time at peak and switch C",
      "Parameterised as a log-normal distribution",
      "Parameterised as a log-normal distribution",
      "Assuming a normal error model",
      "Prior used for the individual-level correlation matrix. There is 
      one parameter which at 1 represents a uniform prior. Smaller values indicate strong correlations, and larger values weaker correlations",
      "Cycle threshold intercept adjustment with an intercept set to 0 (i. no 
      adjustment. This is adjusted using `adjustment_formula()`",
      "Cycle threshold gradient adjustment with an intercept set to 1 (i. no 
      adjustment. This is adjusted using `adjustment_formula()`",
      "Incubation period mean on the natural scale",
      "Incubation period standard deviation on the natural scale"
    ),
    distribution = c(
      "Normal (truncated by onset or 0)",
      "Normal (truncated by user specified limit of detection)",
      "Normal",
      "Normal",
      "Normal",
      "Normal",
      "Normal",
      "Normal",
      "Zero truncated normal",
      "Zero truncated normal",
      "Lewandowski-Kurowicka-Joe (LKJ) distribution",
      "",
      "",
      "",
      ""
    ),
    intercept_mean = c(
      5, 10, 0, 1.61, 0, 1.61, 2.3, 1.62, 0.42, 2, 1, NA, NA, NA, NA
    ),
    intercept_sd = c(
      5, 10, 1, 0.5, 1, 0.5, 0.5, 0.06, 0.07, 2, NA, NA, NA, NA, NA
    )
  )

  priors[, `:=`(
    individual_variation_mean = NA_real_,
    individual_variation_sd = NA_real_
  )]

  params_with_var <- c(
    "c_p", "t_p", "t_s", "c_s", "t_clear"
  )
  priors[
    variable %in% params_with_var,
    individual_variation_mean := individual_variation[1]
  ]
  priors[
    variable %in% params_with_var,
    individual_variation_sd := individual_variation[2]
  ]

  if (!include_descriptive) {
    priors <- priors[!is.na(intercept_mean)]
  }
  return(priors[])
}

#' Format observations and model settings for use with stan
#'
#' @param model_opts A list of model options. See [epict_model_opts()]
#' for details
#' 
#' @param inference_opts A list of options to use for inference. See
#' [epict_inference_opts()] for details.
#' 
#' @return A list as required by stan.
#' @inheritParams epict_obs_as_list
#' @inheritParams epict_onset_obs_as_list
#' @inheritParams epict_model_opts
#' @inheritParams epict_inference_opts
#' @inheritParams epict_formula_as_list
#' @inheritParams epict_population_priors_as_list
#' @inheritParams epict_individual_priors_as_list
#' @inheritParams epict_posterior_as_prior
#' @family model
#' @author Sam Abbott
#' @export
epict_convert_to_list <- function(
  obs, 
  piecewise_formula = epict::piecewise_formula(~1, obs),
  adjustment_formula = epict::adjustment_formula(~1, obs),
  priors = epict::epict_priors(),
  model_opts = epict::epict_model_opts(),
  inference_opts = epict::epict_inference_opts()
) {
  data <- enw_obs_as_list(obs)

  data <- c(data, 
    epict_onset_obs_as_list(obs, onsets = (model_opts$onsets == 1))
  )

  data <- c(
    data,
    epict_formula_as_list(
      data,
      piecewise_formula = piecewise_formula,
      adjustment_formula = adjustment_formula
    )
  )

  data <- c(
    data,
    model_opts,
    inference_opts
  )

  data <- c(
    data,
    epict_population_priors_as_list(priors),
    epict_individual_priors_as_list(priors)
  )
  return(data)
}

#' Set up initial conditions for model
#'
#' @param data A list of data as produced by [epict_as_list()].
#'
#' @param scale Numeric defaults to 0.1. The numeric scaling to apply to the
#' standard deviation when sampling. Helpful for reducing numeric instability
#' during warmup.
#' 
#' @return A function that when called returns a list of initial conditions
#' for the package stan models.
#'
#' @family model
#' @importFrom purrr map_dbl
#' @importFrom truncnorm rtruncnorm
#' @author Sam Abbott
#' @export
epict_inits <- function(data, scale = 0.1) {
  function() {
    inits <- list(
      t_inf = purrr::map_dbl(
        1:data$P,
        ~ truncnorm::rtruncnorm(
          1,
          a = -data$onset_time[.],
          mean = -data$onset_time[.] + data$t_inf_p[1],
          sd = data$t_inf_p[2] * scale
        )
      ),
      c_thres = truncnorm::rtruncnorm(
        1,
        a = data$c_lod,
        mean = data$c_lod + data$c_thres_p[1],
        sd = data$c_thres_p[2] * scale
      ),
      c_p_int = rnorm(1, data$c_p_p[1], data$c_p_p[2]),
      t_p_int = rnorm(1, data$t_p_p[1], data$t_p_p[2]),
      t_clear_int = rnorm(1, data$t_clear_p[1], data$t_clear_p[2]),
      ind_var = abs(rnorm(data$K, 0, data$ind_var_sd * scale)),
      ind_eta = matrix(
        rnorm(data$P * data$K, 0, 1), nrow = data$K, ncol = data$P
      ),
      sigma = truncnorm::rtruncnorm(
        1, a = 0, mean = data$sigma_p[1], sd = data$sigma_p[2]
      )
    )

    if (data$switch > 0) {
      inits$c_s_int <- array(rnorm(1, data$c_s_p[1], data$c_s_p[2]))
      inits$t_s_int <- array(rnorm(1, data$t_s_p[1], data$t_s_p[2]))
    }

    if (data$preds > 0) {
      if (data$adj_t_p > 0) {
        inits$beta_t_p <- rnorm(
          data$preds, data$beta_t_p_m, data$beta_t_p_sd * scale
        )
      }
      if (data$adj_t_s > 0) {
        inits$beta_t_s <- rnorm(
          data$preds, data$beta_t_s_m, data$beta_t_s_sd * scale
        )
      }
      if (data$adj_t_clear > 0) {
        inits$beta_t_clear <- rnorm(
          data$preds, data$beta_t_clear_m, data$beta_t_clear_sd * scale
        )
      }
      if (data$adj_c_p > 0) {
        inits$beta_c_p <- rnorm(
          data$preds, data$beta_c_p_m, data$beta_c_p_sd * scale
        )
      }
      if (data$adj_c_s > 0) {
        inits$beta_c_s <- rnorm(
          data$preds, data$beta_c_s_m, data$beta_c_s_sd * scale
        )
      }
      if (data$adj_inc_mean > 0) {
        inits$beta_inc_mean <- rnorm(
          data$preds, data$beta_inc_mean_m, data$beta_inc_mean_sd * scale
        )
      }
      if (data$adj_inc_sd > 0) {
        inits$beta_inc_sd <- rnorm(
          data$preds, data$beta_inc_sd_m, data$beta_inc_sd_sd * scale
        )
      }
    }

    if (data$ct_preds > 0) {
      if (data$adj_ct > 0) {
        inits$beta_ct_shift <- rnorm(
          data$ct_preds, data$beta_ct_shift_m, data$beta_ct_shift_sd * scale
        )
        inits$beta_ct_scale <- rnorm(
          data$ct_preds, data$beta_ct_scale_m, data$beta_ct_scale_sd * scale
        )
      }
    }

    if (data$any_onsets == 1) {
      inits$inc_mean_int <- rnorm(
        1, data$inc_mean_p[1], data$inc_mean_p[2] * scale
     )
      inits$inc_sd_int <- truncnorm::rtruncnorm(
        1, a = 0, mean =  data$inc_sd_p[1], sd = data$inc_sd_p[2] * scale
      )
    }
    return(inits)
  }
}

#' Load and compile the Cycle thresholds model
#'
#' @param model A character string indicating the path to the model.
#' If not supplied the package default model is used.
#'
#' @param include A character string specifying the path to any stan
#' files to include in the model. If missing the package default is used.
#' 
#' @param compile Logical, defaults to `TRUE`. Should the model
#' be loaded and compiled using [cmdstanr::cmdstan_model()].
#'
#' @param threads Logical, defaults to `FALSE`. Should the model compile with
#' support for multi-thread support in chain. Note that this requires the use of
#' the `threads_per_chain` argument when model fitting using [epict()] and is
#' not currently supported in model.
#'
#' @param stanc_options A list of options to pass to the `stanc_options` of
#' [cmdstanr::cmdstan_model()]. By default "01" is passed which specifies simple
#' optimisations should be done by the prior to compilation.
#' 
#' @param verbose Logical, defaults to `TRUE`. Should verbose
#' messages be shown.
#' 
#' @param ... Additional arguments passed to [cmdstanr::cmdstan_model()].
#'
#' @return A `cmdstanr` model.
#'
#' @family model
#' @export
#' @importFrom cmdstanr cmdstan_model
#' @examplesIf interactive()
#' mod <- epict_model()
#' mod
epict_model <- function(model, include,
                        compile = TRUE, threads = FALSE, 
                        stanc_options = list("O1"), verbose = TRUE, ...) {
  if (missing(model)) {
    model <- "stan/model.stan"
    model <- system.file(model, package = "epict")
  }
  if (missing(include)) {
    include <- system.file("stan", package = "epict")
  }

  if (compile) {
    if (verbose) {
      model <- cmdstanr::cmdstan_model(model,
        include_paths = include,
        stanc_options = stanc_options,
        cpp_options = list(
          stan_threads = threads
        ),
        ...
      )
    } else {
      suppressMessages(
        model <- cmdstanr::cmdstan_model(model,
          include_paths = include,
          cpp_options = list(
            stan_threads = threads
          ), ...
        )
      )
    }
  }
  return(model)
}
