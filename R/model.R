#' Default model priors
#'
#' @return A data frame summarising the model priors.
#'
#' @importFrom data.table data.table
#' @family model
#' @export
#' @examples
#' @author Sam Abbott
#' epict_priors()
epict_priors <- function(individual_variation = c(0, 0.05)) {
  priors <- data.table::data.table(
    variable = c(
      "t_inf",
      "c_int",
      "c_p",
      "t_p",
      "c_s",
      "t_s",
      "t_clear",
      "inc_mean",
      "inc_sd",
      "sigma"
    ),
    name = c(
      "Time between infection and first positive test",
      "Latent upper Ct bound",
      "Peak Ct",
      "Time at peak Ct",
      "Switch Ct",
      "Time at switch Ct",
      "Time at clearance of infection",
      "Mean of the incubation period",
      "Standard deviation of the incubation period",
      "Observation standard deviation"
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
      "Assuming a normal error model"
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
      "Zero truncated normal"
    ),
    intercept_mean = c(5, 10, 0, 1.61, 0, 1.61, 2.3, 1.62, 0.42, 0),
    intercept_sd = c(5, 10, 1, 0.5, 1, 0.5, 0.5, 0.06, 0.07, 2)
  )

  priors[, `:=`(
    individual_variation_mean = individual_variation[1],
    individual_variation_sd = individual_variation[2]
  )]

  priors[variable %in% c("sigma", "c_int"), individual_variation_mean := NA]
  priors[variable %in% c("sigma", "c_int"), individual_variation_sd := NA]
  return(priors[])
}


#' Format data for use with stan
#'
#' @return A list as required by stan.
#' @family model
#' @importFrom data.table copy
#' @export
epict_to_stan <- function(obs,
                          ct_model = subject_design(~1, obs),
                          adjustment_model = test_design(~1, obs),
                          priors = epict::epict_priors(),
                          individual_variation = 0.2,
                          individual_correlation = 1,
                          censoring_threshold = 40, switch = TRUE,
                          onsets = TRUE, incubation_period = get_inc_period(),
                          likelihood = TRUE, output_loglik = FALSE) {
  obs <- data.table::copy(obs)
  obs <- obs[order(id)]
  obs[, obs := 1:.N]

  tests_per_id <- obs[, .(n = .N), by = "id"]$n

  stan_data <- list(
    N = obs[, .N],
    P = length(unique(obs$id)),
    id = obs[, id],
    tests_per_id = tests_per_id,
    cum_tests_per_id = cumsum(tests_per_id),
    day_rel = obs[, t],
    ct_value = obs$ct_value,
    ncensored = length(obs[uncensored == 0, obs]),
    censored = obs[uncensored == 0, obs],
    nuncensored = length(obs[uncensored == 1, obs]),
    uncensored = obs[uncensored == 1, obs],
    uncensored_by_test = obs[, uncensored],
    t_e = 0,
    c_int = censoring_threshold,
    c_lod = censoring_threshold,
    K = ifelse(switch, 5, 3),
    ind_var_sd = individual_variation,
    ind_var_m = as.numeric(individual_variation != 0),
    ind_corr = as.numeric(!is.na(individual_correlation) &&
      individual_variation != 0),
    lkj_prior = ifelse(
      is.na(individual_correlation), 0, individual_correlation
    ),
    lmean = incubation_period$inc_mean_p,
    lsd = incubation_period$inc_sd_p,
    likelihood = as.numeric(likelihood),
    output_loglik = as.numeric(output_loglik),
    preds = ncol(ct_model$design) - 1,
    preds_sd = ct_model$preds_sd,
    design = ct_model$design,
    switch = as.numeric(switch),
    adj_t_p = ct_model$params[["t_p"]],
    adj_t_s = min(ct_model$params[["t_s"]], as.numeric(switch)),
    adj_t_clear = ct_model$params[["t_clear"]],
    adj_c_p = ct_model$params[["c_p"]],
    adj_c_s = min(ct_model$params[["c_s"]], as.numeric(switch)),
    adj_inc_mean = ct_model$params[["inc_mean"]],
    adj_inc_sd = ct_model$params[["inc_sd"]],
    adj_ct = as.numeric(
      (ncol(adjustment_model$design) - 1) > 0
    ),
    ct_preds = ncol(adjustment_model$design) - 1,
    ct_preds_sd = adjustment_model$preds_sd,
    ct_design = adjustment_model$design
  )
  if (is.null(obs$onset_time) | !onsets) {
    stan_data <- c(stan_data, list(
      any_onsets = 0,
      onset_avail = rep(0, stan_data$P),
      onset_time = rep(0, stan_data$P),
      onset_window = rep(0, stan_data$P)
    ))
  } else {
    onset_data <- suppressWarnings(
      obs[,
        .(onset_time = min(onset_time, na.rm = TRUE), id),
        by = "id"
      ][
        is.infinite(onset_time), onset_time := NA
      ]
    )
    stan_data <- c(stan_data, list(
      any_onsets = 1,
      onset_avail = as.numeric(!is.na(onset_data$onset_time)),
      nonsets = sum(as.numeric(!is.na(onset_data$onset_time))),
      ids_with_onsets = onset_data[!is.na(onset_time), id],
      onset_time = ifelse(is.na(onset_data$onset_time),
                            0, onset_data$onset_time),
      onset_window = rep(1, stan_data$P)
    ))
  }

  return(stan_data)
}

#' Set up initial conditions for model
#'
#' @param data A list of data as produced by [epict_to_stan()].
#'
#' @return A function that when called returns a list of initial conditions
#' for the package stan models.
#'
#' @family model
#' @importFrom purrr map_dbl
#' @importFrom truncnorm rtruncnorm
#' @export
epict_inits <- function(data) {
  function() {
    inits <- list(
      t_inf = purrr::map_dbl(
        1:data$P,
        ~ truncnorm::rtruncnorm(
          1,
          a = max(-data$onset_time[.], 0),
          mean = max(-data$onset_time[.] + 5, 5), sd = 1
        )
      ),
      c_int = truncnorm::rtruncnorm(
        1,
        a = data$c_lod, mean = data$c_lod + 10, sd = 1
      ),
      c_p_int = rnorm(1, 0, 1),
      t_p_int = rnorm(1, 1.61, 0.5),
      t_clear_mean = rnorm(1, 2.3, 0.5),
      ind_var = abs(rnorm(data$K, 0, data$ind_var_sd * 0.1)),
      ind_eta = matrix(rnorm(data$P * data$K, 0, 0.1), nrow = data$K, ncol = data$P),
      sigma = truncnorm::rtruncnorm(1, a = 0, mean = 5, sd = 0.5)
    )

    if (data$switch > 0) {
      inits$c_s_int <- array(rnorm(1, 0, 1))
      inits$t_s_int <- array(rnorm(1, 1.61, 0.5))
    }

    if (data$preds > 0) {
      if (data$adj_t_p > 0) {
        inits$beta_t_p <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_t_s > 0) {
        inits$beta_t_s <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_t_clear > 0) {
        inits$beta_t_clear <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_c_p > 0) {
        inits$beta_c_p <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_c_s > 0) {
        inits$beta_c_s <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_inc_mean > 0) {
        inits$beta_inc_mean <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_inc_sd > 0) {
        inits$beta_inc_sd <- rnorm(data$preds, 0, 0.01)
      }
    }

    if (data$ct_preds > 0) {
      if (data$adj_ct > 0) {
        inits$beta_ct_shift <- rnorm(data$ct_preds, 0, 0.001)
        inits$beta_ct_scale <- rnorm(data$ct_preds, 0, 0.001)
      }
    }

    if (data$any_onsets == 1) {
      inits$inc_mean <- rnorm(1, data$lmean[1], data$lmean[2] * 0.1)
      inits$inc_sd <- truncnorm::rtruncnorm(
        1,
        a = 0, mean = data$lsd[1], sd = data$lsd[2] * 0.1
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
                        compile = TRUE, threads = FALSE, verbose = TRUE, ...) {
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
        stanc_options = list("O1"),
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