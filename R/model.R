#' Default model priors
#'
#'
#' @return A data frame summarising the model priors.
#'
#' @importFrom data.table data.table
#' @family model
#' @export
#' @examples
#' epict_priors()
epict_priors <- function() {
  data.table::data.table(
    variable = c(
      "eobs_lsd",
      "logmean_int",
      "logsd_int",
      "logmean_sd",
      "logsd_sd",
      "rd_eff_sd",
      "sqrt_phi"
    ),
    description = c(
      "Standard deviation for expected final observations",
      "Log mean intercept for reference date delay",
      "Log standard deviation for the reference date delay",
      "Standard deviation of scaled pooled logmean effects",
      "Standard deviation of scaled pooled logsd effects",
      "Standard deviation of scaled pooled report date effects",
      "One over the square of the reporting overdispersion"
    ),
    distribution = c(
      "Zero truncated normal",
      "Normal",
      "Zero truncated normal",
      "Zero truncated normal",
      "Zero truncated normal",
      "Zero truncated normal",
      "Zero truncated normal"
    ),
    mean = c(0, 1, 0.5, rep(0, 4)),
    sd = rep(1, 7)
  )
}

get_inc_period <- function(inc_mean = c(1.621, 0.0640),
                           inc_sd = c(0.418, 0.0691)) {
  list(
    inc_mean_p = inc_mean,
    inc_sd_p = inc_sd
  )
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
    c_0 = censoring_threshold,
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
    adj_t_lod = ct_model$params[["t_lod"]],
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
      c_0 = truncnorm::rtruncnorm(
        1,
        a = data$c_lod, mean = data$c_lod + 10, sd = 1
      ),
      c_p_mean = rnorm(1, 0, 1),
      t_p_mean = rnorm(1, 1.61, 0.5),
      t_lod_mean = rnorm(1, 2.3, 0.5),
      ind_var = abs(rnorm(data$K, 0, data$ind_var_sd * 0.1)),
      ind_eta = matrix(rnorm(data$P * data$K, 0, 0.1), nrow = data$K, ncol = data$P),
      sigma = truncnorm::rtruncnorm(1, a = 0, mean = 5, sd = 0.5)
    )

    if (data$switch > 0) {
      inits$c_s_mean <- array(rnorm(1, 0, 1))
      inits$t_s_mean <- array(rnorm(1, 1.61, 0.5))
    }

    if (data$preds > 0) {
      if (data$adj_t_p > 0) {
        inits$beta_t_p <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_t_s > 0) {
        inits$beta_t_s <- rnorm(data$preds, 0, 0.01)
      }
      if (data$adj_t_lod > 0) {
        inits$beta_t_lod <- rnorm(data$preds, 0, 0.01)
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