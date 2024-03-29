epict_extract_subjects <- function(dt) {
  subjects <- data.table::as.data.table(dt)
  subjects <- subjects[, .SD[1, ], by = "id"]
  return(subjects[])
}

epic_extract_draws <- function(fit, params = NULL, format = "df") {
  draws <- fit$draws(format = format, variables = params)
  draws <- data.table::as.data.table(draws)
  return(draws[])
}

epict_extract_params <- function(draws, params, by) {
  if (!missing(by)) {
    params <- c(params, by)
  }
  cols <- intersect(colnames(draws), params)
  cols <- c(".iteration", ".draw", ".chain", cols)
  draws <- draws[, ..cols]
  return(draws[])
}

epict_extract_ct_params <- function(draws, params = c(
                                "c_thres", "c_p_int",
                                "c_s_int[1]", "t_p_int",
                                "t_s_int[1]", "t_clear_int"
                              ),
                              intercept = TRUE, by) {
  if (!mean) {
    params <- stringr::str_remove(params, "_int")
    params <- stringr::str_remove(params, "\\[1\\]")
  }
  draws <- epict::epict_extract_params(draws, params = params, by)
  colnames(draws) <- stringr::str_remove(colnames(draws), "_int")
  colnames(draws) <- stringr::str_remove(colnames(draws), "\\[1\\]")
  return(draws[])
}

epict_extract_ip_params <- function(draws, params = c(
                                "inc_mean[1]", "inc_mean",
                                "inc_sd[1]", "inc_sd"
                              ),
                              by) {
  draws <- epict::epict_extract_params(draws, params = params, by)
  colnames(draws) <- purrr::map_chr(
    colnames(draws), ~ stringr::str_split(., "\\[[0-9]\\]")[[1]][1]
  )
  return(draws[])
}

epict_extract_param_draws <- function(draws) {
  draws <- cbind(
    epict::epict_extract_ct_params(draws),
    epict::epict_extract_ip_params(draws)[, .(inc_mean, inc_sd)]
  )
  return(draws[])
}

epict_extract_coeffs <- function(draws, exponentiate = FALSE,
                                 design, variables) {
  beta_cols <- grep("beta_", colnames(draws), value = TRUE)
  cols <- c(".iteration", ".draw", ".chain", beta_cols)
  draws <- data.table::as.data.table(draws)[, ..cols]

  if (exponentiate) {
    draws[, (beta_cols) := lapply(.SD, exp), .SDcols = beta_cols]
  }

  setnames(draws, beta_cols, gsub("beta_", "", beta_cols))

  draws <- epict::epict_melt_draws(draws)

  draws <- draws[
    ,
    coeff := as.numeric(stringr::str_extract(variable, "[0-9]"))
  ][
    ,
    variable := purrr::map_chr(
      variable, ~ stringr::str_split(., "\\[[0-9]\\]")[[1]][1]
    )
  ]

  if (!missing(variables)) {
    draws <- draws[variable %in% variables]
  }

  if (!missing(design)) {
    design <- data.table::data.table(predictor = colnames(design))
    design <- design[!predictor %in% "(Intercept)"]
    design[, coeff := 1:.N]
    draws <- draws[design, on = "coeff"]
  }
  return(draws[])
}

epict_melt_draws <- function(draws, ids = c(".chain", ".iteration", ".draw")) {
  data.table::melt(draws, id.vars = ids)
}

epict_extract_ct_trajectories <- function(fit, variable = "ct",
                                          infection_time = TRUE) {
  dt_draws <- epict::epict_extract_draws(
    fit, params = variable, format = "array"
  )

  obs_out <- dt_draws[,
    c("id", "time") := tstrsplit(variable, ",")
  ][
    ,
    id := stringr::str_remove(id, paste0("ct", "\\["))
  ][
    ,
    time := stringr::str_remove(time, "\\]")
  ][
    ,
    time := as.numeric(time)
  ][
    ,
    id := factor(id)
  ][
    ,
    c("time", "iteration", "chain", "id", "value")
  ][
    order(id, time)
  ]

  if (infection_time) {
    inf_time_draws <- extract_draws(fit, params = "t_inf", format = "array")[
      ,
      id := str_remove(variable, "t_inf\\[")
    ][
      ,
      id := str_remove(id, "\\]")
    ][
      ,
      .(id, inf_time = value, iteration, chain)
    ]

    obs_out <- obs_out[inf_time_draws, on = c("id", "iteration", "chain")]
  }
  obs_out[, time_since_first_pos := time - inf_time]
  cols <- c("chain", "iteration")
  obs_out[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  return(obs_out[])
}

extract_posterior_predictions <- function(fit, obs) {
  dt_draws <- extract_draws(fit, "sim_ct", format = "array")

  simulated_cts <- dt_draws[
    ,
    obs := str_remove(variable, "sim_ct\\[")
  ][
    ,
    obs := str_remove(obs, "\\]")
  ][
    ,
    .(
      obs = as.numeric(obs), sim_ct = value, iteration = as.numeric(iteration),
      chain = as.numeric(chain)
    )
  ][order(obs)]

  if (!missing(obs)) {
    simulated_cts <- merge(
      obs[order(id), obs := 1:.N], simulated_cts,
      by = "obs"
    )
  }
  return(simulated_cts[])
}
