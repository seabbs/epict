epict_summarise_coeff_pp <- function(fit, params, exponentiate = FALSE) {
  params <- epict::select_piecewise_parameters(params)
  params <- names(params[purrr::map_lgl(params, ~ . == 1)])

  draws <- fit$summary(
    variables = paste0("beta_", params)
  )
  draws <- data.table::as.data.table(draws)

  if (exponentiate) {
    beta_cols <- c("mean", "median", "sd", "mad", "q5", "q95")
    draws[, (beta_cols) := lapply(.SD, exp), .SDcols = beta_cols]
  }
  return(draws[])
}

epict_summarise_draws <- function(draws, by = c("id", "time")) {
  out <- data.table::as.data.table(draws)
  out <- out[,
    .(
      median = quantile(value, c(0.5), na.rm = TRUE),
      lo90 = quantile(value, c(0.05), na.rm = TRUE),
      lo60 = quantile(value, c(0.20), na.rm = TRUE),
      lo30 = quantile(value, c(0.35), na.rm = TRUE),
      hi30 = quantile(value, c(0.65), na.rm = TRUE),
      hi60 = quantile(value, c(0.80), na.rm = TRUE),
      hi90 = quantile(value, c(0.95), na.rm = TRUE)
    ),
    by = by
  ]
  return(out[])
}

epict_summarise_effects <- function(draws, design, variables,
                                    exponentiate = TRUE) {
  eff_draws <- epict::epict_extract_coeffs(
    draws,
    exponentiate = exponentiate, design = design, variables
  )

  by <- "variable"
  if (!missing(design)) {
    by <- c(by, "predictor")
  }
  eff_summary <- epict::epict_summarise_draws(eff_draws, by = by)
  return(eff_summary)
}

epict_summarise_adjustment <- function(draws, design) {
  eff_draws <- epict::epict_extract_coeffs(
    draws,
    exponentiate = FALSE, design = design,
    variables = c("ct_shift", "ct_scale")
  )

  eff_draws[variable %in% "ct_scale", value := exp(value)]

  by <- "variable"
  if (!missing(design)) {
    by <- c(by, "predictor")
  }
  eff_summary <- epict::epict_summarise_draws(eff_draws, by = by)
  return(eff_summary)
}

epict_summarise_pp <- function(fit, obs) {
  ct_pp <- epict::epict_extract_pp(fit, obs)
  ct_pp <- epict::epict_summarise_draws(
    ct_pp[, value := sim_ct],
    by = c("id", "t", "obs")
  )
  return(ct_pp)
}
