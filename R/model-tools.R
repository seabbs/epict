params_avail_to_adjust <- function(params = "all") {
  choices <- c("t_p", "t_s", "t_lod", "c_p", "c_s", "inc_mean", "inc_sd")
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
#' @param data A list of stan observation data as produced by
#' [enw_obs_as_data_list()].
#'
#' @param ct_effects A list of fixed and random design matrices
#' defining the date of reference model. Defaults to [enw_formula()]
#' which is an intercept only model.
#'
#' @param adjustment_effects A list of fixed and random design matrices
#' defining the date of reports model. Defaults to [enw_formula()]
#' which is an intercept only model.
#'
#' @return A list as required by stan.
#' @family modeltools
#' @export
epict_formula_as_data_list <- function(data, ct_model, adjustment_model) {
  fdata <- list(
    npmfs = nrow(reference_effects$fixed$design),
    dpmfs = reference_effects$fixed$index,
    neffs = ncol(reference_effects$fixed$design) - 1,
    d_fixed = reference_effects$fixed$design,
    neff_sds = ncol(reference_effects$random$design) - 1,
    d_random = reference_effects$random$design
  )

  # map report date effects to groups and days
  report_date_eff_ind <- matrix(
    report_effects$fixed$index,
    ncol = data$g, nrow = data$t + data$dmax - 1
  )

  # Add report date data
  fdata <- c(fdata, list(
    rd = data$t + data$dmax - 1,
    urds = nrow(report_effects$fixed$design),
    rdlurd = report_date_eff_ind,
    nrd_effs = ncol(report_effects$fixed$design) - 1,
    rd_fixed = report_effects$fixed$design,
    nrd_eff_sds = ncol(report_effects$random$design) - 1,
    rd_random = report_effects$random$design
  ))
  return(fdata)
}

#' Format model options for use with stan
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
epict_opts_as_data_list <- function(pp = FALSE, likelihood = TRUE,
                                    debug = FALSE, output_loglik = FALSE) {

  data <- list(
    debug = as.numeric(debug),
    likelihood = as.numeric(likelihood),
    pp = as.numeric(pp),
    output_loglik = as.numeric(output_loglik)
  )
  return(data)
}

#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param priors DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @family modeltools
#' @importFrom data.table copy
#' @importFrom purrr map
#' @export
#' @examples
#' priors <- epict_priors()
#' epict_priors_as_data_list(priors)
epict_priors_as_data_list <- function(priors) {
  priors <- data.table::copy(priors)
  priors[, variable := paste0(variable, "_p")]
  priors <- priors[, .(variable, mean, sd)]
  priors <- split(priors, by = "variable", keep.by = FALSE)
  priors <- purrr::map(priors, ~ as.vector(t(.)))
  return(priors)
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
#' @param scale Numeric, defaults to 5. Amount to scale posterior standard
#' deviations by.
#'
#' @return A data frame of priors
#' @family modeltools
#' @importFrom data.table setDT
#' @export
epict_posterior_as_prior <- function(fit, priors = epict::epict_priors(),
                                     variables = c(), scale = 5) {
  posteriors <- fit$summary(variables)
  posteriors <- setDT(posteriors)[, sd := sd * scale]
  posteriors <- posteriors[, .(variable, mean, sd)]
  priors <- priors[!(variable %in% variables)]
  priors <- rbind(priors, posteriors, fill = TRUE)
  return(priors[])
}