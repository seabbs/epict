#' Check unproccessed observations meet the package specification
#'
#'Required raw data format:
#'  - id: An integer vector uniquely identifying eahc infection.
#'  - test_id: An integer vector uniquely identiying each test
#'  - ct_value: Numeric cycle threshold value.
#'  - test_date: Date of the test yielding a Ct value. (optional)
#'  - t: Time of test relative to a baseline date. Optional but required
#' if test_date is not present.
#'  - onset_date Date of onset for each infection (optional).
#'  NA if unavailable/asymptomatic
#'  - onset_t Time on onset relative to a baseline date (optional).
#'  - censored: Logical, indicating if the Ct has been censored.
#'
#' @param obs DESCRIPTION.
#' 
#' @param dates_available DESCRIPTION.
#' 
#' @param check_onset DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_check_raw_obs <- function(obs, dates_available = TRUE,
                                check_onset = FALSE) {
  cols <- c("id", "test_id", "ct_value", "censored")
  if (dates_available) {
    cols <- c(cols, "test_date")
  }else{
    cols <- c(cols, "t")
  }
  if (!check_onsets & dates_available) {
    cols <- c(cols, "onset_date")
  }else if (!check_onsets & !dates_available){
    cols <- c(cols, "onset_t")
  }
  return(obs[])
}
#' Make time relative to first test in the data
#'
#' FUNCTION_DESCRIPTION
#'
#' @param obs DESCRIPTION.
#' 
#' @param baseline_date DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_make_time_rel <- function(obs, baseline_date = min(obs$onset_date,
                                                         na.rm = TRUE)) {
  rel_obs <- data.table::copy(rel_obs)
  rel_obs[, t := as.numeric(test_date - baseline_date, units = "days")]
  if (!is.null(obs$onset_date)) {
    rel_obs[, onset_t := as.numeric(onset_date - baseline_date, units = "days")]
  }
  return(rel_obs[])
}
#' Drop Ct values with no data
#'
#' FUNCTION_DESCRIPTION
#'
#' @param obs DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_drop_na_ct <- function(obs) {
 fil_obs <- data.table::copy(obs)[!is.na(ct_value)]
 return(fil_obs[])
}
#' Filter infection IDs based on characteristics
#'
#' FUNCTION_DESCRIPTION
#'
#' @param obs DESCRIPTION.
#' 
#' @param min_uncensored_tests DESCRIPTION.
#' 
#' @param min_days_with_uncensored DESCRIPTION.
#'
#' @param invert Logical, defaults to `FALSE`. Should the filtering
#' requirements be inverted.
#' 
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_filter_ids <- function(obs, min_uncensored_tests = 2,
                             min_days_with_uncensored = 2, invert = FALSE) {
  fil_obs <- data.table::copy(obs)
  uncensored_by_id <- obs[
    censored == FALSE, .(uncensored_tests = .N), by = "id"
  ]
  fil_obs <- fil_obs[uncensored_by_id, on = "id"]
  if (!invert)  {
    fil_obs <- fil_obs[uncensored_tests >= min_uncensored_tests]
  }else{
    fil_obs <- fil_obs[uncensored_tests < min_uncensored_tests]
  }
  
  days_with_uncensored <- fil_obs[
    censored == FALSE, .(days_with_uncensored = .N), by = "id"
  ]
  fil_obs <- fil_obs[days_with_uncensored, on = "id"]
  fil_obs <- fil_obs[days_with_uncensored >= min_days_with_uncensored]
  if (!invert)  {
    fil_obs <- fil_obs[days_with_uncensored >= min_days_with_uncensored]
  }else{
    fil_obs <- fil_obs[days_with_uncensored < min_days_with_uncensored]
  }
  return(fil_obs[])
} 
#' Make time relative to first uncensored test per ID
#'
#' FUNCTION_DESCRIPTION
#'
#' @param obs DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_make_time_rel_to_first_uncensored <- function(obs) {
  rel_obs <- data.table::copy(rel_obs)
  first_uncensored <- rel_obs[
    censored == FALSE,
    .(t_first_uncensored = min(t, na.rm = TRUE)), by = "id"
  ]
  rel_obs <- rel_obs[first_uncensored, on = "id"]
  rel_obs[, t_rel_uncensored := t - t_first_uncensored]
  if (!is.null(obs$onset_date)) {
    rel_obs[, onset_t_rel_uncensored := onset_t - t_first_uncensored]
  }
  return(rel_obs[])
}
#' Flag potentially spurious IDs
#'
#' FUNCTION_DESCRIPTION
#'
#' @param max_t_rel_uncensored DESCRIPTION.
#' 
#' @param max_onset_t_rel_uncensored DESCRIPTION.
#' 
#' @param flag DESCRIPTION.
#' 
#' @param drop DESCRIPTION.
#' 
#' @param return_spurious DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_flag_spurious_ids <- function(max_t_rel_uncensored = 60,
                                    max_onset_t_rel_uncensored = 15,
                                    flag = TRUE, drop = TRUE,
                                    return_spurious = FALSE) {
  fil_obs <- data.table::copy(obs)
  spur_obs <- data.table::copy(obs)
  spur_obs[abs(t_rel_uncensored) > max_t_rel_uncensored]
  if (!is.null(obs$onset_t_rel_uncensored)) {
    onset_obs <- obs[abs(onset_t_rel_uncensored) >= max_onset_t_rel_uncensored]
    spur_obs <- rbind(spur_obs, onset_obs)
    spur_obs <- unique(spur_obs)
  }
  if (flag & nrow(spur_obs) > 0) {
    message("The following tests have been identified has spurious")
    print(spur_obs)
  }
  if (return_spurious) {
    return(spur_obs[])
  }
  if (drop & nrow(spur_obs) > 0) {
    message("Spurious tests have been dropped")
    fil_obs <- fil_obs[!spur_obs, on = "test_id"]
  }
  return(fil_obs[])
}
#' Clean up factor levels
#'
#' FUNCTION_DESCRIPTION
#'
#' @param vars DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
#' @importFrom forcats fct_drop
epict_clean_factors <- function(vars = c()) {
  clean_obs <- data.table::copy(obs)
  if (length(vars) > 0) {
    clean_obs <- clean_obs[,
     (vars) := lapply(.SD, as.factor), .SDcols = vars
    ]
    clean_obs <- clean_obs[,
     (vars) := lapply(.SD, forcats::fct_drop), .SDcols = vars
    ]
  }
  return(clean_obs[])
}
#' Check processed observations meet package specifications
#' 
#' Required processed data format:
#'  - id: An integer vector uniquely identifying eahc infection.
#'  - test_id: An integer vector uniquely identiying each test
#'  - ct_value: Numeric cycle threshold value.
#'  - t: Relative (to a baseline) time of the test yielding a Ct value.
#'  - t_rel_uncensored: Time of test relative to the first uncensored Ct value
#' for that id.
#'  - onset_t: Relative (to a baseline) time of onset for each infection
#'  - onset_t_rel_uncensored: Time of onset relative to the first uncensored Ct
#' value for that id. (optional). NA if unavailable/asymptomatic.
#'  - censored: Logical, indicating if the Ct has been censored.
#'
#' @param obs DESCRIPTION.
#' 
#' @param check_onsets DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @author Sam Abbott
#' @export
epict_check_obs <- function(obs, check_onsets = FALSE) {
  cols <- c("id", "test_id", "t", "t_rel_uncensored", "ct_value",
            "censored")
  if (check_onsets){
    cols <- c(cols, "onset_t", "onset_t_rel_uncensored")
  }
  return(obs[])
}
