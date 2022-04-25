#' Check unproccessed observations meet the package specification
#'
#' Checks the available raw data for required variables and returns informative
#' errors, warnings, and messages about the structure of the observations. It
#' can be used both on datasets referenced using dates and datasets referenced
#' using relative time.
#'
#' @param obs A data.frame with the following variables:
#'  - `id`: An integer vector uniquely identifying each infection.
#'  - `test_id`: An integer vector uniquely identiying each test
#'  - `ct_value`: Numeric cycle threshold value.
#'  - `test_date`: Date of the test yielding a Ct value. (optional)
#'  - `t`: Time of test relative to a baseline date. Optional but required
#' if test_date is not present.
#'  - `onset_date`:e Date of onset for each infection (optional).
#'  `NA` if unavailable/asymptomatic
#'  - `onset_t`: Time on onset relative to a baseline date (optional).
#'  - `censored`: Logical, indicating if the Ct has been censored. 
#' 
#' @param dates_available Logical, defaults to `TRUE`. Are dates available
#' in the observed data. If `FALSE` it is assumed that relative times should be
#' available instead.
#' 
#' @param check_onset Logical, defaults to `FALSE`. Should observations be
#' checked for the presence of symptom onset data.
#'
#' @return Input observations are returned.
#' @importFrom data.table as.data.table
#' @family preprocess
#' @author Sam Abbott
#' @export
epict_check_raw_obs <- function(obs, dates_available = TRUE,
                                check_onset = FALSE) {
  obs <- data.table::as.data.table(obs)
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
#' Make time relative to first date of the test.
#'
#' @param obs A `data.frame` with the following variables:
#'  - `test_date`: Date of the test yielding a Ct value.
#'  - `onset_date`: Date of onset for each infection (optional).
#' 
#' @param baseline_date A date to use as the baseline date for all other
#' times in the dataset. By default this is the minimum test date in the 
#' dataset.
#'
#' @return A `data.table` with a `t` variable containing the numeric time in
#' days since the `baseline_date`.
#' @importFrom data.table as.data.table
#' @family preprocess
#' @author Sam Abbott
#' @export
#' @examples
#' obs <- data.frame(test_date = as.Date(c("2022-04-22", "2022-04-19")))
#' 
#' # Using default baseline date
#' epict_make_time_rel(obs)
#' 
#' # Using user specified baseline date
#' epict_make_time_rel(obs, baseline_date = as.Date("2022-04-01"))
#' 
#' # When onsets are present
#' obs$onset_date <- obs$test_date
#' epict_make_time_rel(obs)
epict_make_time_rel <- function(obs, baseline_date = min(obs$test_date,
                                                         na.rm = TRUE)) {
  rel_obs <- data.table::as.data.table(obs)
  rel_obs[, t := as.numeric(test_date - baseline_date, units = "days")]
  if (!is.null(obs$onset_date)) {
    rel_obs[, onset_t := as.numeric(onset_date - baseline_date, units = "days")]
  }
  return(rel_obs[])
}
#' Drop Ct values with no data
#'
#' @param obs A `data.frame` with a numeric `ct_value` variable.
#'
#' @return The input `data.table` with `NA` values dropped.
#' @family preprocess
#' @author Sam Abbott
#' @export
#' @examples
#' obs <- data.frame(ct_value = c(1, NA, 2))
#' epict_drop_na_ct(obs)
epict_drop_na_ct <- function(obs) {
 fil_obs <- data.table::as.data.table(obs)[!is.na(ct_value)]
 return(fil_obs[])
}
#' Make time relative to first uncensored test per ID
#'
#' @param obs A data.frame with the following variables:
#'  - `id`: An integer vector uniquely identifying each infection.
#'  - `t`: Time of test relative to a baseline date.
#'  - `onset_t`: Time on onset relative to a baseline date (optional).
#'  - `censored`: Logical, indicating if the Ct has been censored. 
#'
#' @return A `data.table` with times relative to the first uncensored test per
#' ID.
#' @importFrom data.table as.data.table
#' @family preprocess
#' @author Sam Abbott
#' @export
#' @examples
#' obs <- data.frame(
#'  id = c(rep(1, 3), rep(2, 4)), t = c(2, 3, 4, 1, 1, 2, 8), 
#'  onset_t = c(rep(2, 3), rep(4, 4)),
#'  censored = c(TRUE, FALSE, FALSE, TRUE, rep(FALSE, 3))
#' )
#' 
#' epict_make_time_rel_to_first_uncensored(obs)
epict_make_time_rel_to_first_uncensored <- function(obs) {
  rel_obs <- data.table::as.data.table(obs)
  first_uncensored <- rel_obs[
    censored == FALSE,
    .(t_first_uncensored = min(t, na.rm = TRUE)), by = "id"
  ]
  rel_obs <- rel_obs[first_uncensored, on = "id"]
  rel_obs[, t_rel_uncensored := t - t_first_uncensored]
  if (!is.null(obs$onset_t)) {
    rel_obs[, onset_t_rel_uncensored := onset_t - t_first_uncensored]
  }
  return(rel_obs[])
}
#' Flag potentially spurious observations
#'
#' @param obs A data.frame with the following variables:
#'  - `t_rel_uncensored`: Time of test relative to the first uncensored Ct value
#' for that id.
#'  - `onset_t_rel_uncensored`: Time of onset relative to the first uncensored
#'  Ct value for that id. (optional). NA if unavailable/asymptomatic.
#'  
#' @param max_t_rel_uncensored Numeric defaults to 60 days. Flags the maximum
#' absolute relative time that is considered plausible for tests to be
#' spaced by per ID.
#' 
#' @param max_onset_t_rel_uncensored  Numeric defaults to 15 days. Flags the
#' maximum absolute relative time that is considered plausible for symptom
#' onset from the first uncensored test per ID.
#' 
#' @param flag Logical, defaults to `TRUE`. Should spurious tests be flagged
#' 
#' @param drop Logical, defaults to `TRUE`. Should spurious tests
#' and onsets be dropped (set to `NA` for onsets and filtered out for tests).
#' 
#' @param return_spurious Logical, defaults to `FALSE`. Rather than returning
#' observations should observations flagged as spurious be returned.
#'
#' @return A `data.table` of observations
#' @importFrom data.table as.data.table
#' @family preprocess
#' @author Sam Abbott
#' @export
#' @examples
#' obs <- data.frame(t_rel_uncensored = c(0, 2, 10, 60, 100, 30))
#' 
#' # Run with defaults
#' fil_obs <- epict_flag_spurious_obs(obs)
#' fil_obs
#' epict_flag_spurious_obs(fil_obs)
#' 
#' # Add onsets and repeat
#' obs$onset_t_rel_uncensored <- c(0, 15, 2, NA, 40, 5)
#' fil_obs <- epict_flag_spurious_obs(obs)
#' fil_obs
#' epict_flag_spurious_obs(fil_obs)
#' 
#' # Return spurious observations
#' epict_flag_spurious_obs(obs, return_spurious = TRUE)
#' 
#' # Flag spurious results but don't drop
#' epict_flag_spurious_obs(obs, drop = FALSE)
#' 
#' # Siltently drop
#' epict_flag_spurious_obs(obs, flag = FALSE)
epict_flag_spurious_obs <- function(obs, max_t_rel_uncensored = 60,
                                    max_onset_t_rel_uncensored = 15,
                                    flag = TRUE, drop = TRUE,
                                    return_spurious = FALSE) {
  fil_obs <- data.table::as.data.table(obs)
  spur_tests <- data.table::as.data.table(obs)
  spur_tests <- spur_tests[abs(t_rel_uncensored) > max_t_rel_uncensored]
  if (!is.null(obs$onset_t_rel_uncensored)) {
    onset_obs <- fil_obs[
      abs(onset_t_rel_uncensored) >= max_onset_t_rel_uncensored
    ]
    spur_obs <- rbind(spur_tests, onset_obs)
    spur_obs <- unique(spur_obs)
  }else{
    spur_obs <- spur_tests
  }
  if (flag & nrow(spur_obs) > 0 & !return_spurious) {
    message("The following tests have been identified as having spurious data")
    print(spur_obs)
  }
  if (return_spurious) {
    return(spur_obs[])
  }
  if (drop & nrow(spur_obs) > 0) {
    message("Spurious tests have been dropped")
    fil_obs <- fil_obs[abs(t_rel_uncensored) <= max_t_rel_uncensored]
    if (!is.null(obs$onset_t_rel_uncensored)) {
      fil_obs[
        abs(onset_t_rel_uncensored) >= max_onset_t_rel_uncensored,onset_t_rel_uncensored := NA
      ]
    }
  }
  return(fil_obs[])
}
#' Filter infection IDs based on characteristics
#'
#'
#' @param obs A data.frame with the following variables:
#'  - `id`: An integer vector uniquely identifying each infection.
#'  - `t`: Time of test relative to a baseline date.
#'  - `censored`: Logical, indicating if the Ct has been censored. 
#' 
#' @param min_uncensored_tests Numeric defaults to 2. The minimum number of
#' uncensored tests an ID must have in order to be included in the processed
#' dataset.
#' 
#' @param min_days_with_uncensored Numeric defaults to 2. The minimum number of
#' days tests per ID must span in order to be included in the processed dataset.
#'
#' @param invert Logical, defaults to `FALSE`. Should the filtering
#' requirements be inverted. This switches from exluding to including and
#' can be used for debugging and data exploration.
#' 
#' @return A `data.table` with IDs failing/meeting the specified criteria
#' removed.
#' @family preprocess
#' @author Sam Abbott
#' @importFrom data.table as.data.table uniqueN
#' @export
#' @examples
#' obs <- data.frame(
#'  censored = c(TRUE, rep(FALSE, 3), TRUE, rep(FALSE, 4)),
#'  id = c(1, 1, 1, 2, 3, 3, 4, 4, 4), t = c(0, 0, 1, 1, 1, 2, 1, 1, 1)
#' )
#' 
#' # Use defaults
#' epict_filter_ids(obs)
#' 
#' # Invert with defaults
#' epict_filter_ids(obs, invert = TRUE)
#' 
#' # Remove IDs with less than 3 tests
#' epict_filter_ids(obs, min_uncensored_tests = 3, min_days_with_uncensored = 1)
#' 
#' #' # Remove IDs with less than 3 tests - this will be all available data
#' epict_filter_ids(obs, min_days_with_uncensored = 3)
epict_filter_ids <- function(obs, min_uncensored_tests = 2,
                             min_days_with_uncensored = 2, invert = FALSE) {
  obs <- data.table::as.data.table(obs)
  uncensored_by_id <- obs[
    censored == FALSE, .(uncensored_tests = .N), by = "id"
  ]
  fil_obs <- obs[uncensored_by_id, on = "id"]
  if (!invert)  {
    fil_obs <- fil_obs[uncensored_tests >= min_uncensored_tests]
  }else{
    fil_obs <- fil_obs[uncensored_tests < min_uncensored_tests]
  }
  days_with_uncensored <- fil_obs[
    censored == FALSE, .(days_with_uncensored = uniqueN(t)), by = "id"
  ]
  fil_obs <- fil_obs[days_with_uncensored, on = "id"]
  if (!invert)  {
    fil_obs <- fil_obs[days_with_uncensored >= min_days_with_uncensored]
  }else{
    fil_obs <- fil_obs[days_with_uncensored < min_days_with_uncensored]
  }
  return(fil_obs[])
} 
#' Clean up factor levels
#' 
#' This utility function drops empty factor levels from target factors as well
#' as converting them to factors from character variables.
#'
#' @param vars A character vector of variables. Defaults to empty.
#'
#' @return A `data.table` with empty factors dropped and character vectors
#' transformed to factors as specified.
#' @family preprocess
#' @author Sam Abbott
#' @export
#' @importFrom data.table as.data.table
#' @importFrom forcats fct_drop
#' @examples
#' obs <- data.frame(
#'  m = c("fa", "aefwe", "efe"),
#'  c = factor(c("fa", "asas", "asa"), levels = c("fa"))
#')
#' summary(obs)
#' 
#' # Default without specifying variables
#' summary(epict_clean_factors(obs))
#' 
#' # Specify variables
#' summary(epict_clean_factors(obs, vars = c("m", "c")))
epict_clean_factors <- function(obs, vars = c()) {
  clean_obs <- data.table::as.data.table(obs)
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
#' Check proccessed observations meet the package specification
#'
#' Checks the available proccessed observastionsfor required variables and
#' returns informative errors, warnings, and messages about the structure of
#' the observations.
#'
#' @param obs A data.frame with the following variables:
#'  - `id`: An integer vector uniquely identifying eahc infection.
#'  - `test_id`: An integer vector uniquely identiying each test
#'  - `ct_value`: Numeric cycle threshold value.
#'  - `t`: Relative (to a baseline) time of the test yielding a Ct value.
#'  - `t_rel_uncensored`: Time of test relative to the first uncensored Ct value
#' for that id.
#'  - `onset_t`: Relative (to a baseline) time of onset for each infection
#'  - `onset_t_rel_uncensored`: Time of onset relative to the first uncensored
#' Ct value for that id. (optional). NA if unavailable/asymptomatic.
#'  - `censored`: Logical, indicating if the Ct has been censored.
#' 
#' @return A `data.table` of observations ready for use in [epict()] and other 
#' package functions.
#' @inheritParams epict_check_raw_obs
#' @family preprocess
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
