% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{epict_check_raw_obs}
\alias{epict_check_raw_obs}
\title{Check unproccessed observations meet the package specification}
\usage{
epict_check_raw_obs(obs, dates_available = TRUE, check_onset = FALSE)
}
\arguments{
\item{obs}{A data.frame with the following variables:
\itemize{
\item \code{id}: An integer vector uniquely identifying each infection.
\item \code{test_id}: An integer vector uniquely identifying each test
\item \code{ct_value}: Numeric cycle threshold value.
\item \code{test_date}: Date of the test yielding a Ct value. (optional)
\item \code{t}: Time of test relative to a baseline date. Optional but required
if test_date is not present.
\item \code{onset_date}:e Date of onset for each infection (optional).
\code{NA} if unavailable/asymptomatic
\item \code{onset_t}: Time on onset relative to a baseline date (optional).
\item \code{censored}: Logical, indicating if the Ct has been censored.
}}

\item{dates_available}{Logical, defaults to \code{TRUE}. Are dates available
in the observed data. If \code{FALSE} it is assumed that relative times should be
available instead.}

\item{check_onset}{Logical, defaults to \code{FALSE}. Should observations be
checked for the presence of symptom onset data.}
}
\value{
Input observations are returned.
}
\description{
Checks the available raw data for required variables and returns informative
errors, warnings, and messages about the structure of the observations. It
can be used both on datasets referenced using dates and datasets referenced
using relative time.
}
\seealso{
Preprocessing functions
\code{\link{epict_check_obs}()},
\code{\link{epict_clean_factors}()},
\code{\link{epict_drop_na_ct}()},
\code{\link{epict_filter_ids}()},
\code{\link{epict_flag_spurious_obs}()},
\code{\link{epict_make_time_rel_to_first_uncensored}()},
\code{\link{epict_make_time_rel}()}
}
\author{
Sam Abbott
}
\concept{preprocess}
