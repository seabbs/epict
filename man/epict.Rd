% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epict.R
\name{epict}
\alias{epict}
\title{Model Cycle Thresholds}
\usage{
epict(
  obs,
  model = epict_model(),
  piecewise_formula = epict::piecewise_formula(~1, obs),
  adjustment_formula = epict::adjustment_formula(~1, obs),
  priors = epict::epict_priors(),
  model_opts = epict::epict_model_opts(),
  inference_opts = epict::epict_inference_opts(),
  convert_to_list = epict_convert_to_list,
  inits = epict_inits,
  ...
)
}
\arguments{
\item{obs}{A data.frame with the following variables:
\itemize{
\item \code{id}: An integer vector uniquely identifying eahc infection.
\item \code{test_id}: An integer vector uniquely identiying each test
\item \code{ct_value}: Numeric cycle threshold value.
\item \code{t}: Relative (to a baseline) time of the test yielding a Ct value.
\item \code{t_rel_uncensored}: Time of test relative to the first uncensored Ct value
for that id.
\item \code{onset_t}: Relative (to a baseline) time of onset for each infection
\item \code{onset_t_rel_uncensored}: Time of onset relative to the first uncensored
Ct value for that id. (optional). NA if unavailable/asymptomatic.
\item \code{censored}: Logical, indicating if the Ct has been censored.
}}

\item{model}{A compiled \code{cmdstanr} model object. By default this uses
\code{\link[=epict_model]{epict_model()}}.}

\item{piecewise_formula}{A list describing the piecewise linear cycle threshold
formula as described by \code{\link[=piecewise_formula]{piecewise_formula()}}.}

\item{adjustment_formula}{A list describing the cycle threshold linear
adjustment formula (shift and scale) as described by \code{\link[=adjustment_formula]{adjustment_formula()}}.}

\item{priors}{A data.table of population-level priors as produced
by \code{\link[=epict_priors]{epict_priors()}}.}

\item{model_opts}{A list of model options. See \code{\link[=epict_model_opts]{epict_model_opts()}}
for details}

\item{inference_opts}{A list of options to use for inference. See
\code{\link[=epict_inference_opts]{epict_inference_opts()}} for details.}

\item{convert_to_list}{Convert observations and model settings into
a list for use in stan. Defaults to using \code{\link[=epict_convert_to_list]{epict_convert_to_list()}}.}

\item{inits}{A function that combined with a list of data returns
sample initial conditions for use during model fitting based
on model priors. Defaults to \code{\link[=epict_inits]{epict_inits()}}}

\item{...}{Additional arguments passed to \link{cmdstanr}.}
}
\value{
A \link{cmdstanr} model fit containing posterior samples
}
\description{
Provides a user friendly interface around package functionality
to model cycle thresholds, and symptom onsets from observed preprocessed
data, a piecewise linear model, and a linear adjustment model.
}
\concept{epict}
