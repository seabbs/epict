% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-tools.R
\name{epict_inference_opts}
\alias{epict_inference_opts}
\title{Format inference options for use with stan}
\usage{
epict_inference_opts(
  pp = FALSE,
  likelihood = TRUE,
  debug = FALSE,
  output_loglik = FALSE
)
}
\arguments{
\item{pp}{Logical, defaults to \code{FALSE}. Should posterior predictions be made
for observed data. Useful for evaluating the performance of the model.}

\item{likelihood}{Logical, defaults to \code{TRUE}. Should the likelihood be
included in the model}

\item{debug}{Logical, defaults to \code{FALSE}. Should within model debug
information be returned.}

\item{output_loglik}{Logical, defaults to \code{FALSE}. Should the
log-likelihood be output. Disabling this will speed up fitting
if evaluating the model fit is not required.}
}
\value{
A list as required by stan.
}
\description{
Format inference options for use with stan
}
\examples{
epict_inference_opts()
}
\seealso{
Functions used to design and setup models
\code{\link{adjustment_formula}()},
\code{\link{epict_formula_as_list}()},
\code{\link{epict_individual_priors_as_list}()},
\code{\link{epict_model_opts}()},
\code{\link{epict_obs_as_list}()},
\code{\link{epict_onset_obs_as_list}()},
\code{\link{epict_population_priors_as_list}()},
\code{\link{epict_posterior_as_prior}()},
\code{\link{piecewise_formula}()},
\code{\link{select_piecewise_parameters}()}
}
\author{
Sam Abbott
}
\concept{modeltools}
