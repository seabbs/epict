% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-tools.R
\name{epict_individual_priors_as_list}
\alias{epict_individual_priors_as_list}
\title{Format individual-level priors for use in stan}
\usage{
epict_individual_priors_as_list(priors)
}
\arguments{
\item{priors}{A data.table of individual-level priors as produced
by \code{\link[=epict_priors]{epict_priors()}}.}
}
\value{
A list of named priors with the mean and standard deivation
for each.
}
\description{
Format individual-level priors for use in stan
}
\examples{
priors <- epict_priors()
epict_individual_priors_as_list(priors)
}
\seealso{
Functions used to design and setup models
\code{\link{adjustment_formula}()},
\code{\link{epict_formula_as_list}()},
\code{\link{epict_inference_opts}()},
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
