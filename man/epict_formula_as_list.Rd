% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-tools.R
\name{epict_formula_as_list}
\alias{epict_formula_as_list}
\title{Format formula data for use with stan}
\usage{
epict_formula_as_list(
  piecewise_formula,
  adjustment_formula,
  model_opts = epict::model_opts()
)
}
\arguments{
\item{piecewise_formula}{A list describing the piecewise linear cycle threshold
formula as described by \code{\link[=piecewise_formula]{piecewise_formula()}}.}

\item{adjustment_formula}{A list describing the cycle threshold linear
adjustment formula (shift and scale) as described by \code{\link[=adjustment_formula]{adjustment_formula()}}.}

\item{model_opts}{A list of Cycle threshold model options as returned
by \code{\link[=epict_ct_model_opts]{epict_ct_model_opts()}}.}
}
\value{
A list as required by stan.
}
\description{
The package stan code (accessed using \code{\link[=epict_model]{epict_model()}}) has
extensive more extensive documentation.
}
\seealso{
Functions used to design and setup models
\code{\link{adjustment_formula}()},
\code{\link{epict_individual_priors_as_list}()},
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
