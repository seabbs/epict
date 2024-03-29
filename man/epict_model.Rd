% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{epict_model}
\alias{epict_model}
\title{Load and compile the Cycle thresholds model}
\usage{
epict_model(
  model,
  include,
  compile = TRUE,
  threads = FALSE,
  stanc_options = list("O1"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{model}{A character string indicating the path to the model.
If not supplied the package default model is used.}

\item{include}{A character string specifying the path to any stan
files to include in the model. If missing the package default is used.}

\item{compile}{Logical, defaults to \code{TRUE}. Should the model
be loaded and compiled using \code{\link[cmdstanr:cmdstan_model]{cmdstanr::cmdstan_model()}}.}

\item{threads}{Logical, defaults to \code{FALSE}. Should the model compile with
support for multi-thread support in chain. Note that this requires the use of
the \code{threads_per_chain} argument when model fitting using \code{\link[=epict]{epict()}} and is
not currently supported in model.}

\item{stanc_options}{A list of options to pass to the \code{stanc_options} of
\code{\link[cmdstanr:cmdstan_model]{cmdstanr::cmdstan_model()}}. By default "01" is passed which specifies simple
optimisations should be done by the prior to compilation.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose
messages be shown.}

\item{...}{Additional arguments passed to \code{\link[cmdstanr:cmdstan_model]{cmdstanr::cmdstan_model()}}.}
}
\value{
A \code{cmdstanr} model.
}
\description{
Load and compile the Cycle thresholds model
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
mod <- epict_model()
mod
\dontshow{\}) # examplesIf}
}
\seealso{
Functions used for modelling
\code{\link{epict_convert_to_list}()},
\code{\link{epict_inits}()},
\code{\link{epict_priors}()}
}
\concept{model}
