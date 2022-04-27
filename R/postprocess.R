#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param draws DESCRIPTION.
#' @param reverse DESCRIPTION.
#' @param variables DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
#' @family postprocess
#' @author Sam Abbott
#' @examples
#' # ADD_EXAMPLES_HERE
update_variable_labels <- function(
  draws, reverse = FALSE, 
  variables = epict::epict_priors(include_descriptive = TRUE)
) {
  draws <- data.table::as.data.table(draws)
  params <- variables$variable

  clean_params <- variables$name

  if (reverse) {
    params <- rev(params)
    clean_params <- rev(clean_params)
  }

  draws <- draws[variable %in% params][
    ,
    variable := factor(
      variable,
      levels = params,
      labels = clean_params
    )
  ]
  return(draws[])
}
