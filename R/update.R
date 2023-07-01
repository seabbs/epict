#' Update Variable Labels
#'
#' This function is used to update the labels of variables in a given dataset.
#' The labels are updated according to a specified set of new labels.
#' If 'reverse' is set to TRUE, the order of variables and their corresponding
#' labels are reversed before the update.
#'
#' @param draws A dataset (data.frame or data.table) containing the variables
#' to be relabelled.
#' @param reverse A logical value indicating whether to reverse the order of
#' variables and their corresponding labels before updating (Default is FALSE).
#' @param variables A data.frame providing the mapping from old variable names
#' to new names. By default, it uses
#' `epict::epict_priors(include_descriptive = TRUE)`. The data.frame should
#' contain 'variable' and 'name' columns where 'variable' is the old variable
#' name and 'name' is the new variable name.
#'
#' @return A data.table with updated variable labels according to the input
#' mapping. The returned data.table contains only the rows where variables were
#' found and updated.
#' @export
#' @importFrom data.table as.data.table
#' @examples
#' data <- data.frame(variable = c("v1", "v2", "v3"), value = c(1, 2, 3))
#' variables_mapping <- data.frame(
#'  variable= c("v1", "v2", "v3"),
#'  name = c("Variable 1", "Variable 2", "Variable 3")
#' )
#'
#' updated_data <- update_variable_labels(data, FALSE, variables_mapping)
#' updated_data
update_variable_labels <- function(
  draws, reverse = FALSE,
  variables = epict::epict_priors(include_descriptive = TRUE)
) {
  draws <- as.data.table(draws)
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
