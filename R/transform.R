transform_to_model <- function(draws) {
  draws <- data.table::as.data.table(draws)

  if (is.null(draws[["t_s"]])) {
    draws[, t_s := -Inf]
  }

  draws[
    ,
    `:=`(
      t_p = exp(t_p),
      t_s = exp(t_s),
      t_clear = exp(t_clear),
      c_thres = c_thres
    )
  ]

  if (is.null(draws[["c_s"]])) {
    draws[, c_s := c_thres]
  }else{
    draws[, c_s := c_thres * plogis(c_s)]
  }
  
  draws[
    ,
    c_p := c_s * plogis(c_p)
  ]
  return(draws[])
}

transform_to_natural <- function(draws) {
  draws <- transform_to_model(draws)
  draws[
    ,
    t_clear := t_p + t_s + t_clear
  ][
    ,
    t_s := t_p + t_s,
  ]
  return(draws[])
}

transform_ip_to_natural <- function(draws) {
  draws <- draws[, `:=`(
    nat_inc_mean = exp(inc_mean + (inc_sd^2) / 2),
    nat_inc_sd = sqrt((exp(inc_sd^2) - 1) * exp(2 * inc_mean + inc_sd^2))
  )]
  return(draws[])
}
