epict_theme <- function(flip = FALSE, legend_arg = FALSE) {
  epict_theme <- list(
    theme_bw(11),
    theme(strip.placement = "outside"),
    geom_vline(aes(xintercept = -Inf)),
    guides(color = guide_legend(override.aes = list(fill = NA)))
  )

  if (!flip) {
    append(epict_theme, geom_hline(aes(yintercept = -Inf)))
  }

  if (!legend_arg) {
    append(
      epict_theme,
      theme(legend.title = element_blank())
    )
  }

  return(epict_theme)
}

plot_obs <- function(obs, ct_traj, pp, traj_alpha = 0.05, onsets = TRUE,
                     reverse = TRUE, clod = min(obs$ct_value, na.rm = TRUE),
                     samples = 10, ...) {
  if (!missing(pp)) {
    obs <- cbind(
      obs[order(id)],
      data.table::copy(pp)[, c("t", "id", "obs") := NULL]
    )
  }

  plot <- ggplot2::ggplot(obs) +
    ggplot2::aes(x = t, y = ct_value, ...) +
    ggplot2::scale_colour_brewer(palette = "Dark2")

  if (!is.null(obs$onset_time) & onsets) {
    plot <- plot +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = onset_time), linetype = 2, alpha = 0.8
      )
  }

  if (!is.null(clod)) {
    plot <- plot +
      geom_hline(yintercept = clod, linetype = 3, alpha = 0.8)
  }

  if (!missing(pp)) {
    plot <- plot +
      geom_linerange(
        aes(ymin = lo90, ymax = hi90, y = NULL),
        size = 1.1, alpha = 0.15
      ) +
      geom_linerange(
        aes(ymin = lo60, ymax = hi60, y = NULL),
        size = 1.1, alpha = 0.15
      ) +
      geom_linerange(
        aes(ymin = lo30, ymax = hi30, y = NULL),
        size = 1.1, alpha = 0.15
      )
  }

  plot <- plot +
    geom_point(alpha = 0.6)

  if (!missing(ct_traj)) {
    plot <- plot +
      geom_line(
        data = ct_traj[iteration <= ceiling(samples / max(chain))],
        aes(
          y = value, x = time_since_first_pos,
          group = interaction(iteration, chain)
        ),
        colour = "black", alpha = traj_alpha
      )
  }

  if (reverse) {
    plot<- plot +
      ggplot2::scale_y_reverse()
  }

  plot <- plot +
    epict_theme() +
    theme(legend.position = "bottom") +
    labs(
      x = "Days since first positive test",
      y = "CT value"
    )
  return(plot)
}

plot_ct_pp <- function(pp, sum_pp, onsets = TRUE, clod = 40, alpha = 0.05,
                       reverse = TRUE, ...) {
  plot <- ggplot(pp) +
    aes(x = t, y = ct_value, group = interaction(.iteration, .chain), ...)

  if (!is.null(pp$onset_time) & onsets) {
    plot <- plot +
      geom_vline(aes(xintercept = onset_time), linetype = 2, alpha = 0.8)
  }

  if (!is.null(clod)) {
    plot <- plot +
      geom_hline(yintercept = clod, linetype = 3, alpha = 0.8)
  }

  if (!missing(sum_pp)) {
    plot <- plot +
      geom_ribbon(
        data = sum_pp,
        aes(ymin = lo90, ymax = hi90, y = NULL, group = NULL, col = NULL),
        alpha = 0.15
      ) +
      geom_ribbon(
        data = sum_pp,
        aes(ymin = lo60, ymax = hi60, y = NULL, group = NULL, col = NULL),
        alpha = 0.15
      ) +
      geom_ribbon(
        data = sum_pp,
        aes(ymin = lo30, ymax = hi30, y = NULL, group = NULL, col = NULL),
        alpha = 0.15
      )
  }

  plot <- plot +
    geom_line(alpha = alpha)

  if (reverse) {
    plot <- plot +
     ggplot2::scale_y_reverse()
  }

  plot <- plot +
    epict_theme() +
    labs(
      x = "Days since infection",
      y = "CT value"
    )
  return(plot)
}

plot_ip_pp <- function(pp, sum_pp, onsets = TRUE, alpha = 0.05, ...) {
  pp <- data.table::copy(pp)[, ct_value := value]
  plot <- plot_ct_pp(
    pp, sum_pp,
    onsets = TRUE, alpha = alpha, clod = NULL, ...
  ) +
    labs(y = "Probability of symptom onset")
  return(plot)
}

plot_density <- function(draws, ...) {
  plot <- ggplot(draws) +
    geom_density(aes(x = value, y = ..scaled.., ...), alpha = 0.2) +
    epict_theme() +
    labs(x = "", y = "Density")
  return(plot)
}

plot_ct_summary <- function(draws, time_range = seq(0, 60, by = 0.25),
                            samples = 100, by = c(), traj_alpha = 0.05,
                            reverse = TRUE, simulated_samples = 1000, ...) {
  pop_draws <- extract_ct_params(draws, by = by, mean = FALSE)

  pop_ct_draws <- pop_draws[.draw <= simulated_samples] |>
    transform_to_model() |>
    simulate_cts(time_range = time_range, obs_noise = FALSE)

  sum_cols <- c("value", "t", by)
  pop_ct_sum <- summarise_draws(
    pop_ct_draws[, value := ct_value][, ..sum_cols],
    by = setdiff(sum_cols, "value")
  )

  ct_pp_plot <- plot_ct_pp(
    pop_ct_draws[.draw <= samples], pop_ct_sum,
    alpha = traj_alpha, reverse = reverse, ...
  ) +
    guides(col = guide_none(), fill = guide_none())

  no_switch <- is.null(pop_draws[["t_s"]])

  param_pp <- pop_draws |>
    transform_to_natural() |>
    melt_draws(ids = c(".chain", ".iteration", ".draw", by))

  if (no_switch) {
    param_pp <- param_pp[!variable %in% c("t_s", "c_s")]
  }
  param_pp_plot <- param_pp |>
    update_variable_labels() |>
    plot_density(...) +
    ggplot2::facet_wrap(~variable, nrow = 2, scales = "free_x") +
    guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

  plot <- (param_pp_plot / ct_pp_plot) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  return(plot)
}

plot_ip_summary <- function(draws, time_range = seq(0, 20, by = 0.25),
                            samples = 100, by = c(), traj_alpha = 0.05,
                            simulated_samples = 1000, ...) {
  ip_draws <- extract_ip_params(draws, by = by)

  pop_ip_draws <- ip_draws[.draw <= simulated_samples] |>
    simulate_ips(time_range = time_range)

  sum_cols <- c("value", "t", by)
  pop_ip_sum <- summarise_draws(
    pop_ip_draws[, ..sum_cols],
    by = setdiff(sum_cols, "value")
  )

  ip_pp_plot <- plot_ip_pp(
    pop_ip_draws[.draw <= samples], pop_ip_sum,
    alpha = traj_alpha, ...
  ) +
    guides(col = guide_none(), fill = guide_none())

  param_pp_plot <- ip_draws |>
    transform_ip_to_natural() |>
    melt_draws(ids = c(".chain", ".iteration", ".draw", by)) |>
    update_variable_labels() |>
    plot_density(...) +
    ggplot2::facet_wrap(~variable, nrow = 2, scales = "free_x") +
    guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

  plot <- (param_pp_plot / ip_pp_plot) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  return(plot)
}

plot_summary <- function(draws, ct_time_range = seq(0, 60, by = 0.25),
                         ip_time_range = seq(0, 20, by = 0.25), samples = 100,
                         by = c(), traj_alpha = 0.05, simulated_samples = 1000,
                         reverse = TRUE, ...) {
  ct_pp <- plot_ct_summary(
    draws,
    time_range = ct_time_range, samples = samples, by = by,
    simulated_samples = simulated_samples, traj_alpha = traj_alpha,
    reverse = TRUE, ...
  )

  ip_pp <- plot_ip_summary(
    draws,
    time_range = ip_time_range, samples = samples, by = by,
    simulated_samples = simulated_samples, traj_alpha = traj_alpha,
    ...
  )

  parameter_pp <- (
    (ct_pp) | (
      ip_pp & guides(col = guide_none(), fill = guide_none())
    )
  ) +
    patchwork::plot_layout(
      widths = c(3, 2), guides = "collect"
    ) &
    theme(legend.position = "bottom")
  return(parameter_pp)
}

plot_effects <- function(effects, position = "identity", trans = "log", ...) {
  eff_plot <- ggplot(effects) +
    aes(y = variable, ...) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_linerange(
      aes(xmin = lo90, xmax = hi90),
      size = 3, alpha = 0.2,
      position = position
    ) +
    geom_linerange(
      aes(xmin = lo60, xmax = hi60),
      size = 3, alpha = 0.2,
      position = position
    ) +
    geom_linerange(
      aes(xmin = lo30, xmax = hi30),
      size = 3, alpha = 0.2,
      position = position
    ) +
    epict_theme() +
    theme(legend.position = "bottom") +
    scale_x_continuous(trans = trans) +
    labs(x = "Effect size", y = "Variable modified")
  return(eff_plot)
}

plot_effect_summary <- function(draws, ct_design, adjustment_design, variables,
                                variable_labels = function(dt, ...) {
                                  return(dt)
                                }, ...) {
  ct_plot <- draws |>
    summarise_effects(design = ct_design, variables = variables) |>
    update_predictor_labels() |>
    variable_labels(reverse = TRUE) |>
    plot_effects(...)

  adjustment_plot <- draws |>
    summarise_adjustment(
      design = adjustment_design
    ) |>
    update_predictor_labels() |>
    variable_labels(reverse = TRUE) |>
    plot_effects(
      trans = "identity", ...
    )

  summary <- (ct_plot + adjustment_plot) +
    patchwork::plot_layout(heights = c(4, 1))

  return(summary)
}
