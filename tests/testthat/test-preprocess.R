test_that("`epict_clean_factors()` works", {
  obs <- data.frame(
    m = c("fa", "aefwe", "efe"), 
    c = factor(c("fa", "asas", "asa"), levels = c("fa", "ph"))
  )
  ucobs <- epict_clean_factors(obs)
  expect_equal(obs, as.data.frame(ucobs))
  cobs <- epict_clean_factors(obs, vars = c("m", "c"))
  expect_true(is.factor(cobs$m))
  expect_true(is.factor(cobs$c))
  expect_equal(levels(ucobs$c), c("fa", "ph"))
  expect_equal(levels(cobs$c), "fa")
})

test_that("`epict_drop_na_ct()` works", {
  obs <- data.frame(ct_value = c(1, NA, 2))
  expect_equal(epict_drop_na_ct(obs)$ct_value, c(1,2))
})

test_that("`epict_filter_ids()` works", {
  obs <- data.frame(
    censored = c(TRUE, rep(FALSE, 3), TRUE, rep(FALSE, 4)),
    id = c(1, 1, 1, 2, 3, 3, 4, 4, 4),
    t = c(0, 0, 1, 1, 1, 2, 1, 1, 1)
  )
  expect_snapshot(epict_filter_ids(obs))
  expect_snapshot(epict_filter_ids(obs, invert = TRUE))
  expect_snapshot(
    epict_filter_ids(
      obs, min_uncensored_tests = 3, min_days_with_uncensored = 1)
  )
  expect_equal(nrow(epict_filter_ids(obs, min_days_with_uncensored = 3)), 0)
})

test_that("`epict_flag_spurious_obs()` works", {
  obs <- data.frame(t_rel_uncensored = c(0, 2, 10, 60, 100, 30))
  fil_obs <- suppressMessages(epict_flag_spurious_obs(obs, flag = FALSE))
  expect_equal(fil_obs$t_rel_uncensored, c(0, 2, 10, 60, 30))
  expect_equal(
    epict_flag_spurious_obs(fil_obs)$t_rel_uncensored, c(0, 2, 10, 60, 30)
  )
  obs$onset_t_rel_uncensored <- c(0, 15, 2, NA, 40, 5)
  fil_obs <- suppressMessages(epict_flag_spurious_obs(obs, flag = FALSE))
  expect_equal(fil_obs$onset_t_rel_uncensored, c(0, NA, 2, NA, 5))
  expect_equal(
    epict_flag_spurious_obs(fil_obs)$onset_t_rel_uncensored, c(0, NA, 2, NA, 5)
  )
  truth_obs <- data.table::data.table(
    t_rel_uncensored = c(100, 2), onset_t_rel_uncensored = c(40, 15)
  )
  expect_equal(epict_flag_spurious_obs(obs, return_spurious = TRUE),  truth_obs)
  expect_equal(
    suppressMessages(epict_flag_spurious_obs(obs, drop = FALSE)),
    data.table::as.data.table(obs)
  )
  expect_snapshot(epict_flag_spurious_obs(obs, flag = FALSE))
})

test_that("`epict_make_time_rel()` works", {
  obs <- data.frame(test_date = as.Date(c("2022-04-22", "2022-04-19")))
  tobs <- data.table::as.data.table(obs)[, t := c(3, 0)]
  expect_equal(epict_make_time_rel(obs), tobs)
  expect_equal(
    epict_make_time_rel(obs, baseline_date = as.Date("2022-04-01")),
    data.table::copy(tobs)[, t := c(21, 18)] 
  )
  obs$onset_date <- obs$test_date
  expect_equal(epict_make_time_rel(obs)$onset_t, c(3, 0))
})

test_that("`epict_make_time_rel_to_first_uncensored()` works", {
  obs <- data.frame(
    id = c(rep(1, 3), rep(2, 4)),
    t = c(2, 3, 4, 1, 1, 2, 8),
    onset_t = c(rep(2, 3), rep(4, 4)),
    censored = c(TRUE, FALSE, FALSE, TRUE, rep(FALSE, 3))
  )
  expect_snapshot(epict_make_time_rel_to_first_uncensored(obs))
})
