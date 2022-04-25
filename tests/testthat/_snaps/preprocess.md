# `epict_filter_ids()` works

    Code
      epict_filter_ids(obs)
    Output
         censored id t uncensored_tests days_with_uncensored
      1:     TRUE  1 0                2                    2
      2:    FALSE  1 0                2                    2
      3:    FALSE  1 1                2                    2

---

    Code
      epict_filter_ids(obs, invert = TRUE)
    Output
         censored id t uncensored_tests days_with_uncensored
      1:    FALSE  2 1                1                    1
      2:     TRUE  3 1                1                    1
      3:    FALSE  3 2                1                    1

---

    Code
      epict_filter_ids(obs, min_uncensored_tests = 3, min_days_with_uncensored = 1)
    Output
         censored id t uncensored_tests days_with_uncensored
      1:    FALSE  4 1                3                    1
      2:    FALSE  4 1                3                    1
      3:    FALSE  4 1                3                    1

# `epict_flag_spurious_obs()` works

    Code
      epict_flag_spurious_obs(obs, flag = FALSE)
    Message <simpleMessage>
      Spurious tests have been dropped
    Output
         t_rel_uncensored onset_t_rel_uncensored
      1:                0                      0
      2:                2                     NA
      3:               10                      2
      4:               60                     NA
      5:               30                      5

# `epict_make_time_rel_to_first_uncensored()` works

    Code
      epict_make_time_rel_to_first_uncensored(obs)
    Output
         id t onset_t censored t_first_uncensored t_rel_uncensored
      1:  1 2       2     TRUE                  3               -1
      2:  1 3       2    FALSE                  3                0
      3:  1 4       2    FALSE                  3                1
      4:  2 1       4     TRUE                  1                0
      5:  2 1       4    FALSE                  1                0
      6:  2 2       4    FALSE                  1                1
      7:  2 8       4    FALSE                  1                7
         onset_t_rel_uncensored
      1:                     -1
      2:                     -1
      3:                     -1
      4:                      3
      5:                      3
      6:                      3
      7:                      3

