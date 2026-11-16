# Validation helpers --------------------------------------------------------------

# Pseudo real-time test. For each test year t: annual GSP is known to t-1,
# quarterly data to June t. Every method is re-estimated on that vintage
# (Method A's model choice uses only out-of-sample errors observable at t-1;
# Method B's parameters are re-estimated). The implied annual total for year t
# is compared with the published figure.
# Method C in the real-time test: a shorter chain, posterior mean reconciled once,
# vintages run in parallel.
MFVAR_RT <- list(n_iter = 1500, burn = 500, thin = 4)
pseudo_realtime_c <- function(inp, test_years, cores = max(1, min(5, parallel::detectCores() - 1))) {
  actual <- inp$a_state |> rename(actual = gsp)
  res <- parallel::mclapply(test_years, function(t) {
    vin <- vintage_inputs(inp, t - 1)
    method_c(vin, utils::modifyList(MFVAR_DEFAULT, MFVAR_RT), reconcile_draws = FALSE, cores = 1)$est |>
      filter(fy_of(date) == t) |> group_by(basis, state) |> summarise(implied = sum(estimate), .groups = "drop") |>
      mutate(method = "method_c", fy = t)
  }, mc.cores = cores)
  bad <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("Method C real-time run failed: ", res[[which(bad)[1]]])
  bind_rows(res) |> left_join(actual, by = c("state", "basis", "fy")) |>
    mutate(model = NA_character_, err_pct = 100 * (implied / actual - 1))
}

pseudo_realtime <- function(inp, oos, test_years, sel_window) {
  actual <- inp$a_state |> rename(actual = gsp)
  purrr::map_dfr(test_years, function(t) {
    msg("  pseudo real-time vintage: annual data to FY", t - 1, "-> estimating FY", t)
    vin <- vintage_inputs(inp, t - 1)
    models <- select_models(oos, t - 1, sel_window) |> filter(rank == 1)
    a <- method_a(vin, models)$est
    b <- method_b(vin)$est
    all <- bind_rows(
      a |> transmute(state, basis, date, estimate, method = "method_a"),
      a |> transmute(state, basis, date, estimate = prelim, method = "method_a_prelim"),
      b |> transmute(state, basis, date, estimate, method = "method_b"),
      run_baselines(vin)
    )
    all |> filter(fy_of(date) == t) |>
      group_by(method, basis, state) |> summarise(implied = sum(estimate), .groups = "drop") |>
      mutate(fy = t) |>
      left_join(actual, by = c("state", "basis", "fy")) |>
      left_join(models |> transmute(basis, state, method = "method_a", model), by = c("basis", "state", "method"))
  }) |> mutate(err_pct = 100 * (implied / actual - 1))
}

summarise_rt <- function(rt) {
  by_state <- rt |> group_by(basis, method, state) |>
    summarise(mae = mean(abs(err_pct)), rmse = sqrt(mean(err_pct^2)), bias = mean(err_pct),
              max_abs = max(abs(err_pct)), n = n(), .groups = "drop")
  # all-state summaries: simple mean of state MAEs, and the error of the sum of states
  overall <- by_state |> group_by(basis, method) |>
    summarise(state = "Mean of 8", mae = mean(mae), rmse = sqrt(mean(rmse^2)), bias = mean(bias),
              max_abs = max(max_abs), n = sum(n), .groups = "drop")
  nat <- rt |> group_by(basis, method, fy) |> summarise(implied = sum(implied), actual = sum(actual), .groups = "drop") |>
    mutate(err_pct = 100 * (implied / actual - 1)) |> group_by(basis, method) |>
    summarise(state = "Sum of states", mae = mean(abs(err_pct)), rmse = sqrt(mean(err_pct^2)), bias = mean(err_pct),
              max_abs = max(abs(err_pct)), n = n(), .groups = "drop")
  bind_rows(by_state, overall, nat)
}

# Quarterly growth (per cent) by series
add_growth <- function(est) {
  est |> group_by(across(any_of(c("method", "basis", "state")))) |> arrange(date, .by_group = TRUE) |>
    mutate(growth = 100 * (estimate / lag(estimate) - 1)) |> ungroup()
}

# Turning points in the level: quarters where quarterly growth changes sign
turning_points <- function(g) which(sign(g[-1]) != sign(g[-length(g)]) & sign(g[-1]) != 0) + 1

# Compare Queensland with the Queensland Treasury series
compare_qld <- function(est, qld) {
  ours <- add_growth(filter(est, state == "QLD")) |> select(method, basis, date, growth)
  qt <- qld |> group_by(basis) |> arrange(date, .by_group = TRUE) |>
    mutate(qt_growth = 100 * (qld_treasury / lag(qld_treasury) - 1)) |> ungroup() |> select(basis, date, qt_growth)
  j <- inner_join(ours, qt, by = c("basis", "date")) |> filter(!is.na(growth), !is.na(qt_growth))
  stats <- j |> group_by(basis, method) |> arrange(date, .by_group = TRUE) |>
    summarise(
      from = min(date), to = max(date), n = n(),
      correlation = cor(growth, qt_growth),
      mean_abs_diff = mean(abs(growth - qt_growth)),
      sd_qt_growth = sd(qt_growth),
      same_sign_pct = 100 * mean(sign(growth) == sign(qt_growth)),
      same_accel_pct = 100 * mean(sign(diff(growth)) == sign(diff(qt_growth))),
      qt_turning_points = length(turning_points(qt_growth)),
      matched_within_1q = {
        tq <- turning_points(qt_growth); to <- turning_points(growth)
        sum(vapply(tq, function(i) any(abs(to - i) <= 1 & sign(growth[to]) == sign(qt_growth[i])), logical(1)))
      },
      our_turning_points = length(turning_points(growth)),
      .groups = "drop")
  list(stats = stats, series = j)
}
