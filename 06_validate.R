# 06_validate.R ----------------------------------------------------------------
# Assemble output/quarterly_gsp.csv, run the automated constraint tests, the
# pseudo real-time test, the Queensland Treasury cross-check, the Method A vs
# Method B comparison and the sanity checks.

source("R/utils.R")
source("R/constraints.R")
source("R/baselines.R")
source("R/method_a.R")
source("R/method_b.R")
source("R/method_c.R")
source("R/validate.R")

inp <- load_inputs()
q_target <- read_proc("q_target")
est_a <- read_proc("est_method_a"); est_b <- read_proc("est_method_b"); est_bl <- read_proc("est_baselines")
a_diag <- read_proc("method_a_diag")
# Method C (mixed-frequency VAR) is optional: included when 05b_method_c.R has been run
has_c <- file.exists(file.path(DIR_PROC, "est_method_c.rds"))
est_c <- if (has_c) read_proc("est_method_c") else NULL

# 1. Output file ------------------------------------------------------------------
bench <- distinct(q_target, basis, date, benchmarked)
out <- bind_rows(
  est_a |> select(state, date, basis, method, estimate),
  est_b |> select(state, date, basis, method, estimate, lower_68, upper_68, lower_90, upper_90),
  if (has_c) est_c |> select(state, date, basis, method, estimate, lower_68, upper_68, lower_90, upper_90),
  est_bl |> select(state, date, basis, method, estimate)
) |>
  left_join(bench, by = c("basis", "date")) |>
  transmute(state, quarter = date, price_basis = basis, method, estimate,
            lower_68, upper_68, lower_90, upper_90, benchmarked) |>
  arrange(price_basis, method, state, quarter)
readr::write_csv(out, file.path(DIR_OUT, "quarterly_gsp.csv"))
msg("Wrote output/quarterly_gsp.csv:", nrow(out), "rows")
est <- out |> transmute(state, date = quarter, basis = price_basis, method, estimate, lower_90, upper_90, benchmarked)

# 2. Automated tests ---------------------------------------------------------------
Sys.setenv(GSP_ROOT = getwd())
tt <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE))
tests <- tibble::tibble(test = tt$test, passed = tt$passed, failed = tt$failed, error = tt$error)
msg("testthat:", sum(tests$passed), "expectations passed,", sum(tests$failed), "failed")

# 3. Pseudo real-time test ----------------------------------------------------------
test_years <- (inp$last_fy - a_diag$N_TEST + 1):inp$last_fy
msg("Pseudo real-time test, years", min(test_years), "to", max(test_years))
rt <- pseudo_realtime(inp, read_proc("oos_prelim"), test_years, a_diag$SEL_WINDOW)
if (has_c) {
  msg("Pseudo real-time test for Method C (", MFVAR_RT$n_iter, "iterations per vintage and price basis, in parallel)")
  # Method C's real-time run takes ~20 minutes; reuse the saved one when nothing it
  # depends on has changed (its inputs, settings and code).
  key_obj <- list(inp$a_state, inp$q_nat, select(inp$panel, state, basis, date, sfd, hours), MFVAR_DEFAULT, MFVAR_RT,
                  test_years, readLines("R/method_c.R"), readLines("R/validate.R"))
  tf <- tempfile(); saveRDS(key_obj, tf); key <- unname(tools::md5sum(tf)); unlink(tf)
  cache <- file.path(DIR_PROC, "rt_method_c_cache.rds")
  if (file.exists(cache) && identical(readRDS(cache)$key, key)) {
    msg("  reusing saved Method C real-time results (inputs unchanged)")
    rt_c <- readRDS(cache)$rt
  } else {
    rt_c <- pseudo_realtime_c(inp, test_years)
    saveRDS(list(key = key, rt = rt_c), cache)
  }
  rt <- bind_rows(rt, rt_c)
}
rt_sum <- summarise_rt(rt)
readr::write_csv(rt, file.path(DIR_OUT, "validation_realtime_errors.csv"))
readr::write_csv(rt_sum, file.path(DIR_OUT, "validation_realtime_summary.csv"))
cat("\nPseudo real-time mean absolute error of the implied annual GSP (% of actual), by state and method:\n")
print(as.data.frame(rt_sum |> select(basis, method, state, mae) |>
                      pivot_wider(names_from = method, values_from = mae) |>
                      mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

# 4. Queensland Treasury cross-check -------------------------------------------------
qld <- read_proc("qld_treasury")
cq <- compare_qld(est, qld$data)
readr::write_csv(cq$stats, file.path(DIR_OUT, "validation_qld_treasury.csv"))
cat("\nQueensland quarterly growth vs Queensland Treasury (", basename(qld$file), "):\n", sep = "")
print(as.data.frame(cq$stats |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

# context: how far has the annual picture moved since Queensland Treasury's last issue?
# (their financial-year sums vs today's ABS annual GSP, annual growth rates)
qld_annual <- qld$data |> mutate(fy = fy_of(date)) |> group_by(basis, fy) |>
  summarise(qt = sum(qld_treasury), n = n(), .groups = "drop") |> filter(n == 4) |>
  inner_join(filter(inp$a_state, state == "QLD"), by = c("basis", "fy")) |>
  group_by(basis) |> arrange(fy, .by_group = TRUE) |>
  mutate(g_qt = 100 * (qt / lag(qt) - 1), g_abs = 100 * (gsp / lag(gsp) - 1)) |>
  summarise(years = sum(!is.na(g_qt)), annual_growth_correlation = cor(g_qt, g_abs, use = "complete.obs"),
            annual_growth_mean_abs_diff = mean(abs(g_qt - g_abs), na.rm = TRUE), .groups = "drop")
cat("\nAnnual growth: Queensland Treasury (their vintage) vs current ABS annual GSP:\n"); print(as.data.frame(qld_annual), digits = 3)
readr::write_csv(qld_annual, file.path(DIR_OUT, "validation_qld_treasury_annual_context.csv"))

# 5. Method A vs Method B -------------------------------------------------------------
DIVERGE_PPT <- 1 # quarterly growth difference treated as material
g <- add_growth(est)
ab <- g |> filter(method %in% c("method_a", "method_b")) |>
  select(state, basis, date, method, growth, estimate, benchmarked) |>
  pivot_wider(names_from = method, values_from = c(growth, estimate)) |>
  left_join(est |> filter(method == "method_b") |> select(state, basis, date, lower_90, upper_90),
            by = c("state", "basis", "date")) |>
  mutate(diff_ppt = growth_method_a - growth_method_b,
         level_diff_pct = 100 * (estimate_method_a / estimate_method_b - 1),
         a_inside_b90 = estimate_method_a >= lower_90 & estimate_method_a <= upper_90,
         diverges = abs(diff_ppt) > DIVERGE_PPT)
ab_sum <- ab |> filter(!is.na(diff_ppt)) |> group_by(basis, state) |>
  summarise(correlation = cor(growth_method_a, growth_method_b), mean_abs_diff_ppt = mean(abs(diff_ppt)),
            max_abs_diff_ppt = max(abs(diff_ppt)), sd_growth_a = sd(growth_method_a), sd_growth_b = sd(growth_method_b),
            n_diverge = sum(diverges), pct_diverge = 100 * mean(diverges),
            mean_abs_level_diff_pct = mean(abs(level_diff_pct)), a_inside_b90_pct = 100 * mean(a_inside_b90), .groups = "drop")
ab_flags <- ab |> filter(diverges) |>
  transmute(basis, state, quarter = date, growth_method_a, growth_method_b, diff_ppt, benchmarked) |> arrange(basis, state, quarter)
readr::write_csv(ab_sum, file.path(DIR_OUT, "validation_a_vs_b_summary.csv"))
readr::write_csv(ab_flags, file.path(DIR_OUT, "validation_a_vs_b_divergent_quarters.csv"))
cat("\nMethod A vs Method B quarterly growth (divergence = more than", DIVERGE_PPT, "ppt apart):\n")
print(as.data.frame(ab_sum |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

# 6. Sanity checks ---------------------------------------------------------------------
GROWTH_FLAG <- 5 # per cent per quarter
flags <- g |> filter(!is.na(growth), abs(growth) > GROWTH_FLAG) |>
  transmute(method, basis, state, quarter = date, growth, benchmarked) |> arrange(method, basis, state, quarter)
readr::write_csv(flags, file.path(DIR_OUT, "sanity_growth_flags.csv"))
flag_sum <- flags |> count(method, basis, state, name = "n_flagged") |>
  pivot_wider(names_from = method, values_from = n_flagged, values_fill = 0)
cat("\nQuarters with growth beyond +/-", GROWTH_FLAG, "% (count by state):\n"); print(as.data.frame(flag_sum), row.names = FALSE)
# context for the flags: SFD growth in the same quarter
sfd_g <- inp$panel |> group_by(basis, state) |> arrange(date, .by_group = TRUE) |>
  mutate(sfd_growth = 100 * (sfd / lag(sfd) - 1)) |> ungroup() |> select(basis, state, date, sfd, sfd_growth)

gap_share <- est |> filter(method %in% c("method_a", "method_b")) |>
  left_join(sfd_g, by = c("basis", "state", "date")) |>
  mutate(gap_share_pct = 100 * (estimate - sfd) / estimate)
gap_sum <- gap_share |> group_by(method, basis, state) |>
  summarise(mean = mean(gap_share_pct), min = min(gap_share_pct), max = max(gap_share_pct),
            latest = last(gap_share_pct, order_by = date),
            max_abs_qtr_change_ppt = max(abs(diff(gap_share_pct[order(date)]))), .groups = "drop")
readr::write_csv(gap_sum, file.path(DIR_OUT, "sanity_gap_share.csv"))
cat("\nImplied gap (GSP - SFD) as % of GSP:\n")
print(as.data.frame(gap_sum |> mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

save_proc(list(tests = tests, rt = rt, rt_sum = rt_sum, qld = cq, qld_annual = qld_annual, qld_file = basename(qld$file),
               ab = ab, ab_sum = ab_sum, ab_flags = ab_flags, flags = flags, flag_sum = flag_sum,
               gap_share = gap_share, gap_sum = gap_sum, growth = g, sfd_g = sfd_g,
               DIVERGE_PPT = DIVERGE_PPT, GROWTH_FLAG = GROWTH_FLAG, test_years = test_years), "validation")
msg("06_validate.R done")
