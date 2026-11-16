# 04b_bottom_up_test.R ------------------------------------------------------------
# Does the bottom-up industry indicator help? Run after 04_method_a.R.
#  1. Fixed preliminary models: each indicator set with the lagged-share indicator
#     (prodgap) against the same set with the bottom-up indicator (bugap).
#  2. Method A with model selection restricted to the old indicator sets, against
#     selection from all sets (the production setting).
#  3. Method B with prodgap against Method B with bugap.
# 2 and 3 use the pseudo real-time test of 06_validate.R (same years; every
# vintage rebuilds the bottom-up indicator from the annual data it would have had).

source("R/utils.R")
source("R/constraints.R")
source("R/method_a.R")
source("R/method_b.R")
source("R/validate.R")

inp <- load_inputs()
a_diag <- read_proc("method_a_diag")
oos <- read_proc("oos_prelim")
test_years <- (inp$last_fy - a_diag$N_TEST + 1):inp$last_fy
spec_of <- function(m) sub("^[^:]+:", "", m)
td_of <- function(m) sub(":.*$", "", m)

# 1. Fixed models ----------------------------------------------------------------
pairs <- tibble::tribble(~old, ~new,
  "s1_prodgap", "s7_bugap", "s2_prodgap_nx", "s8_bugap_nx",
  "s5_prodgap_nx_hrs", "s9_bugap_nx_hrs", "s6_prod_sfd_nx", "s10_bu_sfd_nx")
fixed <- oos |> filter(fy %in% test_years) |> mutate(spec = spec_of(model), td = td_of(model)) |>
  group_by(basis, td, spec, state) |> summarise(mae = mean(abs(err_pct)), .groups = "drop")
fixed_cmp <- pairs |> inner_join(fixed, by = c("old" = "spec"), relationship = "many-to-many") |> rename(mae_old = mae) |>
  inner_join(fixed |> rename(mae_new = mae), by = c("new" = "spec", "basis", "td", "state")) |>
  group_by(basis, td, old, new) |>
  summarise(mae_old = mean(mae_old), mae_new = mean(mae_new), states_better = sum(mae_new < mae_old), .groups = "drop") |>
  mutate(change_pct = 100 * (mae_new / mae_old - 1))
readr::write_csv(fixed_cmp, file.path(DIR_OUT, "bottom_up_test_fixed_models.csv"))
cat("\nFixed preliminary models, mean of 8 states' MAE (%), ", min(test_years), "-", max(test_years), ":\n", sep = "")
print(as.data.frame(fixed_cmp |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

# 2 and 3. Pseudo real-time, both versions of each method ---------------------------
actual <- inp$a_state |> rename(actual = gsp)
msg("Pseudo real-time comparison,", length(test_years), "vintages in parallel")
res <- parallel::mclapply(test_years, function(t) {
  vin <- vintage_inputs(inp, t - 1)
  m_old <- select_models(oos, t - 1, a_diag$SEL_WINDOW, specs = SPECS_NO_BU) |> filter(rank == 1)
  m_all <- select_models(oos, t - 1, a_diag$SEL_WINDOW, specs = names(SPECS)) |> filter(rank == 1)
  bind_rows(
    method_a(vin, m_old)$est |> mutate(variant = "A: without bottom-up"),
    method_a(vin, m_all)$est |> mutate(variant = "A: with bottom-up"),
    method_b(vin, "prodgap")$est |> mutate(variant = "B: prodgap"),
    method_b(vin, "bugap")$est |> mutate(variant = "B: bugap")
  ) |> filter(fy_of(date) == t) |> group_by(variant, basis, state) |>
    summarise(implied = sum(estimate), .groups = "drop") |>
    mutate(fy = t,
           a_models_bu = sum(spec_of(m_all$model) %in% setdiff(names(SPECS), SPECS_NO_BU)))
}, mc.cores = max(1, min(5, parallel::detectCores() - 1)))
bad <- vapply(res, inherits, logical(1), what = "try-error")
if (any(bad)) stop(res[[which(bad)[1]]])
rt <- bind_rows(res) |> left_join(actual, by = c("state", "basis", "fy")) |> mutate(err_pct = 100 * (implied / actual - 1))
rt_sum <- rt |> group_by(basis, variant, state) |> summarise(mae = mean(abs(err_pct)), .groups = "drop")
rt_tab <- bind_rows(rt_sum, rt_sum |> group_by(basis, variant) |> summarise(state = "Mean of 8", mae = mean(mae), .groups = "drop")) |>
  pivot_wider(names_from = variant, values_from = mae)
readr::write_csv(rt, file.path(DIR_OUT, "bottom_up_test_realtime_errors.csv"))
readr::write_csv(rt_tab, file.path(DIR_OUT, "bottom_up_test_realtime_summary.csv"))
cat("\nPseudo real-time MAE of implied annual GSP (%):\n")
print(as.data.frame(rt_tab |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
cat("\nHow many of the 16 state/basis models Method A chose that use the bottom-up indicator, by vintage:\n")
print(distinct(rt, fy, a_models_bu))

# Queensland Treasury comparison of quarterly growth, full-sample estimates
qld <- read_proc("qld_treasury")
full <- bind_rows(
  method_b(inp, "prodgap")$est |> mutate(method = "B: prodgap"),
  method_b(inp, "bugap")$est |> mutate(method = "B: bugap"),
  method_a(inp, select_models(oos, inp$last_fy, a_diag$SEL_WINDOW, specs = SPECS_NO_BU) |> filter(rank == 1))$est |> mutate(method = "A: without bottom-up"),
  method_a(inp, select_models(oos, inp$last_fy, a_diag$SEL_WINDOW, specs = names(SPECS)) |> filter(rank == 1))$est |> mutate(method = "A: with bottom-up")
) |> select(state, date, basis, method, estimate)
cq <- compare_qld(full, qld$data)$stats
readr::write_csv(cq, file.path(DIR_OUT, "bottom_up_test_qld.csv"))
cat("\nQueensland quarterly growth vs Queensland Treasury:\n")
print(as.data.frame(cq |> select(basis, method, correlation, mean_abs_diff, qt_turning_points, matched_within_1q) |>
                      mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
save_proc(list(fixed = fixed_cmp, rt = rt, rt_tab = rt_tab, qld = cq), "bottom_up_test")
msg("04b_bottom_up_test.R done")
