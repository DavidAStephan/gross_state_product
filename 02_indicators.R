# 02_indicators.R --------------------------------------------------------------
# Build the production-side indicator and the quarterly indicator panel, and
# the annual gap (GSP - SFD) that the models explain.

source("R/utils.R")
source("R/indicators.R")

q_state <- read_proc("q_state"); q_nat <- read_proc("q_nat")
a_state <- read_proc("a_state"); a_ind <- read_proc("a_ind"); q_ind <- read_proc("q_ind")

prod <- build_prod_indicator(a_ind, q_ind, q_nat)
panel <- build_panel(q_state, prod)
stopifnot(!anyNA(panel[, c("sfd", "prod", "prodgap", "xg", "mg", "coe", "hours")]),
          nrow(panel) == length(STATES) * 2 * length(unique(q_nat$date)))

# Annual gap = annual GSP - sum of the four quarterly SA SFD figures. Defining
# it against the quarterly SFD actually used (not the annual SFD published in
# 5220) means quarterly SFD + quarterly gap sums exactly to annual GSP even
# though SA SFD does not sum exactly to annual SFD.
a_gap <- panel |> group_by(state, basis, fy) |> summarise(sfd_q_sum = sum(sfd), nq = n(), .groups = "drop") |>
  inner_join(a_state, by = c("state", "basis", "fy")) |>
  mutate(gap = gsp - sfd_q_sum, gap_share = gap / gsp)
stopifnot(all(a_gap$nq == 4))

save_proc(prod, "prod"); save_proc(panel, "panel"); save_proc(a_gap, "a_gap")

# Diagnostics -------------------------------------------------------------------
# (1) How well does the indicator track annual GSP growth? (2) how big is the
# step that chaining removes?
diag <- panel |> group_by(state, basis, fy) |>
  summarise(prod = sum(prod), sfd = sum(sfd), xg = sum(xg), .groups = "drop") |>
  inner_join(select(a_state, state, basis, fy, gsp), by = c("state", "basis", "fy")) |>
  group_by(state, basis) |> arrange(fy, .by_group = TRUE) |>
  mutate(across(c(prod, sfd, gsp), ~ 100 * (.x / lag(.x) - 1), .names = "g_{.col}")) |>
  summarise(cor_gsp_prod = cor(g_gsp, g_prod, use = "complete.obs"),
            cor_gsp_sfd = cor(g_gsp, g_sfd, use = "complete.obs"),
            mae_prod = mean(abs(g_gsp - g_prod), na.rm = TRUE),
            mae_sfd = mean(abs(g_gsp - g_sfd), na.rm = TRUE), .groups = "drop")
step <- prod |> filter(!is.na(prod_step), fyq_of(date) == 1) |> group_by(state, basis) |> arrange(date, .by_group = TRUE) |>
  mutate(step_pct = 100 * ((prod_step / lag(prod_step)) / (prod / lag(prod)) - 1)) |>
  summarise(mean_abs_sep_step_pct = mean(abs(step_pct), na.rm = TRUE), .groups = "drop")
ind_diag <- left_join(diag, step, by = c("state", "basis")) |>
  left_join(a_gap |> group_by(state, basis) |> summarise(mean_gap_share = mean(gap_share), min_gap_share = min(gap_share),
                                                           max_gap_share = max(gap_share), .groups = "drop"), by = c("state", "basis"))
readr::write_csv(ind_diag, file.path(DIR_OUT, "indicator_diagnostics.csv"))
save_proc(ind_diag, "ind_diag")
cat("\nAnnual growth: correlation and mean abs difference (ppt) of GSP with the production indicator and with SFD;\n",
    "mean abs September-quarter step removed by chaining (%); gap as a share of GSP\n", sep = "")
print(as.data.frame(ind_diag |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
msg("02_indicators.R done")
