# 01_data.R -------------------------------------------------------------------
# Download ABS data, build tidy inputs, write the data dictionary, and run the
# seasonal-adjustment and constraint-consistency checks.
#
# Set GSP_REFRESH=false to reuse the spreadsheets already in data/raw.

source("R/utils.R")
source("R/data_helpers.R")
source("R/constraints.R")
refresh <- !identical(tolower(Sys.getenv("GSP_REFRESH", "true")), "false")
dir.create(DIR_OUT, showWarnings = FALSE)

# 1. Downloads -----------------------------------------------------------------
msg("Downloading ABS tables (refresh =", refresh, ")")
na <- get_abs("5206.0", c(1, 2, 3, 5, 6, 26:33, 44, 45), refresh) # quarterly national accounts
sa <- get_abs("5220.0", 1:9, refresh)                              # annual state accounts
tr <- get_abs("5368.0", c("15a", "15b"), refresh)                  # merchandise trade by state (monthly)
lf <- get_abs("6202.0", "17", refresh)                             # monthly hours worked by state

st_tab <- setNames(sprintf("52060%02d_", 26:33), STATES)  # 5206 state detail tables
sa_tab <- setNames(sprintf("52200%02d_", 2:9), STATES)    # 5220 state tables

# 2. National quarterly aggregates ---------------------------------------------
nat_spec <- tibble::tribble(
  ~var,  ~basis,    ~table,     ~series,
  "gdp", "cvm",     "5206001_", "^Gross domestic product: Chain volume measures ;$",
  "gdp", "current", "5206001_", "^Gross domestic product: Current prices ;$",
  "dfd", "cvm",     "5206002_", "^Domestic final demand ;$",
  "dfd", "current", "5206003_", "^Domestic final demand ;$"
)
q_nat <- purrr::pmap_dfr(nat_spec, function(var, basis, table, series) {
  pick(na, table, series, concept = var, state = "AUS", basis = basis) |>
    mutate(var = var, basis = basis)
}) |> pivot_wider(names_from = var, values_from = value)

# Original (unadjusted) GDP, used only for the SA-sums check
q_nat_orig <- purrr::pmap_dfr(nat_spec[1:2, ], function(var, basis, table, series) {
  pick(na, table, series, type = "Original", concept = "gdp (original, check only)",
       state = "AUS", basis = basis) |> mutate(basis = basis)
}) |> rename(gdp_orig = value)

# Trade deflators (goods-specific IPDs are not in 5206; goods & services used)
ipd <- bind_rows(
  pick(na, "5206005_", "^Exports of goods and services ;$", concept = "export IPD", state = "AUS") |> mutate(var = "px"),
  pick(na, "5206005_", "^Imports of goods and services ;$", concept = "import IPD", state = "AUS") |> mutate(var = "pm")
) |> pivot_wider(names_from = var, values_from = value)

# 3. State quarterly series ----------------------------------------------------
msg("Building state quarterly series")
q_state <- purrr::map_dfr(STATES, function(s) {
  nm <- STATE_NAMES[[s]]
  sfd_cvm <- pick(na, st_tab[[s]], "^STATE FINAL DEMAND: Chain volume measures ;$", concept = "sfd", state = s, basis = "cvm")
  sfd_cur <- pick(na, st_tab[[s]], "^STATE FINAL DEMAND: Current prices ;$", concept = "sfd", state = s, basis = "current")
  coe <- pick(na, "5206044_", paste0("^", nm, " ;  Total  ;  Compensation of employees ;$"),
              concept = "coe", state = s, basis = "current")
  # Merchandise trade: monthly, original, current prices -> STL SA -> quarters
  x <- pick(tr, "5368015a", paste0("^", nm, " ;$"), type = "Original", concept = "goods exports (SA by STL here)", state = s, basis = "current", freq = "M")
  m <- pick(tr, "5368015b", paste0("^", nm, " ;$"), type = "Original", concept = "goods imports (SA by STL here)", state = s, basis = "current", freq = "M")
  xq <- m_to_q(stl_sa(x$date, x$value)); mq <- m_to_q(stl_sa(m$date, m$value))
  # Hours: SA for the six states; ABS publishes trend only for NT and ACT
  h_type <- if (s %in% c("NT", "ACT")) "Trend" else "Seasonally Adjusted"
  h <- pick(lf, "62020017", paste0("^Monthly hours worked in all jobs ;  > ", nm, " ;$"), type = h_type,
            concept = "hours worked", state = s, freq = "M")
  hq <- m_to_q(h)
  sfd_cvm |> rename(sfd_cvm = value) |>
    full_join(rename(sfd_cur, sfd_current = value), by = "date") |>
    full_join(rename(coe, coe_current = value), by = "date") |>
    full_join(rename(xq, xg_current = value), by = "date") |>
    full_join(rename(mq, mg_current = value), by = "date") |>
    full_join(rename(hq, hours = value), by = "date") |>
    mutate(state = s)
}) |>
  left_join(ipd, by = "date") |>
  mutate(
    xg_cvm = xg_current / px * 100, mg_cvm = mg_current / pm * 100,
    # volume proxy for labour income: COE deflated by the state's SFD deflator
    coe_cvm = coe_current / (sfd_current / sfd_cvm)
  ) |>
  select(-px, -pm) |>
  arrange(state, date)

# 4. Annual state accounts -----------------------------------------------------
msg("Building annual state accounts")
a_state <- purrr::map_dfr(STATES, function(s) {
  nm <- STATE_NAMES[[s]]
  purrr::map_dfr(BASES, function(b) {
    lab <- if (b == "cvm") "Chain volume measures" else "Current prices"
    gsp <- pick(sa, "5220001_", paste0("^", nm, " ;  Gross state product: ", lab, " ;$"),
                type = "Original", concept = "gsp", state = s, basis = b, freq = "A")
    sfd <- pick(sa, sa_tab[[s]], paste0("STATE FINAL DEMAND: ", lab, " ;$"),
                type = "Original", concept = "sfd (annual, check only)", state = s, basis = b, freq = "A")
    gsp |> rename(gsp = value) |> left_join(rename(sfd, sfd_annual = value), by = "date") |>
      mutate(state = s, basis = b)
  })
}) |> mutate(fy = fy_of(date)) |> select(state, basis, fy, gsp, sfd_annual)

a_nat <- purrr::map_dfr(BASES, function(b) {
  lab <- if (b == "cvm") "Chain volume measures" else "Current prices"
  pick(sa, "5220001_", paste0("^Australia ;  Gross domestic product: ", lab, " ;$"),
       type = "Original", concept = "gdp (annual)", state = "AUS", basis = b, freq = "A") |>
    mutate(basis = b)
}) |> mutate(fy = fy_of(date)) |> select(basis, fy, gdp_annual = value)

# 5. Industry value added ------------------------------------------------------
# 19 ANZSIC divisions + ownership of dwellings. Annual by state (5220 tables
# 2-9); national quarterly SA (5206 table 6 CVM; table 45 current prices, which
# only starts in September quarter 2002).
ind_names <- sa |>
  filter(grepl(sa_tab[["NSW"]], table_no), grepl(";  Industry gross value added: Current prices ;$", series),
         !grepl("^Total all industries", series)) |>
  distinct(series) |> pull(series) |> sub(pattern = " ;  Industry gross value added.*$", replacement = "")
stopifnot(length(ind_names) == 20)
rx <- function(x) gsub("([()])", "\\\\\\1", x)

a_ind <- purrr::map_dfr(STATES, function(s) purrr::map_dfr(BASES, function(b) purrr::map_dfr(ind_names, function(i) {
  lab <- if (b == "cvm") "Chain volume measures" else "Current prices"
  pick(sa, sa_tab[[s]], paste0("^", rx(i), " ;  Industry gross value added: ", lab, " ;$"),
       type = "Original", concept = paste("industry GVA:", i), state = s, basis = b, freq = "A") |>
    mutate(state = s, basis = b, industry = i)
}))) |> mutate(fy = fy_of(date)) |> select(state, basis, industry, fy, gva = value)

q_ind <- purrr::map_dfr(ind_names, function(i) {
  bind_rows(
    pick(na, "5206006_", paste0("^", rx(i), " ;$"), concept = paste("national industry GVA:", i), state = "AUS", basis = "cvm") |>
      mutate(basis = "cvm"),
    pick(na, "5206045_", paste0("^", rx(i), " ;  Gross value added at basic prices ;$"),
         concept = paste("national industry GVA:", i), state = "AUS", basis = "current") |> mutate(basis = "current")
  ) |> mutate(industry = i)
}) |> select(basis, industry, date, gva = value)

# 6. Common sample -------------------------------------------------------------
# First benchmark year is 1989-90, so the quarterly sample starts Sep 1989 and
# ends at the latest quarter for which national GDP and all SFD are published.
q_start <- as.Date("1989-09-01")
q_end <- min(max(q_nat$date), max(q_state$date[!is.na(q_state$sfd_cvm)]))
q_nat <- filter(q_nat, date >= q_start, date <= q_end) |> mutate(fy = fy_of(date))
q_state <- filter(q_state, date >= q_start, date <= q_end) |> mutate(fy = fy_of(date))
q_ind <- filter(q_ind, date >= as.Date("1989-06-01"), date <= q_end)
stopifnot(nrow(q_nat) == 2 * length(unique(q_nat$date)))
core <- c("sfd_cvm", "sfd_current", "coe_current", "xg_current", "mg_current", "hours")
miss <- colSums(is.na(q_state[, core]))
if (any(miss > 0)) {
  stop_data(paste("missing quarterly state observations in:", paste(names(miss)[miss > 0], collapse = ", ")),
            "the corresponding ABS table (see output/data_dictionary.csv)")
}
last_bench_fy <- max(a_state$fy)
msg("Quarterly sample:", format(q_start), "to", format(q_end),
    "| annual benchmarks: FY1990 to FY", last_bench_fy,
    "| unbenchmarked quarters:", sum(unique(q_nat$date) > as.Date(sprintf("%d-06-01", last_bench_fy))))

# 7. CHECK: do seasonally adjusted quarters sum to the annual figures? ----------
msg("Checking seasonally adjusted quarters against annual totals")
full_fy <- q_nat |> count(basis, fy) |> filter(n == 4) |> distinct(fy) |> pull(fy)
chk_nat <- q_nat |> filter(fy %in% full_fy) |>
  left_join(q_nat_orig, by = c("date", "basis")) |>
  group_by(basis, fy) |> summarise(sa_sum = sum(gdp), orig_sum = sum(gdp_orig), .groups = "drop") |>
  inner_join(a_nat, by = c("basis", "fy")) |>
  transmute(series = "GDP", state = "AUS", basis, fy, sa_sum, orig_sum, annual = gdp_annual,
            sa_vs_annual_pct = 100 * (sa_sum / annual - 1), orig_vs_annual_pct = 100 * (orig_sum / annual - 1))
chk_sfd <- q_state |> filter(fy %in% full_fy) |>
  pivot_longer(c(sfd_cvm, sfd_current), names_to = "basis", names_prefix = "sfd_", values_to = "sfd") |>
  group_by(state, basis, fy) |> summarise(sa_sum = sum(sfd), .groups = "drop") |>
  inner_join(a_state, by = c("state", "basis", "fy")) |>
  transmute(series = "SFD", state, basis, fy, sa_sum, orig_sum = NA_real_, annual = sfd_annual,
            sa_vs_annual_pct = 100 * (sa_sum / annual - 1), orig_vs_annual_pct = NA_real_)
chk_sa <- bind_rows(chk_nat, chk_sfd)
readr::write_csv(chk_sa, file.path(DIR_OUT, "check_sa_sums.csv"))
chk_sa_summary <- chk_sa |> group_by(series, basis) |>
  summarise(max_abs_pct = max(abs(sa_vs_annual_pct)), mean_abs_pct = mean(abs(sa_vs_annual_pct)),
            worst_fy = fy[which.max(abs(sa_vs_annual_pct))], worst_state = state[which.max(abs(sa_vs_annual_pct))],
            n_over_tol = sum(abs(sa_vs_annual_pct) > 100 * TOL), n = n(), .groups = "drop")
print(as.data.frame(chk_sa_summary), digits = 3)

# 8. Constraint targets and their mutual consistency ---------------------------
# See build_target() in R/constraints.R. k = (sum of state annual GSP) / (sum of
# the four national SA quarters) is the adjustment needed for the two constraint
# sets to agree. It decomposes into
#   k_nonadd = sum of state annual GSP / national annual GDP (CVM non-additivity;
#              ~1 in current prices), and
#   k_sa     = national annual GDP / sum of four SA quarters (SA-sum or vintage gap).
target_method <- Sys.getenv("GSP_TARGET", "denton") # "denton" (default) or "prorata"
targets_fy <- a_state |> group_by(basis, fy) |> summarise(states_sum = sum(gsp), .groups = "drop") |>
  inner_join(a_nat, by = c("basis", "fy")) |>
  inner_join(q_nat |> group_by(basis, fy) |> summarise(q_sum = sum(gdp), .groups = "drop"), by = c("basis", "fy")) |>
  group_by(basis) |> arrange(fy, .by_group = TRUE) |>
  mutate(k_nonadd = states_sum / gdp_annual, k_sa = gdp_annual / q_sum, k = states_sum / q_sum,
         raw_inconsistency_pct = 100 * (k - 1), k_step_pct = 100 * (k / lag(k) - 1)) |> ungroup()
readr::write_csv(targets_fy, file.path(DIR_OUT, "check_constraint_consistency.csv"))
cat("\nConsistency of the two constraint sets BEFORE adjustment\n(sum of state annual GSP vs sum of the 4 national SA quarters, per cent):\n")
print(as.data.frame(targets_fy |> group_by(basis) |>
  summarise(max_abs_raw = max(abs(raw_inconsistency_pct)), mean_abs_raw = mean(abs(raw_inconsistency_pct)),
            max_abs_nonadd = 100 * max(abs(k_nonadd - 1)), max_abs_sa = 100 * max(abs(k_sa - 1)),
            max_abs_k_step = max(abs(k_step_pct), na.rm = TRUE), mean_abs_k_step = mean(abs(k_step_pct), na.rm = TRUE))),
  digits = 3)
q_target <- build_target(q_nat, a_state, a_nat, method = target_method)
cat("\nNational target method:", target_method, "| max abs difference from the per-year pro-rata target (per cent):\n")
print(as.data.frame(q_target |> group_by(basis) |> summarise(max_abs_diff_pct = 100 * max(abs(k / k_prorata - 1)))), digits = 3)

# 9. Save ----------------------------------------------------------------------
save_proc(q_nat, "q_nat"); save_proc(q_state, "q_state"); save_proc(a_state, "a_state")
save_proc(a_nat, "a_nat"); save_proc(a_ind, "a_ind"); save_proc(q_ind, "q_ind")
save_proc(q_target, "q_target"); save_proc(targets_fy, "targets_fy")
save_proc(target_method, "target_method")
save_proc(chk_sa, "chk_sa"); save_proc(list(q_start = q_start, q_end = q_end, last_bench_fy = last_bench_fy,
                                            run_date = Sys.Date()), "meta")
dd <- dict_get()
readr::write_csv(dd, file.path(DIR_OUT, "data_dictionary.csv"))
save_proc(dd, "data_dictionary")
msg(nrow(dd), "ABS series recorded in output/data_dictionary.csv")

# 10. Queensland Treasury quarterly state accounts (cross-check only) -----------
# Done last so that a failure here leaves the ABS inputs saved. 06_validate.R
# needs this file and will stop with the same message if it is absent.
qld_file <- get_qld_workbook()
qld <- read_qld(qld_file)
save_proc(list(file = qld_file, data = qld), "qld_treasury")
msg("01_data.R done. Queensland workbook:", basename(qld_file), "| quarters:", format(min(qld$date)), "to", format(max(qld$date)))
