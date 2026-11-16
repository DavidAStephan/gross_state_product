# 07_dashboard_export.R ---------------------------------------------------------
# Optional: writes the small JSON files behind the results dashboard
# (output/dashboard/). Each file is {rows, meta}; meta carries this script.
# Run after 06_validate.R.

source("R/utils.R")
dir.create(file.path(DIR_OUT, "dashboard"), showWarnings = FALSE)
v <- read_proc("validation")
out <- readr::read_csv(file.path(DIR_OUT, "quarterly_gsp.csv"), show_col_types = FALSE)
me <- paste(readLines("07_dashboard_export.R"), collapse = "\n")
write_ds <- function(rows, name, input) {
  meta <- list(list(label = "07_dashboard_export.R", value = me, type = "code", lang = "r"),
               list(label = "Input", value = input),
               list(label = "ABS data to", value = format(max(out$quarter), "%B %Y")))
  jsonlite::write_json(list(rows = rows, meta = meta), file.path(DIR_OUT, "dashboard", paste0(name, ".json")),
                       auto_unbox = TRUE, digits = 6, na = "null")
}
models <- setdiff(unique(out$method), c("baseline_prorata", "baseline_denton"))

gsp <- out |> filter(method %in% models) |> group_by(state, price_basis, method) |> arrange(quarter, .by_group = TRUE) |>
  mutate(growth = round(100 * (estimate / lag(estimate) - 1), 3)) |> ungroup() |>
  transmute(key = paste(state, price_basis, method, quarter, sep = "|"), state, quarter = as.character(quarter),
            basis = price_basis, method, estimate = round(estimate, 1), growth,
            lower_90 = round(lower_90, 1), upper_90 = round(upper_90, 1), benchmarked)
write_ds(gsp, "gsp", "output/quarterly_gsp.csv")

rt <- v$rt_sum |> filter(method != "method_a_prelim") |>
  transmute(key = paste(basis, method, state, sep = "|"), basis, method, state,
            mae = round(mae, 3), rmse = round(rmse, 3), bias = round(bias, 3), max_abs = round(max_abs, 3))
write_ds(rt, "realtime", "output/validation_realtime_summary.csv")

qld <- v$qld$stats |> transmute(key = paste(basis, method, sep = "|"), basis, method, to = as.character(to), quarters = n,
                                correlation = round(correlation, 3), mean_abs_diff = round(mean_abs_diff, 3),
                                treasury_turning_points = qt_turning_points, matched_within_1q)
write_ds(qld, "qld", "output/validation_qld_treasury.csv")

# pairwise agreement of quarterly growth between the model-based methods
g <- gsp |> select(state, basis, quarter, method, growth) |> filter(!is.na(growth))
pairs <- utils::combn(sort(models), 2, simplify = FALSE)
agree <- purrr::map_dfr(pairs, function(p) {
  inner_join(filter(g, method == p[1]), filter(g, method == p[2]), by = c("state", "basis", "quarter")) |>
    group_by(basis, state) |>
    summarise(pair = paste(p, collapse = " vs "), correlation = round(cor(growth.x, growth.y), 3),
              mean_abs_diff = round(mean(abs(growth.x - growth.y)), 3),
              quarters_over_1ppt = sum(abs(growth.x - growth.y) > 1), quarters = n(), .groups = "drop")
}) |> mutate(key = paste(basis, state, pair, sep = "|"))
write_ds(agree, "agreement", "output/quarterly_gsp.csv")
msg("07_dashboard_export.R done:", paste(list.files(file.path(DIR_OUT, "dashboard")), collapse = ", "))
