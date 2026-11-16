# 04_method_a.R ----------------------------------------------------------------
# Method A: temporal disaggregation of the annual gap on quarterly indicators
# (model chosen per state by out-of-sample performance), then Di Fonzo-Marini
# two-way reconciliation of GSP levels.

source("R/utils.R")
source("R/constraints.R")
source("R/method_a.R")

inp <- load_inputs()
N_TEST <- 10      # pseudo real-time test years (used again in 06_validate.R)
SEL_WINDOW <- 10  # years of out-of-sample errors used to choose a model

# 1. Out-of-sample errors of every candidate preliminary model -------------------
# Target years: enough history that each of the last N_TEST test years can pick
# its model from the SEL_WINDOW years before it.
target_years <- (inp$last_fy - N_TEST - SEL_WINDOW + 1):inp$last_fy
msg("Out-of-sample errors for", nrow(model_grid()), "candidate models,", length(target_years), "target years")
oos <- oos_prelim_errors(inp, target_years)
save_proc(oos, "oos_prelim")

# 2. Model selection for the production run ---------------------------------------
sel_all <- select_models(oos, inp$last_fy, SEL_WINDOW)
chosen <- sel_all |> filter(rank == 1) |> select(basis, state, model, rmse, mae)
# context: how the chosen model compares with plain Fernandez on prodgap, and the spread
ctx <- sel_all |> group_by(basis, state) |>
  summarise(rmse_fernandez_s1 = rmse[model == "fernandez:s1_prodgap"],
            rmse_median_model = median(rmse), rmse_worst_model = max(rmse), .groups = "drop")
chosen <- left_join(chosen, ctx, by = c("basis", "state"))
cat("\nChosen preliminary model per state (out-of-sample RMSE of implied annual GSP, % of GSP):\n")
print(as.data.frame(chosen |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
readr::write_csv(sel_all, file.path(DIR_OUT, "method_a_model_selection_all.csv"))
readr::write_csv(chosen, file.path(DIR_OUT, "method_a_model_chosen.csv"))

# 3. Estimate and reconcile ---------------------------------------------------------
res <- method_a(inp, chosen, keep_fits = TRUE)
est <- res$est |> mutate(method = "method_a")
cat("\nReconciliation diagnostics:\n"); print(as.data.frame(res$recon), digits = 3)

# coefficients of the chosen models
coefs <- purrr::map_dfr(BASES, function(b) purrr::map_dfr(STATES, function(s) {
  f <- res$fits[[b]][[s]]
  if (inherits(f, "td")) {
    sm <- summary(f)$coefficients
    tibble::tibble(basis = b, state = s, term = rownames(sm), estimate = sm[, 1], t_value = sm[, 3],
                   rho = if (is.null(f$rho)) NA_real_ else f$rho)
  } else {
    tibble::tibble(basis = b, state = s, term = names(f$beta), estimate = f$beta, t_value = NA_real_, rho = NA_real_)
  }
}))
readr::write_csv(coefs, file.path(DIR_OUT, "method_a_coefficients.csv"))

# size of the reconciliation adjustment
adj <- est |> mutate(adj_pct = 100 * (estimate / prelim - 1)) |> group_by(basis, state) |>
  summarise(mean_abs_adj_pct = mean(abs(adj_pct)), max_abs_adj_pct = max(abs(adj_pct)), .groups = "drop") |>
  left_join(res$weights, by = c("basis", "state"))
cat("\nReconciliation adjustment to preliminary levels (%), and the residual s.d. that sets the weights:\n")
print(as.data.frame(adj |> mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
readr::write_csv(adj, file.path(DIR_OUT, "method_a_reconciliation_adjustment.csv"))

# 4. Constraint check ------------------------------------------------------------------
q_target <- read_proc("q_target")
chk <- check_constraints(est, inp$a_state, q_target)
cat("\nConstraint check (max relative error):\n"); print(as.data.frame(chk), digits = 3)
stopifnot(chk$temporal_max_rel < TOL, chk$cross_max_rel < TOL)

save_proc(est, "est_method_a"); save_proc(chosen, "method_a_chosen")
save_proc(list(weights = res$weights, recon = res$recon, coefs = coefs, adj = adj,
               N_TEST = N_TEST, SEL_WINDOW = SEL_WINDOW), "method_a_diag")
msg("04_method_a.R done")
