# 05_method_b.R ----------------------------------------------------------------
# Method B: joint state-space model (KFAS). Smoothed quarterly GSP with 68% and
# 90% intervals; the unbenchmarked quarters are the model's nowcasts.

source("R/utils.R")
source("R/constraints.R")
source("R/method_a.R") # for state_X()
source("R/method_b.R")

inp <- load_inputs()
msg("Estimating the state-space model (two price bases)")
res <- method_b(inp)
est <- res$est |> mutate(method = "method_b")

cat("\nEstimated parameters (units: per cent of the state's quarterly SFD):\n")
print(as.data.frame(res$pars |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\nLog-likelihood with and without the common factor:\n"); print(as.data.frame(res$fit))
if (any(res$pars$stage1_convergence != 0) || any(res$fit$stage2_convergence != 0))
  warning("An optimiser did not report convergence - inspect output/method_b_parameters.csv")
readr::write_csv(res$pars, file.path(DIR_OUT, "method_b_parameters.csv"))

# Constraint check: must hold exactly because both constraints have zero
# measurement error (including the June cross-sectional rows left out as redundant)
q_target <- read_proc("q_target")
chk <- check_constraints(est, inp$a_state, q_target)
cat("\nConstraint check (max relative error):\n"); print(as.data.frame(chk), digits = 3)
stopifnot(chk$temporal_max_rel < TOL, chk$cross_max_rel < TOL)

# Width of the intervals, and the nowcast quarters
unc <- est |> mutate(se_pct = 100 * se / estimate) |> group_by(basis, state, benchmarked) |>
  summarise(mean_se_pct = mean(se_pct), max_se_pct = max(se_pct), .groups = "drop")
cat("\nStandard error of quarterly GSP level, % of level (benchmarked vs nowcast quarters):\n")
print(as.data.frame(unc |> pivot_wider(names_from = benchmarked, values_from = c(mean_se_pct, max_se_pct)) |>
                      mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
readr::write_csv(unc, file.path(DIR_OUT, "method_b_uncertainty.csv"))

save_proc(est, "est_method_b"); save_proc(list(pars = res$pars, fit = res$fit, unc = unc), "method_b_diag")
msg("05_method_b.R done")
