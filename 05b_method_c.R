# 05b_method_c.R ---------------------------------------------------------------
# Method C (alternate): mixed-frequency VAR with stochastic volatility and a
# Dirichlet-Laplace prior, following Koop, McIntyre, Mitchell and Poon (2020).
# See R/method_c.R for the model and the few places it departs from the paper.
# Each posterior draw is reconciled to the two accounting constraints.
# Slow: about 15 minutes (the two price bases run in parallel).

source("R/utils.R")
source("R/constraints.R")
source("R/method_a.R") # reconcile_dfm()
source("R/method_c.R")

inp <- load_inputs()
msg("MF-VAR-SV: Gibbs sampler,", MFVAR_DEFAULT$n_iter, "iterations per price basis")
res <- method_c(inp)
est <- res$est |> mutate(method = "method_c")

cat("\nSampler diagnostics (growth in per cent per quarter):\n")
print(as.data.frame(res$diag |> mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
readr::write_csv(res$diag, file.path(DIR_OUT, "method_c_diagnostics.csv"))

q_target <- read_proc("q_target")
chk <- check_constraints(est, inp$a_state, q_target)
cat("\nConstraint check (max relative error):\n"); print(as.data.frame(chk), digits = 3)
stopifnot(chk$temporal_max_rel < TOL, chk$cross_max_rel < TOL)

readr::write_csv(res$vol, file.path(DIR_OUT, "method_c_volatility.csv"))
save_proc(est, "est_method_c"); save_proc(res$diag, "method_c_diag"); save_proc(res$vol, "method_c_vol")
msg("05b_method_c.R done")
