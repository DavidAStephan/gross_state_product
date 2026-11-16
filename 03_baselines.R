# 03_baselines.R ---------------------------------------------------------------
# Quick baselines for comparison: pro-rata on SFD and Denton on SFD, each
# scaled to the national total.

source("R/utils.R")
source("R/constraints.R")
source("R/baselines.R")

inp <- load_inputs()
est <- run_baselines(inp)
save_proc(est, "est_baselines")

q_target <- read_proc("q_target")
chk <- est |> group_by(method) |> group_modify(~ check_constraints(.x, inp$a_state, q_target)) |> ungroup()
cat("\nConstraint check for baselines (max relative error). The temporal constraint is\n",
    "NOT expected to hold exactly after scaling to the national total:\n", sep = "")
print(as.data.frame(chk), digits = 3)
readr::write_csv(chk, file.path(DIR_OUT, "check_baselines.csv"))
msg("03_baselines.R done")
