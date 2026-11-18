# Automated checks on the published output and on the reconciliation routine.
# Run from 06_validate.R (which sets GSP_ROOT), or with
#   GSP_ROOT=$(pwd) Rscript -e 'testthat::test_dir("tests/testthat")'

root <- Sys.getenv("GSP_ROOT", unset = normalizePath(file.path("..", "..")))
old <- setwd(root)
source("R/utils.R"); source("R/constraints.R"); source("R/method_a.R")
out <- readr::read_csv(file.path(DIR_OUT, "quarterly_gsp.csv"), show_col_types = FALSE)
inp <- load_inputs()
q_target <- read_proc("q_target")
setwd(old)

est <- out |> transmute(state, date = quarter, basis = price_basis, method, estimate)
models <- intersect(c("method_a", "method_b", "method_c"), unique(est$method))

test_that("the two constraint sets are mutually consistent", {
  chk <- q_target |> filter(benchmarked) |> group_by(basis, fy) |> summarise(t = sum(target), .groups = "drop") |>
    inner_join(inp$a_state |> group_by(basis, fy) |> summarise(a = sum(gsp), .groups = "drop"), by = c("basis", "fy"))
  expect_equal(nrow(chk), 2 * length(unique(inp$a_state$fy)))
  expect_lt(max(abs(chk$t - chk$a) / chk$a), 1e-9)
})

for (m in models) for (b in BASES) {
  e <- filter(est, method == m, basis == b)
  chk <- check_constraints(e, inp$a_state, q_target)

  test_that(paste(m, b, "temporal: four quarters sum to annual GSP in every benchmark year"), {
    expect_true(chk$all_4q)
    expect_equal(chk$n_temporal, length(STATES) * length(unique(inp$a_state$fy)))
    expect_lt(chk$temporal_max_rel, TOL)
  })

  test_that(paste(m, b, "cross-sectional: states sum to the national target in every quarter"), {
    expect_true(chk$all_states)
    expect_equal(chk$n_cross, length(unique(q_target$date)))
    expect_lt(chk$cross_max_rel, TOL)
  })

  test_that(paste(m, b, "unbenchmarked quarters satisfy the cross-sectional constraint"), {
    ub <- q_target$date[q_target$basis == b & !q_target$benchmarked]
    # Right after an annual state accounts release and before the next national
    # accounts, every quarter is benchmarked and there is nothing to test here.
    if (length(ub) == 0) skip("no unbenchmarked quarters in this data release")
    s <- e |> filter(date %in% ub) |> group_by(date) |> summarise(s = sum(estimate)) |>
      inner_join(filter(q_target, basis == b), by = "date")
    expect_lt(max(abs(s$s - s$target) / s$target), TOL)
  })

  test_that(paste(m, b, "estimates are complete and positive"), {
    expect_false(anyNA(e$estimate))
    expect_true(all(e$estimate > 0))
    expect_equal(nrow(e), length(STATES) * length(unique(q_target$date)))
  })
}

test_that("Method B intervals are ordered and nested", {
  b <- filter(out, method == "method_b")
  expect_false(anyNA(b[, c("lower_68", "upper_68", "lower_90", "upper_90")]))
  expect_true(all(b$lower_90 <= b$lower_68 & b$lower_68 <= b$estimate &
                    b$estimate <= b$upper_68 & b$upper_68 <= b$upper_90))
})

if ("method_c" %in% models) test_that("Method C intervals are ordered and contain the estimate", {
  cc <- filter(out, method == "method_c")
  expect_false(anyNA(cc[, c("lower_68", "upper_68", "lower_90", "upper_90")]))
  # posterior quantiles of reconciled draws; the estimate is their mean
  expect_true(all(cc$lower_90 <= cc$lower_68 & cc$lower_68 <= cc$upper_68 & cc$upper_68 <= cc$upper_90))
  expect_true(all(cc$lower_90 <= cc$estimate & cc$estimate <= cc$upper_90))
})

test_that("reconciliation handles the redundant constraints and rejects inconsistent ones", {
  old <- setwd(root); on.exit(setwd(old))
  pb <- prep_basis(inp, "current")
  P <- baseline_like <- pb$SFD * (pb$GSPa / (pb$C %*% pb$SFD))[match(pmin(pb$fy, max(pb$years)), pb$years), ]
  r <- suppressWarnings(reconcile_dfm(P, pb$GSPa, pb$target, pb$fy, pb$years))
  expect_equal(r$n_dropped, length(pb$years))       # one redundant constraint per benchmark year
  expect_equal(r$rank, r$n_constraints - r$n_dropped)
  expect_lt(r$err_temporal, 1e-8); expect_lt(r$err_cross, 1e-8)
  bad <- pb$target; bad[1] <- bad[1] * 1.01          # break the consistency of the two sets
  expect_error(suppressWarnings(reconcile_dfm(P, pb$GSPa, bad, pb$fy, pb$years)), "inconsistent")
})
