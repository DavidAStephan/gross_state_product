# Constraint targets and checks -------------------------------------------------

# National quarterly target for the cross-sectional constraint.
#
# The quarterly national SA GDP series cannot be used as-is because
#   (a) CVMs are not additive: the states' annual CVMs do not sum to national
#       annual CVM, and
#   (b) the SA quarters do not always sum to the annual figure (current prices;
#       and the latest year when the quarterly accounts are a newer vintage than
#       the annual state accounts).
# Both are removed by benchmarking national GDP to the annual sum of states:
#   method = "prorata": target = GDP_q x k_fy, k_fy = sum of states / sum of the
#       four national quarters (the literal definition in the brief; k steps at
#       each July).
#   method = "denton":  proportional Denton (Cholette) benchmark of GDP_q to the
#       same annual totals. Annual sums are identical to "prorata", but the
#       ratio moves smoothly between years, so no step is pushed into
#       September-quarter growth.
# Quarters after the last benchmark year carry the last ratio forward.
#
# q_nat: date, basis, gdp, fy. a_state: state, basis, fy, gsp (benchmark years only).
build_target <- function(q_nat, a_state, a_nat = NULL, method = c("denton", "prorata")) {
  method <- match.arg(method)
  purrr::map_dfr(unique(q_nat$basis), function(b) {
    q <- q_nat[q_nat$basis == b, ] |> arrange(date)
    stopifnot(fyq_of(q$date[1]) == 1)
    a <- a_state[a_state$basis == b, ] |> group_by(fy) |>
      summarise(states_sum = sum(gsp), n = n(), .groups = "drop") |> arrange(fy)
    stopifnot(all(a$n == length(STATES)), a$fy[1] == q$fy[1], all(diff(a$fy) == 1))
    qs <- q |> group_by(fy) |> summarise(q_sum = sum(gdp), nq = n(), .groups = "drop")
    a <- left_join(a, qs, by = "fy")
    stopifnot(all(a$nq == 4))
    a$k <- a$states_sum / a$q_sum
    q$benchmarked <- q$fy %in% a$fy
    if (method == "prorata") {
      q$k <- a$k[match(q$fy, a$fy)]
      q <- tidyr::fill(q, k, .direction = "down")
    } else {
      hf <- stats::ts(q$gdp, start = c(q$fy[1], 1), frequency = 4)
      lf <- stats::ts(a$states_sum, start = a$fy[1], frequency = 1)
      fit <- tempdisagg::td(lf ~ 0 + hf, method = "denton-cholette", conversion = "sum")
      q$k <- as.numeric(stats::predict(fit)) / q$gdp
    }
    q$target <- q$gdp * q$k
    # the two constraint sets must now be mutually consistent
    chk <- q[q$benchmarked, ] |> group_by(fy) |> summarise(t = sum(target), .groups = "drop") |>
      left_join(a, by = "fy")
    stopifnot(max(rel_err(chk$t, chk$states_sum)) < 1e-9)
    q$k_prorata <- tidyr::fill(tibble::tibble(k = a$k[match(q$fy, a$fy)]), k, .direction = "down")$k
    if (!is.null(a_nat)) {
      an <- a_nat[a_nat$basis == b, ]
      q$k_nonadd <- (a$states_sum / an$gdp_annual[match(a$fy, an$fy)])[match(q$fy, a$fy)]
    }
    q
  })
}

source("R/indicators.R")

# Inputs --------------------------------------------------------------------------
# All estimators take the same `inp` list: quarterly panel, quarterly national
# GDP, annual GSP benchmarks and the last benchmark year.
load_inputs <- function() {
  a_state <- read_proc("a_state")
  a_ind <- read_proc("a_ind"); q_ind <- read_proc("q_ind")
  last_fy <- max(a_state$fy)
  list(panel = add_bottom_up(read_proc("panel"), a_ind, q_ind, last_fy), q_nat = read_proc("q_nat"),
       a_state = select(a_state, state, basis, fy, gsp), last_fy = last_fy,
       a_ind = a_ind, q_ind = q_ind, target_method = read_proc("target_method"))
}

# Pseudo real-time vintage: annual GSP up to and including `last_fy`, quarterly
# data up to June of `last_fy + 1` (so exactly one unbenchmarked year).
# Caveat: this truncates today's data vintage; it does not undo ABS revisions.
vintage_inputs <- function(inp, last_fy) {
  qmax <- as.Date(sprintf("%d-06-01", last_fy + 1))
  inp$a_state <- filter(inp$a_state, fy <= last_fy)
  inp$panel <- filter(inp$panel, date <= qmax)
  inp$q_nat <- filter(inp$q_nat, date <= qmax)
  inp$last_fy <- last_fy
  # the bottom-up indicator uses annual industry data, so rebuild it from this vintage's years only
  if (!is.null(inp$a_ind)) {
    inp$a_ind <- filter(inp$a_ind, fy <= last_fy)
    inp$panel <- add_bottom_up(inp$panel, inp$a_ind, inp$q_ind, last_fy)
  }
  inp
}

# Matrices for one price basis (rows = quarters or years, columns = STATES)
prep_basis <- function(inp, b) {
  panel <- inp$panel[inp$panel$basis == b, ]
  dates <- sort(unique(panel$date))
  a <- inp$a_state[inp$a_state$basis == b, ]
  years <- sort(unique(a$fy))
  tgt <- build_target(inp$q_nat[inp$q_nat$basis == b, ], a, method = inp$target_method)
  stopifnot(identical(tgt$date, dates))
  fy <- fy_of(dates)
  SFD <- to_wide(panel, "sfd", dates)
  GSPa <- as.matrix(a |> select(fy, state, gsp) |> pivot_wider(names_from = state, values_from = gsp) |>
                      arrange(fy) |> select(all_of(STATES)))
  rownames(GSPa) <- years
  C <- outer(years, fy, "==") * 1 # annual aggregation matrix (Y x T)
  rownames(C) <- years
  list(basis = b, dates = dates, fy = fy, years = years, bench = fy %in% years,
       panel = panel, SFD = SFD, GSPa = GSPa, C = C,
       gap_a = GSPa - C %*% SFD,          # annual gap to be distributed
       target = tgt$target,               # national quarterly target for sum of states
       natgap = tgt$target - rowSums(SFD), # national quarterly gap
       q_target = tgt)
}

mat_to_long <- function(M, dates, b, name = "estimate") {
  d <- tibble::tibble(date = rep(dates, times = ncol(M)), state = rep(colnames(M), each = nrow(M)),
                      basis = b, v = as.numeric(M))
  names(d)[4] <- name
  d
}

# Check both accounting constraints on a long table of estimates
# est: state, date, basis, estimate. Returns one row per basis.
check_constraints <- function(est, a_state, q_target) {
  est <- mutate(est, fy = fy_of(date))
  temporal <- est |> group_by(state, basis, fy) |>
    summarise(q_sum = sum(estimate), nq = n(), .groups = "drop") |>
    inner_join(a_state, by = c("state", "basis", "fy")) |>
    mutate(rel = rel_err(q_sum, gsp))
  cross <- est |> group_by(basis, date) |> summarise(s_sum = sum(estimate), ns = n(), .groups = "drop") |>
    inner_join(select(q_target, basis, date, target), by = c("basis", "date")) |>
    mutate(rel = rel_err(s_sum, target))
  full_join(
    temporal |> group_by(basis) |> summarise(temporal_max_rel = max(rel), n_temporal = n(), all_4q = all(nq == 4)),
    cross |> group_by(basis) |> summarise(cross_max_rel = max(rel), n_cross = n(), all_states = all(ns == length(STATES))),
    by = "basis"
  )
}
