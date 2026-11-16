# Production-side indicator -----------------------------------------------------
#
# For each state: national quarterly industry GVA weighted by the state's share
# of each industry. Shares for all four quarters of financial year y are the
# state's shares in year y-1 (the latest annual state accounts available when
# year y is being estimated), so there is no look-ahead in unbenchmarked
# quarters.
#
# Judgement call: when the shares are updated each September quarter the level
# of a naive weighted sum steps by the amount the state out- or under-performed
# the nation in the previous year. That step would be read as September-quarter
# growth. The indicator is therefore CHAINED: growth into quarter q is
#     sum_i share[i, fy(q)-1] * N[i, q] / sum_i share[i, fy(q)-1] * N[i, q-1]
# (same weights in numerator and denominator), cumulated from the first quarter.
# The unchained version is also returned (prod_step) so the step can be seen.
#
# Current prices: national quarterly industry GVA in current prices (5206 table
# 45) only starts in September quarter 2002. Before that the current-price
# indicator is back-cast with the growth of the volume indicator times the
# growth of the national GDP implicit price deflator.

build_prod_indicator <- function(a_ind, q_ind, q_nat) {
  shares <- a_ind |>
    group_by(basis, industry, fy) |>
    mutate(share = gva / sum(gva)) |>
    ungroup() |>
    select(state, basis, industry, fy, share)
  if (anyNA(shares$share)) stop_data("missing state industry GVA in 5220 tables 2-9", "ABS 5220.0 tables 2 to 9")
  first_fy <- min(shares$fy)

  one <- function(b) {
    q <- q_ind[q_ind$basis == b, ] |>
      arrange(industry, date) |>
      group_by(industry) |> mutate(gva_lag = lag(gva)) |> ungroup() |>
      mutate(fy = fy_of(date), share_fy = pmax(fy - 1L, first_fy))
    q |>
      inner_join(shares[shares$basis == b, ], by = c("industry", "share_fy" = "fy"), relationship = "many-to-many") |>
      group_by(state, date) |>
      summarise(prod_step = sum(share * gva), denom = sum(share * gva_lag), .groups = "drop") |>
      arrange(state, date) |>
      group_by(state) |>
      mutate(g = prod_step / denom, g = ifelse(row_number() == 1 | is.na(g), 1, g),
             prod = first(prod_step) * cumprod(g)) |>
      ungroup() |>
      select(state, date, prod, prod_step, g)
  }
  cvm <- one("cvm") |> filter(date >= min(q_nat$date))
  cur <- one("current") # starts at the first quarter of table 45

  # back-cast current prices before table 45 begins
  ipd <- q_nat |> select(date, basis, gdp) |> pivot_wider(names_from = basis, values_from = gdp) |>
    arrange(date) |> transmute(date, ipd = current / cvm)
  splice_date <- min(cur$date)
  back <- cvm |> left_join(ipd, by = "date") |> group_by(state) |> arrange(date, .by_group = TRUE) |>
    mutate(nominal_proxy = prod * ipd) |> ungroup() |> select(state, date, nominal_proxy)
  cur_full <- back |> left_join(select(cur, state, date, prod, prod_step), by = c("state", "date")) |>
    group_by(state) |>
    mutate(link = prod[date == splice_date] / nominal_proxy[date == splice_date],
           backcast = date < splice_date,
           prod = ifelse(backcast, nominal_proxy * link, prod),
           prod_step = ifelse(backcast, NA_real_, prod_step)) |>
    ungroup() |> select(state, date, prod, prod_step, backcast)

  bind_rows(
    cvm |> transmute(state, date, basis = "cvm", prod, prod_step, backcast = FALSE),
    cur_full |> mutate(basis = "current")
  ) |> arrange(basis, state, date)
}

# Assemble the quarterly indicator panel used by every method
build_panel <- function(q_state, prod) {
  long <- q_state |>
    pivot_longer(matches("_(cvm|current)$"), names_to = c(".value", "basis"), names_pattern = "(.*)_(cvm|current)$")
  long |>
    inner_join(select(prod, state, date, basis, prod, prod_step), by = c("state", "date", "basis")) |>
    mutate(prodgap = prod - sfd, nx = xg - mg, fy = fy_of(date)) |>
    select(state, basis, date, fy, sfd, prod, prod_step, prodgap, xg, mg, nx, coe, hours) |>
    arrange(basis, state, date)
}

# Bottom-up production indicator ---------------------------------------------------
#
# For each state, industry and price basis: proportional Denton (Cholette,
# first differences) of the state's ANNUAL industry value added on the NATIONAL
# quarterly value added of that industry. Summing the 20 industries gives a
# quarterly series whose four quarters equal the state's annual industry total
# in every benchmark year. After the last benchmark year the state/national ratio
# is held at its last value (that is how Denton extrapolates), so the nowcast
# quarters use no annual information that would not have been available.
#
# Differences from the lagged-share indicator `prod`:
#  * within benchmark years it uses that year's own industry structure;
#  * each industry's state/national ratio moves smoothly between years rather
#    than being reset each September.
#
# It is built from whatever annual data the inputs contain, so a pseudo real-time
# vintage must rebuild it from truncated annual data (vintage_inputs() does).
# Current prices before September 2002 (no national current-price industry data):
# each industry's indicator is its volume series rescaled to join the
# current-price series at September 2002; the annual current-price benchmarks
# then carry prices, interpolated smoothly within each year.
# Taxes less subsidies on products are not included (value added only).

denton_prop <- function(p, a, C) {
  Tn <- length(p); Y <- nrow(C)
  D <- diff(diag(Tn)); Cp <- sweep(C, 2, p, `*`)
  K <- rbind(cbind(2 * crossprod(D), t(Cp)), cbind(Cp, matrix(0, Y, Y)))
  p * solve(K, c(rep(0, Tn), a))[seq_len(Tn)]
}

build_bottom_up <- function(a_ind, q_ind, dates, last_fy) {
  fy <- fy_of(dates)
  purrr::map_dfr(BASES, function(b) {
    nat <- function(bb) {
      w <- q_ind[q_ind$basis == bb, ] |> select(date, industry, gva) |>
        pivot_wider(names_from = industry, values_from = gva) |> arrange(date)
      m <- as.matrix(w[, -1]); rownames(m) <- as.character(w$date); m
    }
    N <- nat(b)
    if (b == "current") {
      V <- nat("cvm")[as.character(dates), colnames(N), drop = FALSE]
      first <- rownames(N)[1]
      link <- N[first, ] / V[first, ]
      N2 <- sweep(V, 2, link, `*`)
      have <- intersect(rownames(N), as.character(dates))
      N2[have, ] <- N[have, ]
      N <- N2
    }
    N <- N[as.character(dates), , drop = FALSE]
    if (anyNA(N)) stop_data("missing national quarterly industry value added", "ABS 5206.0 tables 6 and 45")
    a <- a_ind[a_ind$basis == b & a_ind$fy <= last_fy & a_ind$fy >= min(fy), ]
    years <- sort(unique(a$fy))
    C <- outer(years, fy, "==") * 1
    purrr::map_dfr(STATES, function(s) {
      A <- a[a$state == s, ] |> select(fy, industry, gva) |>
        pivot_wider(names_from = industry, values_from = gva) |> arrange(fy)
      A <- as.matrix(A[, colnames(N)])
      bu <- rowSums(vapply(colnames(N), function(i) denton_prop(N[, i], A[, i], C), numeric(length(dates))))
      tibble::tibble(state = s, basis = b, date = dates, bu = bu)
    })
  })
}

add_bottom_up <- function(panel, a_ind, q_ind, last_fy) {
  bu <- build_bottom_up(a_ind, q_ind, sort(unique(panel$date)), last_fy)
  panel |> select(-any_of(c("bu", "bugap"))) |>
    left_join(bu, by = c("state", "basis", "date")) |>
    mutate(bugap = bu - sfd)
}
