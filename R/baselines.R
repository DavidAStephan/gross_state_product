# Baselines (for comparison only) ------------------------------------------------
# Both produce a preliminary quarterly GSP from SFD alone and then scale every
# quarter proportionally so the states sum to the national target. The scaling
# means the cross-sectional constraint holds exactly but the temporal one only
# approximately - these are yardsticks, not candidates.

scale_to_national <- function(P, target) P * (target / rowSums(P))

# Pro-rata: GSP_q = SFD_q x (annual GSP / annual SFD); last ratio carried forward
baseline_prorata <- function(pb) {
  ratio <- pb$GSPa / (pb$C %*% pb$SFD)
  idx <- match(pmin(pb$fy, max(pb$years)), pb$years)
  P <- pb$SFD * ratio[idx, , drop = FALSE]
  list(prelim = P, estimate = scale_to_national(P, pb$target))
}

# Proportional Denton (Cholette variant) of annual GSP on quarterly SFD
baseline_denton <- function(pb) {
  P <- sapply(STATES, function(s) {
    hf <- stats::ts(pb$SFD[, s], start = c(pb$fy[1], 1), frequency = 4)
    lf <- stats::ts(pb$GSPa[, s], start = pb$years[1], frequency = 1)
    as.numeric(stats::predict(tempdisagg::td(lf ~ 0 + hf, method = "denton-cholette", conversion = "sum")))
  })
  dimnames(P) <- dimnames(pb$SFD)
  list(prelim = P, estimate = scale_to_national(P, pb$target))
}

run_baselines <- function(inp) {
  purrr::map_dfr(BASES, function(b) {
    pb <- prep_basis(inp, b)
    bind_rows(
      mat_to_long(baseline_prorata(pb)$estimate, pb$dates, b) |> mutate(method = "baseline_prorata"),
      mat_to_long(baseline_denton(pb)$estimate, pb$dates, b) |> mutate(method = "baseline_denton")
    )
  })
}
