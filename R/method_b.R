# Method B: joint linear Gaussian state-space model (KFAS) ------------------------
#
# Everything is modelled in LEVELS, expressed as per cent of the state's own
# quarterly SFD so that the model is homoskedastic over 35 years of growth:
#     g[s,t] = 100 * gap[s,t] / SFD[s,t]              (the unknown)
#     z1[s,t] = 100 * (production indicator - SFD) / SFD
#     z2[s,t] = 100 * (goods exports - goods imports) / SFD
#
# State equations (per state s)
#     g[s,t]  = g[s,t-1] + lambda_s f_t + eta[s,t]     local level + common factor f_t ~ N(0,1)
#     m1[s,t] = m1[s,t-1] + zeta1[s,t]                  persistent error of indicator 1
#     m2[s,t] = m2[s,t-1] + zeta2[s,t]                  persistent error of indicator 2
#     c[s,t]  = psi_t c[s,t-1] + a[s,t] g[s,t]          cumulator; psi_t = 0 in September quarters,
#                                                       a[s,t] = SFD[s,t] / annual SFD[s,fy]
# Measurement equations
#     z1[s,t] = b1_s g[s,t] + m1[s,t] + e1[s,t]         indicators load on the gap
#     z2[s,t] = b2_s g[s,t] + m2[s,t] + e2[s,t]
#     100 * annual gap / annual SFD = c[s,t]            June quarters of benchmark years, NO error
#     100 * national gap_t / sum SFD_t = sum_s omega[s,t] g[s,t]   every quarter, NO error
# (omega = state share of summed SFD). Because the last two have zero
# measurement error, the smoothed g satisfy both accounting constraints exactly.
#
# Estimation is two-stage maximum likelihood:
#   1. state by state (no cross-sectional row): b1, b2 and the five standard deviations;
#   2. jointly: the common-factor loadings, holding each state's total gap
#      innovation variance at its stage-1 value (lambda_s = tau_s tanh(theta_s)).
# Intervals are conditional on these estimated parameters and treat SFD,
# national GDP and the annual benchmarks as known.

suppressPackageStartupMessages(library(KFAS))

# Data for the state-space model, one price basis. A dummy first period (all
# observations missing) is prepended so the cumulator can start at exactly zero.
# First indicator of the gap: "prodgap" (lagged-share production indicator) or
# "bugap" (bottom-up industry indicator). Chosen by the comparison in
# 04b_bottom_up_test.R; see output/bottom_up_test_*.csv.
MB_IND1 <- "prodgap"

ssm_data <- function(pb, ind1 = MB_IND1) {
  Tn <- length(pb$dates)
  Z1 <- 100 * sapply(STATES, function(s) state_X(pb, s, ind1)) / pb$SFD
  Z2 <- 100 * sapply(STATES, function(s) state_X(pb, s, "nx")) / pb$SFD
  yr <- sort(unique(pb$fy))
  sfd_fy <- rowsum(pb$SFD, pb$fy)[as.character(pb$fy), , drop = FALSE]
  A <- pb$SFD / sfd_fy                       # cumulator weights (sum to 1 within a full year)
  YA <- matrix(NA_real_, Tn, length(STATES), dimnames = list(NULL, STATES))
  june_b <- which(fyq_of(pb$dates) == 4 & pb$bench)
  YA[june_b, ] <- 100 * pb$gap_a[as.character(pb$fy[june_b]), ] / (pb$C %*% pb$SFD)[as.character(pb$fy[june_b]), ]
  # Rank deficiency: in the June quarter of a benchmark year the cross-sectional
  # observation is implied exactly by the eight annual observations and the three
  # earlier cross-sectional observations of that year. A redundant zero-error
  # observation makes the likelihood numerically erratic, so it is set to
  # missing; the smoothed estimates still satisfy it (checked in 05_method_b.R).
  ycross <- 100 * pb$natgap / rowSums(pb$SFD)
  ycross[june_b] <- NA_real_
  list(Z1 = Z1, Z2 = Z2, A = A, YA = YA, psi = as.numeric(fyq_of(pb$dates) != 1),
       omega = pb$SFD / rowSums(pb$SFD), ycross = ycross, Tn = Tn)
}

# Build the KFAS model for a subset of states `idx` (1 state in stage 1, all 8
# in stage 2). par: list(b1, b2, sd_eta, lam, sd_m1, sd_m2, sd_e1, sd_e2).
build_ssm <- function(d, idx, par, cross = TRUE) {
  S <- length(idx); n <- d$Tn + 1; m <- 4 * S; k <- 3 * S
  ig <- seq_len(S); ic <- S + ig; im1 <- 2 * S + ig; im2 <- 3 * S + ig
  pad <- function(M) rbind(NA_real_, M[, idx, drop = FALSE])
  y <- cbind(pad(d$Z1), pad(d$Z2), pad(d$YA))
  if (cross) y <- cbind(y, c(NA_real_, d$ycross))
  p <- ncol(y)
  A <- rbind(d$A[, idx, drop = FALSE], 0); psi <- c(d$psi, 0) # value used in transition j -> j+1

  Z <- array(0, c(p, m, n))
  for (i in ig) {
    Z[i, ig[i], ] <- par$b1[i];     Z[i, im1[i], ] <- 1
    Z[S + i, ig[i], ] <- par$b2[i]; Z[S + i, im2[i], ] <- 1
    Z[2 * S + i, ic[i], ] <- 1
  }
  if (cross) Z[3 * S + 1, ig, ] <- t(rbind(0, d$omega[, idx, drop = FALSE]))
  Tt <- array(0, c(m, m, n)); Rt <- array(0, c(m, k, n))
  for (i in ig) {
    Tt[ig[i], ig[i], ] <- 1; Tt[im1[i], im1[i], ] <- 1; Tt[im2[i], im2[i], ] <- 1
    Tt[ic[i], ic[i], ] <- psi; Tt[ic[i], ig[i], ] <- A[, i]
    Rt[ig[i], i, ] <- 1; Rt[ic[i], i, ] <- A[, i]
    Rt[im1[i], S + i, ] <- 1; Rt[im2[i], 2 * S + i, ] <- 1
  }
  P1inf <- diag(m); P1inf[cbind(ic, ic)] <- 0
  model <- SSModel(y ~ -1 + SSMcustom(Z = Z, T = Tt, R = Rt, Q = diag(k), a1 = rep(0, m),
                                      P1 = matrix(0, m, m), P1inf = P1inf, n = n),
                   H = diag(p))
  model <- set_ssm_par(model, par, S, cross)
  attr(model, "layout") <- list(S = S, ig = ig, ic = ic, cross = cross)
  model
}

set_ssm_par <- function(model, par, S, cross) {
  ig <- seq_len(S)
  for (i in ig) { model$Z[i, i, ] <- par$b1[i]; model$Z[S + i, i, ] <- par$b2[i] }
  Q <- matrix(0, 3 * S, 3 * S)
  Q[ig, ig] <- diag(par$sd_eta^2, S) + tcrossprod(par$lam)
  diag(Q)[S + ig] <- par$sd_m1^2; diag(Q)[2 * S + ig] <- par$sd_m2^2
  model$Q[, , 1] <- Q
  h <- c(par$sd_e1^2, par$sd_e2^2, rep(0, S), if (cross) 0)
  model$H[, , 1] <- diag(h, length(h))
  model
}

safe_loglik <- function(model) {
  ll <- tryCatch(logLik(model), error = function(e) -Inf)
  if (!is.finite(ll)) -1e10 else ll
}

# Stage 1: one state at a time
fit_stage1 <- function(d, i) {
  # starting values: slope of annual changes in the indicator on annual changes in the gap
  ya <- d$YA[, i]; jn <- which(!is.na(ya))
  za1 <- sapply(jn, function(j) sum(d$A[(j - 3):j, i] * d$Z1[(j - 3):j, i]))
  za2 <- sapply(jn, function(j) sum(d$A[(j - 3):j, i] * d$Z2[(j - 3):j, i]))
  b10 <- stats::coef(stats::lm(diff(za1) ~ 0 + diff(ya[jn])))[[1]]
  b20 <- stats::coef(stats::lm(diff(za2) ~ 0 + diff(ya[jn])))[[1]]
  s_g <- stats::sd(diff(ya[jn])) / 2
  start <- c(b10, b20, log(c(s_g, stats::sd(diff(d$Z1[, i])), stats::sd(diff(d$Z2[, i])),
                            stats::sd(diff(d$Z1[, i])) / 2, stats::sd(diff(d$Z2[, i])) / 2)))
  topar <- function(th) list(b1 = th[1], b2 = th[2], sd_eta = exp(th[3]), lam = 0,
                             sd_m1 = exp(th[4]), sd_m2 = exp(th[5]), sd_e1 = exp(th[6]), sd_e2 = exp(th[7]))
  model <- build_ssm(d, i, topar(start), cross = FALSE)
  nll <- function(th) {
    if (any(abs(th[3:7]) > 15)) return(1e10)
    -safe_loglik(set_ssm_par(model, topar(th), 1, FALSE))
  }
  opt <- stats::optim(start, nll, method = "BFGS", control = list(maxit = 300, reltol = 1e-9))
  c(topar(opt$par), list(loglik = -opt$value, convergence = opt$convergence))
}

# Stage 2: common-factor loadings in the joint model
fit_stage2 <- function(d, s1) {
  S <- length(STATES)
  base <- list(b1 = sapply(s1, `[[`, "b1"), b2 = sapply(s1, `[[`, "b2"),
               sd_m1 = sapply(s1, `[[`, "sd_m1"), sd_m2 = sapply(s1, `[[`, "sd_m2"),
               sd_e1 = sapply(s1, `[[`, "sd_e1"), sd_e2 = sapply(s1, `[[`, "sd_e2"))
  tau <- sapply(s1, `[[`, "sd_eta")
  topar <- function(th) c(base, list(lam = tau * tanh(th), sd_eta = tau * sqrt(1 - tanh(th)^2)))
  model <- build_ssm(d, seq_len(S), topar(rep(0.3, S)), cross = TRUE)
  nll <- function(th) -safe_loglik(set_ssm_par(model, topar(th), S, TRUE))
  opt <- stats::optim(rep(0.3, S), nll, method = "BFGS", control = list(maxit = 200, reltol = 1e-9))
  ll0 <- -nll(rep(0, S))
  list(par = topar(opt$par), model = set_ssm_par(model, topar(opt$par), S, TRUE),
       loglik = -opt$value, loglik_no_factor = ll0, convergence = opt$convergence, tau = tau)
}

# Full Method B for one basis
method_b_basis <- function(pb, ind1 = MB_IND1) {
  d <- ssm_data(pb, ind1)
  s1 <- lapply(seq_along(STATES), function(i) fit_stage1(d, i))
  s2 <- fit_stage2(d, s1)
  ks <- KFS(s2$model, filtering = "none", smoothing = "state")
  L <- attr(s2$model, "layout")
  g <- ks$alphahat[-1, L$ig, drop = FALSE]
  v <- t(apply(ks$V[L$ig, L$ig, -1, drop = FALSE], 3, diag))
  se <- sqrt(pmax(v, 0))
  X <- pb$SFD * (1 + g / 100); SE <- pb$SFD * se / 100
  dimnames(X) <- dimnames(SE) <- dimnames(pb$SFD)
  pars <- tibble::tibble(
    basis = pb$basis, state = STATES,
    ind1 = ind1, b_ind1 = s2$par$b1, b_nx = s2$par$b2, sd_gap_innov = s2$tau,
    factor_loading = s2$par$lam, factor_share = (s2$par$lam / s2$tau)^2,
    sd_m_ind1 = s2$par$sd_m1, sd_m_nx = s2$par$sd_m2, sd_e_ind1 = s2$par$sd_e1, sd_e_nx = s2$par$sd_e2,
    stage1_convergence = sapply(s1, `[[`, "convergence")
  )
  fit <- tibble::tibble(basis = pb$basis, loglik = s2$loglik, loglik_no_factor = s2$loglik_no_factor,
                        stage2_convergence = s2$convergence)
  est <- mat_to_long(X, pb$dates, pb$basis) |>
    left_join(mat_to_long(SE, pb$dates, pb$basis, "se"), by = c("date", "state", "basis")) |>
    mutate(lower_68 = estimate - stats::qnorm(0.84) * se, upper_68 = estimate + stats::qnorm(0.84) * se,
           lower_90 = estimate - stats::qnorm(0.95) * se, upper_90 = estimate + stats::qnorm(0.95) * se,
           benchmarked = rep(pb$bench, times = length(STATES)))
  list(est = est, pars = pars, fit = fit)
}

method_b <- function(inp, ind1 = MB_IND1) {
  out <- lapply(BASES, function(b) method_b_basis(prep_basis(inp, b), ind1))
  list(est = purrr::map_dfr(out, "est"), pars = purrr::map_dfr(out, "pars"), fit = purrr::map_dfr(out, "fit"))
}
