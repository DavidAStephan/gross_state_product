# Method C: mixed-frequency VAR with stochastic volatility (MF-VAR-SV) -------------
#
# Implements the model of Koop, McIntyre, Mitchell and Poon (2020, Journal of
# Applied Econometrics 35(2), "Regional output growth in the United Kingdom:
# more timely and higher frequency estimates from 1970"), equation numbers below
# refer to that paper.
#
# y_t = (national growth, growth of the 8 states)', log differences x 100.
# National growth is observed every quarter; state growth is latent.
#
#  (1) VAR:      y_t = Phi_0 + Phi_1 y_{t-1} + ... + Phi_7 y_{t-7} + G x_t + u_t,  u_t ~ N(0, Sigma_t)
#                x_t: state-specific exogenous indicators - the equation for state i
#                contains only state i's SFD growth and hours-worked growth.
#  (2) temporal: annual growth of state r = 1/4 y_t + 1/2 y_{t-1} + 3/4 y_{t-2} + y_{t-3}
#                + 3/4 y_{t-4} + 1/2 y_{t-5} + 1/4 y_{t-6}, observed each June quarter, no error.
#  (5) cross-sectional: national growth_t = sum_r w_r y_t^r + eta_t, eta ~ N(0, s2_cs),
#                with a tight inverse-gamma prior on s2_cs.
#  (6)-(8) stochastic volatility: Sigma_t^{-1} = L' D_t^{-1} L, L unit lower triangular,
#                D_t = diag(exp(h_t)), h_t a random walk.
#  (9)-(12) Dirichlet-Laplace prior (alpha = 1/2) on all VAR coefficients, and the
#                same prior on the free elements of L.
# MCMC as in the paper's appendix B: coefficients, L, log-volatilities (Kim-Shephard-
# Chib auxiliary mixture), Dirichlet-Laplace scales, then the latent states.
# Hyperparameters follow appendix B.1.3 (nu_h = 5, S_h = 0.01, V_h = 10).
#
# Where this departs from the paper, and why:
#  * Cross-sectional weights are previous-year GSP shares, not 1/R. The paper's
#    equal weights are a first-order approximation suited to regions of similar
#    size; Australian states range from 1% to 31% of the total.
#  * The four national macro variables (inflation, interest rate, exchange rate,
#    oil price) are not included: this is the paper's 'GVA-only' core plus regional
#    indicators.
#  * Latent states are drawn with a Kalman simulation smoother (KFAS) instead of
#    the precision sampler; both draw from the same conditional distribution.
#  * The prior on s2_cs is rescaled to per cent units: IG(1000, 10), prior mean
#    about 0.01 (s.d. 0.1 percentage point).
#  * Added: because (2) and (5) are log approximations, each retained draw is
#    converted to levels and passed through reconcile_dfm(), so that the brief's
#    two level constraints hold exactly for every draw and for the posterior mean.

suppressPackageStartupMessages({ library(KFAS); library(Matrix) })

MFVAR_DEFAULT <- list(p = 7, n_iter = 6000, burn = 2000, thin = 8, alpha_dl = 0.5,
                      cs_prior = c(1000, 10), nu_h = 5, S_h = 0.01, V_h0 = 10, seed = 20261116)

mfvar_data <- function(pb, inp) {
  Tn <- length(pb$dates); S <- length(STATES)
  dl <- function(M) 100 * rbind(NA_real_, diff(log(M)))
  X1 <- dl(pb$SFD)                                           # own SFD growth
  X2 <- dl(sapply(STATES, function(s) state_X(pb, s, "hours"))) # own hours growth
  ya <- matrix(NA_real_, Tn, S, dimnames = list(NULL, STATES))
  for (i in seq_along(pb$years)[-1]) {
    t <- which(pb$fy == pb$years[i] & fyq_of(pb$dates) == 4)
    ya[t, ] <- 100 * log(pb$GSPa[i, ] / pb$GSPa[i - 1, ])
  }
  ynat <- c(NA_real_, 100 * diff(log(pb$target)))
  cur <- inp$a_state[inp$a_state$basis == "current", ]
  sh <- as.matrix(cur |> select(fy, state, gsp) |> pivot_wider(names_from = state, values_from = gsp) |>
                    arrange(fy) |> select(all_of(STATES)))
  sh <- sh / rowSums(sh); yrs <- sort(unique(cur$fy))
  W <- sh[match(pmin(pmax(pb$fy - 1L, min(yrs)), max(yrs)), yrs), , drop = FALSE]
  list(X1 = X1, X2 = X2, ya = ya, ynat = ynat, W = W, Tn = Tn)
}

# State vector: (y_t, ..., y_{t-6}, 1), y_t of dimension n = 9. Observations:
# 8 annual rows, national growth, and the cross-sectional row.
mfvar_ssm <- function(d) {
  S <- length(STATES); n <- S + 1; L <- 7; m <- n * L + 1; Tn <- d$Tn
  y <- cbind(d$ya, d$ynat, d$ynat)
  Z <- array(0, c(S + 2, m, Tn))
  om <- c(1, 2, 3, 4, 3, 2, 1) / 4
  for (s in seq_len(S)) for (j in 0:6) Z[s, j * n + 1 + s, ] <- om[j + 1]
  Z[S + 1, 1, ] <- 1
  Z[S + 2, 1 + seq_len(S), ] <- t(d$W)
  Tt <- array(0, c(m, m, Tn))
  T1 <- matrix(0, m, m)
  T1[cbind(n + seq_len(n * (L - 1)), seq_len(n * (L - 1)))] <- 1
  T1[m, m] <- 1
  Tt[] <- T1
  R <- matrix(0, m, n); R[cbind(seq_len(n), seq_len(n))] <- 1
  SSModel(y ~ -1 + SSMcustom(Z = Z, T = Tt, R = R, Q = array(diag(n), c(n, n, Tn)), a1 = c(rep(0, m - 1), 1),
                             P1 = diag(c(rep(25, m - 1), 0)), P1inf = matrix(0, m, m), n = Tn),
          H = diag(c(rep(0, S + 1), 0.01)))
}

# --- small samplers ---------------------------------------------------------------
rinvgauss1 <- function(mu) {                 # inverse Gaussian(mean mu, shape 1), vectorised
  w <- mu * stats::rnorm(length(mu))^2 / 2
  x <- pmax(mu / (1 + w + sqrt(2 * w + w^2)), 1e-12)  # smaller root, in a form that does not cancel for large mu
  ifelse(stats::runif(length(mu)) <= mu / (mu + x), x, mu^2 / x)
}
# Dirichlet-Laplace conditional updates (appendix B.17-B.20); returns prior variances
dl_update <- function(coef, st, alpha) {
  k <- length(coef); ab <- pmax(abs(coef), 1e-8)
  psi <- 1 / pmax(rinvgauss1(pmin(st$theta * st$tau / ab, 1e8)), 1e-10)
  tau <- GIGrvg::rgig(1, lambda = k * (alpha - 1), chi = 2 * sum(ab / st$theta), psi = 1)
  Rj <- vapply(ab, function(a) GIGrvg::rgig(1, lambda = alpha - 1, chi = 2 * a, psi = 1), numeric(1))
  theta <- pmax(Rj / sum(Rj), 1e-10)
  list(theta = theta, tau = tau, var = pmin(pmax(psi * theta^2 * tau^2, 1e-10), 1e3))
}
KSC <- list(p = c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750),
            m = c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) - 1.2704,
            v = c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261))
# One log-volatility path given residuals e (Kim, Shephard and Chib 1998)
draw_logvol <- function(e, h, om2, V0) {
  N <- length(e); ys <- log(e^2 + 1e-4)
  lp <- sapply(1:7, function(k) log(KSC$p[k]) + stats::dnorm(ys, h + KSC$m[k], sqrt(KSC$v[k]), log = TRUE))
  pr <- exp(lp - apply(lp, 1, max)); kk <- apply(pr, 1, function(q) sample.int(7, 1, prob = q))
  yy <- ys - KSC$m[kk]; vv <- KSC$v[kk]
  mf <- vf <- numeric(N); a <- 0; P <- V0
  for (t in seq_len(N)) {                      # forward filter (h_0 ~ N(0, V0))
    P <- P + om2; K <- P / (P + vv[t]); a <- a + K * (yy[t] - a); P <- (1 - K) * P
    mf[t] <- a; vf[t] <- P
  }
  hn <- numeric(N); hn[N] <- stats::rnorm(1, mf[N], sqrt(vf[N]))
  for (t in (N - 1):1) {                       # backward sampling
    g <- vf[t] / (vf[t] + om2)
    hn[t] <- stats::rnorm(1, mf[t] + g * (hn[t + 1] - mf[t]), sqrt(vf[t] * (1 - g)))
  }
  g0 <- V0 / (V0 + om2); h0 <- stats::rnorm(1, g0 * hn[1], sqrt(V0 * (1 - g0)))
  list(h = hn, h0 = h0)
}

# --- Gibbs sampler for one price basis ----------------------------------------------
mfvar_gibbs <- function(pb, inp, cfg = MFVAR_DEFAULT) {
  set.seed(cfg$seed)
  S <- length(STATES); n <- S + 1; p <- cfg$p; stopifnot(p == 7)
  d <- mfvar_data(pb, inp); Tn <- d$Tn; N <- Tn - 1           # regression rows are quarters 2..T
  model <- mfvar_ssm(d); m <- dim(model$T)[1]
  rows <- 2:Tn
  ex <- lapply(seq_len(n), function(i) if (i == 1) NULL else cbind(d$X1[rows, i - 1], d$X2[rows, i - 1]))
  kc <- 1 + n * p; ki <- kc + vapply(ex, function(e) if (is.null(e)) 0L else ncol(e), integer(1))
  off <- c(0, cumsum(ki)); ktot <- sum(ki)
  beta <- lapply(ki, function(k) numeric(k))
  Lm <- diag(n); h <- matrix(0, N, n); om2 <- rep(0.05, n); s2cs <- 0.01
  dl_b <- list(theta = rep(1 / ktot, ktot), tau = 1, var = rep(0.01, ktot))  # start tight; the prior scales adapt
  bvec <- numeric(ktot); n_reject <- 0
  na <- n * (n - 1) / 2
  dl_a <- list(theta = rep(1 / na, na), tau = 1, var = rep(0.5, na))
  keep <- seq(cfg$burn + cfg$thin, cfg$n_iter, by = cfg$thin)
  G <- array(NA_real_, c(Tn, S, length(keep)))
  vol <- array(NA_real_, c(N, n, length(keep))); trace <- matrix(NA_real_, length(keep), 3,
                                                               dimnames = list(NULL, c("own_sfd", "own_hours", "sd_cs")))
  kept <- 0
  for (it in seq_len(cfg$n_iter)) {
    # (a) latent states given everything else
    Phi <- t(vapply(beta, function(b) b[2:kc], numeric(kc - 1)))            # n x (n p)
    model$T[seq_len(n), seq_len(n * p), ] <- Phi
    dmat <- vapply(seq_len(n), function(i) beta[[i]][1] + if (i == 1) rep(0, N) else as.numeric(ex[[i]] %*% beta[[i]][kc + 1:2]),
                   numeric(N))                                              # N x n: intercept + exogenous, quarters 2..T
    model$T[seq_len(n), m, seq_len(N)] <- t(dmat)
    Linv <- solve(Lm)
    for (t in seq_len(N)) model$Q[, , t] <- Linv %*% (exp(h[t, ]) * t(Linv))
    model$Q[, , Tn] <- model$Q[, , N]
    model$H[S + 2, S + 2, 1] <- s2cs
    al <- simulateSSM(model, type = "states", nsim = 1)[, , 1]
    Y <- al[rows, seq_len(n), drop = FALSE]
    Zc <- cbind(1, al[rows - 1, seq_len(n * p), drop = FALSE])
    Zi <- lapply(seq_len(n), function(i) cbind(Zc, ex[[i]]))

    # (b) VAR coefficients, drawn jointly (B.9-B.10)
    E <- exp(-h)                                                            # N x n
    K <- matrix(0, ktot, ktot); rhs <- numeric(ktot)
    for (i in seq_len(n)) for (j in i:n) {
      sij <- as.numeric(E %*% (Lm[, i] * Lm[, j]))
      blk <- crossprod(Zi[[i]] * sij, Zi[[j]])
      ii <- off[i] + seq_len(ki[i]); jj <- off[j] + seq_len(ki[j])
      K[ii, jj] <- blk; if (j > i) K[jj, ii] <- t(blk)
      rhs[ii] <- rhs[ii] + crossprod(Zi[[i]], sij * Y[, j])
      if (j > i) rhs[jj] <- rhs[jj] + crossprod(Zi[[j]], sij * Y[, i])
    }
    diag(K) <- diag(K) + 1 / dl_b$var
    U <- chol(K)
    bhat <- backsolve(U, rhs, transpose = TRUE)
    # Draws implying an explosive VAR are rejected (the paper is silent on this; without
    # it an occasional explosive draw sends the latent states to infinity). After 100
    # failed tries the previous coefficients are kept for this sweep.
    for (try_i in 1:100) {
      cand <- backsolve(U, bhat + stats::rnorm(ktot))
      cb <- lapply(seq_len(n), function(i) cand[off[i] + seq_len(ki[i])])
      comp <- rbind(t(vapply(cb, function(b) b[2:kc], numeric(kc - 1))), cbind(diag(n * (p - 1)), matrix(0, n * (p - 1), n)))
      if (max(Mod(eigen(comp, only.values = TRUE)$values)) < 0.999) { bvec <- cand; beta <- cb; break }
      n_reject <- n_reject + 1
    }
    dl_b <- dl_update(bvec, dl_b, cfg$alpha_dl)
    u <- vapply(seq_len(n), function(i) Y[, i] - as.numeric(Zi[[i]] %*% beta[[i]]), numeric(N))

    # (c) free elements of L (B.11-B.12), with their own Dirichlet-Laplace prior
    avec <- numeric(na); pos <- 0
    for (i in 2:n) {
      idx <- pos + seq_len(i - 1); Wm <- -u[, seq_len(i - 1), drop = FALSE]; wt <- exp(-h[, i])
      Ka <- crossprod(Wm * wt, Wm) + diag(1 / dl_a$var[idx], i - 1)
      Ua <- chol(Ka)
      ai <- backsolve(Ua, backsolve(Ua, crossprod(Wm, wt * u[, i]), transpose = TRUE) + stats::rnorm(i - 1))
      Lm[i, seq_len(i - 1)] <- ai; avec[idx] <- ai; pos <- pos + i - 1
    }
    dl_a <- dl_update(avec, dl_a, cfg$alpha_dl)

    # (d) log-volatilities and their innovation variances (B.14-B.16)
    e <- u %*% t(Lm)
    for (i in seq_len(n)) {
      dv <- draw_logvol(e[, i], h[, i], om2[i], cfg$V_h0)
      h[, i] <- dv$h
      om2[i] <- 1 / stats::rgamma(1, cfg$nu_h + N / 2, cfg$S_h + 0.5 * sum(diff(c(dv$h0, dv$h))^2))
    }
    # (e) variance of the cross-sectional restriction
    r <- d$ynat[rows] - rowSums(d$W[rows, ] * Y[, -1])
    s2cs <- 1 / stats::rgamma(1, cfg$cs_prior[1] + N / 2, cfg$cs_prior[2] + 0.5 * sum(r^2))

    if (it %in% keep) {
      kept <- kept + 1; G[, , kept] <- al[, 1 + seq_len(S)]
      Sd <- vapply(seq_len(N), function(t) sqrt(diag(Linv %*% (exp(h[t, ]) * t(Linv)))), numeric(n))
      vol[, , kept] <- t(Sd)
      trace[kept, ] <- c(mean(vapply(2:n, function(i) beta[[i]][kc + 1], 0)), mean(vapply(2:n, function(i) beta[[i]][kc + 2], 0)), sqrt(s2cs))
    }
  }
  list(G = G, vol = vol, trace = trace, d = d, reject_rate = n_reject / (n_reject + cfg$n_iter))
}

# Growth draw -> levels anchored on the first benchmark year
growth_to_levels <- function(g, pb) {
  L <- exp(apply(rbind(0, g[-1, , drop = FALSE]) / 100, 2, cumsum))
  first <- pb$fy == pb$years[1]
  sweep(L, 2, pb$GSPa[1, ] / colSums(L[first, , drop = FALSE]), `*`)
}

method_c_basis <- function(pb, inp, cfg = MFVAR_DEFAULT, reconcile_draws = TRUE) {
  fit <- mfvar_gibbs(pb, inp, cfg)
  nd <- dim(fit$G)[3]
  gsd <- apply(fit$G[-1, , , drop = FALSE], 2, function(x) mean(apply(x, 1, stats::sd)))
  w <- 1 / pmax(gsd, 1e-6)^2
  rec <- function(g) suppressWarnings(reconcile_dfm(growth_to_levels(g, pb), pb$GSPa, pb$target, pb$fy, pb$years, w))
  if (reconcile_draws) {
    Xd <- array(NA_real_, dim(fit$G)); adj <- numeric(nd)
    for (k in seq_len(nd)) { r <- rec(fit$G[, , k]); Xd[, , k] <- r$X; adj[k] <- mean(abs(r$ratio - 1)) }
    est <- apply(Xd, c(1, 2), mean)
    q <- function(pr) apply(Xd, c(1, 2), stats::quantile, probs = pr, names = FALSE)
    bands <- list(lower_68 = q(0.16), upper_68 = q(0.84), lower_90 = q(0.05), upper_90 = q(0.95))
  } else {
    r <- rec(apply(fit$G, c(1, 2), mean)); est <- r$X; adj <- mean(abs(r$ratio - 1)); bands <- NULL
  }
  dimnames(est) <- dimnames(pb$SFD)
  out <- mat_to_long(est, pb$dates, pb$basis)
  if (!is.null(bands)) for (nm in names(bands)) {
    M <- bands[[nm]]; dimnames(M) <- dimnames(pb$SFD)
    out <- left_join(out, mat_to_long(M, pb$dates, pb$basis, nm), by = c("date", "state", "basis"))
  }
  out$benchmarked <- rep(pb$bench, times = length(STATES))
  half <- seq_len(nd) <= nd / 2
  diag <- tibble::tibble(basis = pb$basis, state = STATES, draws = nd, post_sd_growth = gsd,
                         # crude convergence check: mean growth in the first vs second half of the kept draws
                         halves_max_abs_diff = apply(abs(apply(fit$G[-1, , half, drop = FALSE], c(1, 2), mean) -
                                                         apply(fit$G[-1, , !half, drop = FALSE], c(1, 2), mean)), 2, max),
                         mean_abs_recon_adj_pct = 100 * mean(adj), own_sfd_coef = mean(fit$trace[, "own_sfd"]),
                         own_hours_coef = mean(fit$trace[, "own_hours"]), sd_cross_section = mean(fit$trace[, "sd_cs"]),
                         explosive_draws_rejected = fit$reject_rate)
  volm <- apply(fit$vol, c(1, 2), mean); colnames(volm) <- c("AUS", STATES)
  vol <- tibble::tibble(date = rep(pb$dates[-1], times = ncol(volm)), series = rep(colnames(volm), each = nrow(volm)),
                        basis = pb$basis, shock_sd = as.numeric(volm))
  list(est = out, diag = diag, vol = vol)
}

method_c <- function(inp, cfg = MFVAR_DEFAULT, reconcile_draws = TRUE, cores = 2) {
  out <- parallel::mclapply(BASES, function(b) method_c_basis(prep_basis(inp, b), inp, cfg, reconcile_draws), mc.cores = cores)
  bad <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("Method C failed: ", out[[which(bad)[1]]])
  list(est = purrr::map_dfr(out, "est"), diag = purrr::map_dfr(out, "diag"), vol = purrr::map_dfr(out, "vol"))
}
