# Method A: indicator regression + movement-preserving reconciliation ------------

suppressPackageStartupMessages({
  library(tempdisagg)
  library(Matrix)
})

# Candidate indicator sets for the gap (GSP - SFD). All are quarterly, in the
# same price basis as the gap. prodgap = production indicator - SFD.
SPECS <- list(
  s1_prodgap        = "prodgap",
  s2_prodgap_nx     = c("prodgap", "nx"),
  s3_prodgap_x_m    = c("prodgap", "xg", "mg"),
  s4_prodgap_nx_coe = c("prodgap", "nx", "coe"),
  s5_prodgap_nx_hrs = c("prodgap", "nx", "hours"),
  s6_prod_sfd_nx    = c("prod", "sfd", "nx"),
  # bottom-up industry indicator (bu = sum of state industries benchmarked to annual
  # state industry value added; bugap = bu - SFD). See build_bottom_up().
  s7_bugap          = "bugap",
  s8_bugap_nx       = c("bugap", "nx"),
  s9_bugap_nx_hrs   = c("bugap", "nx", "hours"),
  s10_bu_sfd_nx     = c("bu", "sfd", "nx"),
  s11_bugap_prodgap_nx = c("bugap", "prodgap", "nx")
)
# specifications that do not use the bottom-up indicator (for comparisons)
SPECS_NO_BU <- c("s1_prodgap", "s2_prodgap_nx", "s3_prodgap_x_m", "s4_prodgap_nx_coe", "s5_prodgap_nx_hrs", "s6_prod_sfd_nx")
TD_METHODS <- c(fernandez = "fernandez", chowlin = "chow-lin-maxlog", litterman = "litterman-maxlog")
POOLED_SPECS <- c("s1_prodgap", "s2_prodgap_nx", "s7_bugap", "s8_bugap_nx")

model_grid <- function() {
  g <- expand.grid(spec = names(SPECS), td = names(TD_METHODS), stringsAsFactors = FALSE)
  g$model <- paste(g$td, g$spec, sep = ":")
  rbind(g, data.frame(spec = POOLED_SPECS, td = "pooled_fernandez",
                      model = paste("pooled_fernandez", POOLED_SPECS, sep = ":")))
}

state_X <- function(pb, s, vars) {
  d <- pb$panel[pb$panel$state == s, ]
  d <- d[match(pb$dates, d$date), ]
  as.matrix(d[, vars, drop = FALSE])
}

# --- Single-state temporal disaggregation (tempdisagg) --------------------------
# Returns the quarterly preliminary gap (extrapolated over unbenchmarked
# quarters), the regression part, and the fit.
td_state <- function(pb, s, vars, method) {
  X <- state_X(pb, s, vars)
  env <- new.env()
  env$y <- stats::ts(pb$gap_a[, s], start = pb$years[1], frequency = 1)
  for (v in vars) env[[v]] <- stats::ts(X[, v] / mean(abs(X[, v])), start = c(pb$fy[1], 1), frequency = 4)
  f <- stats::as.formula(paste("y ~", paste(vars, collapse = " + ")), env = env)
  fit <- tempdisagg::td(f, method = method, conversion = "sum")
  gap <- as.numeric(stats::predict(fit))
  list(gap = gap, reg = as.numeric(fit$p), fit = fit)
}

# --- Pooled Fernandez: common slopes across states -------------------------------
# gap_s / k_s = a_s + b' (indicators_s / k_s) + u_s,  u_s a random walk,
# k_s = the state's mean quarterly GSP (a constant, so annual sums are
# preserved). State intercepts, common slopes b, state-specific innovation
# variances (two-step feasible GLS).
fernandez_sigma <- function(n) {
  D <- diag(n); D[cbind(2:n, 1:(n - 1))] <- -1
  solve(crossprod(D))
}
td_pooled <- function(pb, vars) {
  S <- length(STATES); Tn <- length(pb$dates); Y <- length(pb$years)
  C <- pb$C
  Sig <- fernandez_sigma(Tn)
  Vinv <- solve(C %*% Sig %*% t(C))
  k <- colMeans(pb$GSPa) / 4
  Xs <- lapply(STATES, function(s) state_X(pb, s, vars) / k[[s]])
  ys <- lapply(STATES, function(s) pb$gap_a[, s] / k[[s]])
  K <- length(vars)
  build <- function(i) { # design: S intercept columns + K common slopes
    Z <- matrix(0, Tn, S + K); Z[, i] <- 1; Z[, S + seq_len(K)] <- Xs[[i]]; Z
  }
  Zs <- lapply(seq_len(S), build)
  sig2 <- rep(1, S)
  for (iter in 1:3) {
    A <- matrix(0, S + K, S + K); bvec <- numeric(S + K)
    for (i in seq_len(S)) {
      Za <- C %*% Zs[[i]]
      A <- A + crossprod(Za, Vinv %*% Za) / sig2[i]
      bvec <- bvec + crossprod(Za, Vinv %*% ys[[i]]) / sig2[i]
    }
    beta <- solve(A, bvec)
    for (i in seq_len(S)) {
      e <- ys[[i]] - C %*% Zs[[i]] %*% beta
      sig2[i] <- as.numeric(crossprod(e, Vinv %*% e)) / Y
    }
  }
  L <- Sig %*% t(C) %*% Vinv
  gap <- reg <- matrix(NA_real_, Tn, S, dimnames = list(NULL, STATES))
  for (i in seq_len(S)) {
    p <- as.numeric(Zs[[i]] %*% beta)
    e <- ys[[i]] - C %*% p
    reg[, i] <- p * k[i]
    gap[, i] <- (p + as.numeric(L %*% e)) * k[i]
  }
  list(gap = gap, reg = reg, beta = setNames(as.numeric(beta)[S + seq_len(K)], vars), sig2 = sig2)
}

# --- Preliminary gaps for one basis under a given model choice -------------------
# models: named character vector (names = STATES) of model labels from model_grid()
prelim_gaps <- function(pb, models, keep_fits = FALSE) {
  Tn <- length(pb$dates)
  gap <- reg <- matrix(NA_real_, Tn, length(STATES), dimnames = list(NULL, STATES))
  fits <- list(); pooled_cache <- list()
  for (s in STATES) {
    parts <- strsplit(models[[s]], ":")[[1]]
    vars <- SPECS[[parts[2]]]
    if (parts[1] == "pooled_fernandez") {
      if (is.null(pooled_cache[[parts[2]]])) pooled_cache[[parts[2]]] <- td_pooled(pb, vars)
      gap[, s] <- pooled_cache[[parts[2]]]$gap[, s]; reg[, s] <- pooled_cache[[parts[2]]]$reg[, s]
      if (keep_fits) fits[[s]] <- pooled_cache[[parts[2]]]
    } else {
      r <- td_state(pb, s, vars, TD_METHODS[[parts[1]]])
      gap[, s] <- r$gap; reg[, s] <- r$reg
      if (keep_fits) fits[[s]] <- r$fit
    }
  }
  list(gap = gap, reg = reg, fits = fits)
}

# --- Out-of-sample errors of the PRELIMINARY estimates, all candidate models -----
# For target year t: annual GSP known to t-1, quarterly data to June t. Error =
# 100 x (implied annual GSP_t - actual) / actual. Used for model selection.
oos_prelim_errors <- function(inp, target_years) {
  grid <- model_grid()
  actual <- inp$a_state
  purrr::map_dfr(target_years, function(t) {
    vin <- vintage_inputs(inp, t - 1)
    purrr::map_dfr(BASES, function(b) {
      pb <- prep_basis(vin, b)
      in_t <- pb$fy == t
      act <- actual[actual$basis == b & actual$fy == t, ]
      act <- act$gsp[match(STATES, act$state)]
      sfd_t <- colSums(pb$SFD[in_t, , drop = FALSE])
      pooled <- lapply(setNames(POOLED_SPECS, POOLED_SPECS), function(sp)
        tryCatch(td_pooled(pb, SPECS[[sp]])$gap, error = function(e) NULL))
      purrr::map_dfr(seq_len(nrow(grid)), function(i) {
        g <- grid[i, ]
        fc <- vapply(STATES, function(s) {
          gap <- if (g$td == "pooled_fernandez") pooled[[g$spec]][, s] else
            tryCatch(td_state(pb, s, SPECS[[g$spec]], TD_METHODS[[g$td]])$gap, error = function(e) NULL)
          if (is.null(gap)) NA_real_ else sum(gap[in_t])
        }, numeric(1))
        tibble::tibble(basis = b, state = STATES, model = g$model, fy = t,
                       err_pct = 100 * (sfd_t + fc - act) / act)
      })
    })
  })
}

# Choose, per state and basis, the model with the lowest RMSE over the `window`
# target years ending at `last_fy` (i.e. using only errors observable at that
# vintage).
# Which indicator sets compete in model selection. The bottom-up sets are tested
# (04b_bottom_up_test.R) but excluded by default because they did not improve the
# pseudo real-time errors; set GSP_USE_BU=true to let them compete.
default_specs <- function() {
  if (identical(tolower(Sys.getenv("GSP_USE_BU", "false")), "true")) names(SPECS) else SPECS_NO_BU
}
select_models <- function(oos, last_fy, window = 10, specs = default_specs()) {
  oos |>
    filter(sub("^[^:]+:", "", model) %in% specs) |>
    filter(fy <= last_fy, fy > last_fy - window) |>
    group_by(basis, state, model) |>
    summarise(rmse = sqrt(mean(err_pct^2)), mae = mean(abs(err_pct)), n = n(), ok = !anyNA(err_pct), .groups = "drop") |>
    filter(ok) |>
    group_by(basis, state) |>
    mutate(rank = rank(rmse, ties.method = "first"), n_models = n()) |>
    ungroup()
}

# --- Two-way reconciliation (Di Fonzo & Marini, proportional first differences) --
# Find X (T x S), close to preliminary P in quarterly movements, such that
#   temporal:        sum over the 4 quarters of each benchmark year of X[, s] = GSPa[y, s]
#   cross-sectional: sum over states of X[t, ] = target[t]   (every quarter)
# by minimising   sum_s w_s sum_t ( X[t,s]/P[t,s] - X[t-1,s]/P[t-1,s] )^2 .
# Working variable r = X / P (S*T vector, stacked by state).
#
# Rank deficiency: within each benchmark year the four cross-sectional rows and
# the S temporal rows add up to the same total, so one constraint per benchmark
# year is redundant. After checking the totals really do agree, the last
# state's temporal row is dropped for each benchmark year; the dropped rows are
# then verified on the solution.
reconcile_dfm <- function(P, GSPa, target, fy, years, w = rep(1, ncol(P)), tol = 1e-9) {
  Tn <- nrow(P); S <- ncol(P); Y <- length(years)
  stopifnot(all(P > 0), length(target) == Tn, nrow(GSPa) == Y, length(w) == S, all(w > 0))
  C <- outer(years, fy, "==") * 1
  incons <- max(abs(rowSums(GSPa) - as.numeric(C %*% target)) / rowSums(GSPa))
  if (incons > tol) stop("Constraint sets are mutually inconsistent (max relative gap ", signif(incons, 3), ")")

  w <- w / mean(w)
  D <- Matrix::bandSparse(Tn - 1, Tn, k = c(0, 1), diagonals = list(rep(-1, Tn - 1), rep(1, Tn - 1)))
  DtD <- Matrix::crossprod(D)
  Q <- Matrix::bdiag(lapply(seq_len(S), function(s) w[s] * DtD))

  # temporal rows (scaled so each right-hand side is 1)
  tr <- expand.grid(y = seq_len(Y), s = seq_len(S))
  ti <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    tt <- which(fy == years[tr$y[i]])
    cbind(i, (tr$s[i] - 1) * Tn + tt, P[tt, tr$s[i]] / GSPa[tr$y[i], tr$s[i]])
  }))
  A_temp <- Matrix::sparseMatrix(i = ti[, 1], j = ti[, 2], x = ti[, 3], dims = c(nrow(tr), S * Tn))
  # cross-sectional rows
  A_cross <- Matrix::sparseMatrix(i = rep(seq_len(Tn), S), j = seq_len(S * Tn),
                                  x = as.numeric(P / target), dims = c(Tn, S * Tn))
  keep <- tr$s != S # drop the redundant row: last state, every benchmark year
  A <- rbind(A_temp[keep, , drop = FALSE], A_cross)
  rk <- as.integer(Matrix::rankMatrix(A, method = "qr"))
  if (rk != nrow(A)) stop("Constraint matrix still rank deficient after dropping redundant rows: rank ", rk, " of ", nrow(A))

  m <- nrow(A); n <- S * Tn
  KKT <- rbind(cbind(2 * Q, Matrix::t(A)), cbind(A, Matrix::Matrix(0, m, m, sparse = TRUE)))
  rhs <- c(rep(0, n), rep(1, m))
  sol <- as.numeric(Matrix::solve(KKT, rhs))
  # one step of iterative refinement
  sol <- sol + as.numeric(Matrix::solve(KKT, rhs - as.numeric(KKT %*% sol)))
  r <- matrix(sol[seq_len(n)], Tn, S)
  X <- P * r
  dimnames(X) <- dimnames(P)

  e_temp <- max(abs(C %*% X - GSPa) / GSPa)          # includes the dropped rows
  e_cross <- max(abs(rowSums(X) - target) / target)
  if (max(e_temp, e_cross) > 1e-8) stop("Reconciliation failed constraint check: ", signif(max(e_temp, e_cross), 3))
  list(X = X, ratio = r, n_constraints = nrow(A_temp) + Tn, n_dropped = sum(!keep), rank = rk,
       pre_inconsistency = incons, err_temporal = e_temp, err_cross = e_cross)
}

# Reconciliation weights: a state whose preliminary estimate is less reliable
# should absorb more of the adjustment. Reliability = variance of the quarterly
# change in the distributed regression residual (preliminary gap - regression
# part), in $m. In the proportional criterion a weight of
#   w_s = (mean preliminary level_s)^2 / var_s
# makes the LEVEL adjustment a state absorbs proportional to var_s.
recon_weights <- function(P, gap, reg) {
  v <- apply(gap - reg, 2, function(u) stats::var(diff(u)))
  v <- pmax(v, 1e-8 * colMeans(P)^2)
  list(w = colMeans(P)^2 / v, resid_sd = sqrt(v), rel_sd_pct = 100 * sqrt(v) / colMeans(P))
}

# --- Full Method A for one set of inputs ----------------------------------------
# models: tibble(basis, state, model). Returns long estimates plus diagnostics.
method_a <- function(inp, models, keep_fits = FALSE) {
  out <- lapply(setNames(BASES, BASES), function(b) {
    pb <- prep_basis(inp, b)
    mb <- models[models$basis == b, ]
    mvec <- setNames(mb$model[match(STATES, mb$state)], STATES)
    pg <- prelim_gaps(pb, mvec, keep_fits)
    P <- pb$SFD + pg$gap
    wt <- recon_weights(P, pg$gap, pg$reg)
    rec <- reconcile_dfm(P, pb$GSPa, pb$target, pb$fy, pb$years, wt$w)
    est <- mat_to_long(rec$X, pb$dates, b) |>
      left_join(mat_to_long(P, pb$dates, b, "prelim"), by = c("date", "state", "basis"))
    list(est = est, fits = pg$fits, weights = tibble::tibble(basis = b, state = STATES, model = mvec,
                                                             resid_sd = wt$resid_sd, rel_sd_pct = wt$rel_sd_pct),
         recon = tibble::tibble(basis = b, n_constraints = rec$n_constraints, n_dropped = rec$n_dropped,
                                rank = rec$rank, pre_inconsistency = rec$pre_inconsistency,
                                err_temporal = rec$err_temporal, err_cross = rec$err_cross))
  })
  list(est = purrr::map_dfr(out, "est"), weights = purrr::map_dfr(out, "weights"),
       recon = purrr::map_dfr(out, "recon"), fits = lapply(out, `[[`, "fits"))
}
