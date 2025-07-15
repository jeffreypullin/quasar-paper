#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(lme4)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(rrBLUP)
  library(dplyr)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

grm <- read_tsv(args[1], show_col_types = FALSE)

N <- 100
K <- as.matrix(grm[1:100, 2:101])
indiv_ids <- colnames(K)
N <- length(indiv_ids)
Q <- N

n_sims <- 1000
rhat_vec <- numeric(n_sims)
regenie_vec <- numeric(n_sims)
quasar_vec <- numeric(n_sims)
for (i in 1:n_sims) {

  X <- cbind(rep(1, N), rnorm(N))
  b <- matrix(data = c(0.2, 1.2), nrow = 2, ncol = 1)

  Z <- diag(Q)
  sigmau2 <- 3
  u <- matrix(mvrnorm(n = 1, mu = rep(0, Q), Sigma = K))

  lambda <- 1
  sigma2 <- sigmau2 / lambda
  R <- sigma2 * diag(N)
  e <- matrix(mvrnorm(n = 1, mu = rep(0, N), Sigma = R))

  G <- matrix(0, nrow = Q, ncol = Q)
  for (k in seq_len(Q)) {
    g <- rbinom(N, 2, 0.2)
    g_cov_adj <- (diag(N) - X %*% solve(t(X) %*% X) %*% t(X)) %*% g
    g_norm <- as.vector(scale(g_cov_adj))
    G[, k] <- g_norm
  }
  beta <- matrix(rnorm(Q, 0, 0.2), ncol = 1)

  y <- X %*% b + G %*% beta + Z %*% u + e

  fit_null <- mixed.solve(
    y = y,
    Z = Z,
    K = K,
    X = X,
    method = "REML",
    SE = TRUE,
    return.Hinv = TRUE
  )

  Omega <- diag(N) * fit_null$Ve + K * fit_null$Vu
  Omega_inv <- solve(Omega)

  P <- Omega_inv - Omega_inv %*% X %*% solve(t(X) %*% Omega_inv %*% X) %*% t(X) %*% Omega_inv

  r <- 0
  m <- 100
  for (j in 1:m) {
    g <- G[, j]
    r <- r + (t(g) %*% P %*% g) / (t(g) %*% g)
  }
  rhat_vec[[i]] <- r / m

  y_tilde <- P %*% y
  regenie_vec[[i]] <- (t(y_tilde) %*% y_tilde) / (N - 2)

  quasar_vec[[i]] <- sum(diag(P)) / (N - 2)
}

lmm_data <- tibble(
  rhat = rhat_vec,
  quasar = quasar_vec,
  regenie = regenie_vec
)

p1 <- lmm_data |>
  ggplot(aes(rhat, quasar)) +
  geom_point() +
  geom_abline(linetype = "dotted") +
  labs(
    x = "Variance-ratio approximation",
    y = "quasar trace-based approximation"
  ) +
  theme_jp()

p2 <- lmm_data |>
  ggplot(aes(quasar, regenie)) +
  geom_point() +
  geom_abline(linetype = "dotted") +
  labs(
    x = "quasar trace-based approximation",
    y = "regenie variance approximation"
  ) +
  theme_jp()

n_sims <- 1000
rhat_vec <- numeric(n_sims)
quasar_first_vec <- numeric(n_sims)
quasar_second_vec <- numeric(n_sims)
for (i in 1:n_sims) {

  X <- cbind(rep(1, N), rnorm(N))
  b <- matrix(data = c(-1, 0.2), nrow = 2, ncol = 1)

  Z <- diag(Q)
  sigmau2 <- 3
  u <- matrix(mvrnorm(n = 1, mu = rep(0, Q), Sigma = K))

  lambda <- 1
  sigma2 <- sigmau2 / lambda
  R <- sigma2 * diag(N)
  e <- matrix(mvrnorm(n = 1, mu = rep(0, N), Sigma = R))

  G <- matrix(0, nrow = Q, ncol = Q)
  G_no_norm <- matrix(0, nrow = Q, ncol = Q)
  for (k in seq_len(Q)) {
    g <- rbinom(N, 2, 0.2)
    g_cov_adj <- (diag(N) - X %*% solve(t(X) %*% X) %*% t(X)) %*% g
    g_norm <- as.vector(scale(g_cov_adj))
    G[, k] <- g_norm
    G_no_norm[, k] <- g
  }
  beta <- matrix(rnorm(Q, 0, 0.2), ncol = 1)

  eta <- X %*% b + G %*% beta + Z %*% u + e
  mu <- exp(eta)
  y <- rpois(N, mu)

  dat <- data.frame(
    fix = X[,2],
    indiv = factor(indiv_ids),
    response = y
  )

  formula <- response ~ 1 + fix + (1 | indiv)
  parsed_formula <- lFormula(
    formula = formula,
    data = dat,
    control = glmerControl(
      check.nobs.vs.nlev = "ignore",
      check.nobs.vs.nRE = "ignore"
    )
  )
  relmat <- list(indiv = K)
  relfac <- relmat
  flist <- parsed_formula$reTrms[["flist"]]
  ztlist <- parsed_formula$reTrms[["Ztlist"]]

  asgn <- attr(flist, "assign")
  for (n in seq_along(relmat)) {
    tn <- which(match(names(relmat)[n], names(flist)) == asgn)
    relmat[[n]] <- Matrix(relmat[[n]], sparse = TRUE)
    relfac[[n]] <- chol(relmat[[n]])
    ztlist[[n]] <- relfac[[n]] %*% ztlist[[n]]
  }
  parsed_formula$reTrms[["Ztlist"]] <- ztlist
  parsed_formula$reTrms[["Zt"]] <- ztlist[[1]]

  deviance_function <- mkGlmerDevfun(
    fr = parsed_formula$fr,
    reTrms = parsed_formula$reTrms,
    X = parsed_formula$X,
    family = poisson()
  )
  optimizer_output <- optimizeGlmer(deviance_function)

  fit_glmm <- mkMerMod(
    rho = environment(deviance_function),
    opt = optimizer_output,
    reTrms = parsed_formula$reTrms,
    fr = parsed_formula$fr
  )

  vc <- as.data.frame(VarCorr(fit_glmm))
  sigma2 <- vc[vc$grp == "indiv", "vcov"]

  mu <- getME(fit_glmm, "mu")
  W <- diag(mu)

  Sigma <- solve(W) + sigma2 * K
  Sigma_inv <- solve(Sigma)
  P <- Sigma_inv - Sigma_inv %*% X %*% solve(t(X) %*% Sigma_inv %*% X) %*% t(X) %*% Sigma_inv
  
  G_w_norm <- matrix(0, nrow = Q, ncol = Q)
  for (k in seq_len(Q)) {
    g <- G_no_norm[, k]
    G_w_norm[, k] <- g - X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% g
  }

  r <- 0
  m <- 100
  for (j in 1:m) {
    g <- G_w_norm[, j]
    r <- r + (t(g) %*% P %*% g) / (t(g) %*% W %*% g)
  }
  rhat_vec[[i]] <- r / m
  P_W <- diag(N) - X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  
  e <- sum(diag(P)) / (sum(diag(W %*% P_W)))
  a <- 2 * sum(diag(P %*% W %*% P_W)) / (sum(diag(W %*% P_W)))^2
  b <- (2 * sum(diag(W %*% P_W %*% W %*% P_W)) * sum(diag(P))) / (sum(diag(W %*% P_W)))^3
      
  tr_W <- sum(diag(W))
  tr_P <- sum(diag(P))
  XtWX_inv <- solve(t(X) %*% W %*% X)
  XtW2X <- t(X) %*% W^2 %*% X
  XtW3X <- t(X) %*% W^3 %*% X
  WX <- W %*% X

  tr_WPw <- tr_W - sum(diag(XtW2X %*% XtWX_inv))
  e1 <- tr_P / tr_WPw

  tr_PWPw <- sum(diag(P %*% W)) - sum(diag(((P %*% WX) %*% XtWX_inv) %*% t(WX)))
  a1 <- 2 * tr_PWPw / tr_WPw^2

  tmp <- sum(diag(XtW2X %*% XtWX_inv %*% XtW2X %*% XtWX_inv))
  tr_WPwWPw <- sum(diag(W^2)) - 2 * sum(diag(XtW3X %*% XtWX_inv)) + tmp
  b1 <- (2 * tr_WPwWPw * tr_P) / tr_WPw^3
  
  quasar_first_vec[[i]] <- e
  quasar_second_vec[[i]] <- e - a + b

}

glmm_data <- tibble(
  rhat = rhat_vec, 
  quasar_first = quasar_first_vec,
  quasar_second = quasar_second_vec
)

p3 <- glmm_data |>
  ggplot(aes(rhat, quasar_first)) +
  geom_point() +
  geom_abline(linetype = "dashed") +
  labs(
    x = "Variance-ratio approximation",
    y = "First-order quasar trace-based approximation"
  ) +
  theme_jp()

p4 <- glmm_data |>
  ggplot(aes(rhat, quasar_second)) +
  geom_point() +
  geom_abline(linetype = "dashed") +
  labs(
    x = "Variance-ratio approximation",
    y = "Second-order quasar trace-based approximation"
  ) +
  theme_jp()

p <- p1 + p2 + p3 + p4 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  "plot-trace-approx-sims.pdf",
  p,
  width = 14,
  height = 10
)

