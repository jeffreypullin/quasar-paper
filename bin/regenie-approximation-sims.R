## Model: a specific kind of linear mixed model known as "animal model" by geneticists
## y = mu 1_N + X b + Z u + e = W a + Z u + e
## y is N x 1; X is N x P; Z is N x Q; W is N x (P+1)
## u ~ Norm_Q(0, sigma_u^2 A); e ~ Norm_N(0, sigma^2 I_N)

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(rrBLUP)
  library(lme4)
})

N <- 100
indiv_ids <- sprintf(fmt = paste0("ind%0", floor(log10(N))+1, "i"), 1:N)

X <- cbind(rep(1, N), rnorm(N))
b <- matrix(data = c(4, 2), nrow = 2, ncol = 1)

## ----simul_A-------------------------------------------------------------
Q <- N
n_factors <- 5
n_snps <- 100

snp_ids <- sprintf(fmt = paste0("snp%0", floor(log10(n_snps))+1, "i"), 1:n_snps)
F_mat <- matrix(data=rnorm(n = Q * n_factors), nrow = Q, ncol = n_factors)
M <- mvrnorm(n = n_snps, mu = rep(1, Q), Sigma = F_mat %*% t(F_mat) + diag(Q))
M[M < 0.2] <- 0
M[M >= 0.2 & M <= 1.8] <- 1
M[M > 1.8] <- 2
A <- (1 / n_snps) * t(M) %*% M

Z <- diag(Q)
sigmau2 <- 5
G <- as.matrix(nearPD(sigmau2 * A)$mat)
u <- matrix(mvrnorm(n = 1, mu = rep(0, Q), Sigma = G))

lambda <- 3
sigma2 <- sigmau2 / lambda
R <- sigma2 * diag(N)
e <- matrix(mvrnorm(n = 1, mu = rep(0, N), Sigma = R))

(h2 <- sigmau2 / (sigmau2 + sigma2))

Q <- 100
G <- matrix(0, nrow = 100, ncol = 100)
for (i in 1:100) {
  g <- rbinom(N, 2, 0.2)
  g_cov_adj <- (diag(N) - X %*% solve(t(X) %*% X) %*% t(X)) %*% g
  g_norm <- as.vector(scale(g_cov_adj))
  G[, i] <- g_norm
}
beta <- matrix(rnorm(Q, 0, 0.5), ncol = 1)

# Full model.
y <- X %*% b + G %*% beta + Z %*% u + e

fit_null <- mixed.solve(
  y = y,
  Z = Z,
  K = A,
  X = X,
  method = "REML",
  SE = TRUE,
  return.Hinv = TRUE
)

Omega <- diag(N) * fit_null$Ve + A * fit_null$Vu
Omega_inv <- solve(Omega)

P <- Omega_inv - Omega_inv %*% X %*% solve(t(X) %*% Omega_inv %*% X) %*% t(X) %*% Omega_inv

r <- 0
m <- 100
for (i in 1:m) {
  g <- G[, i]
  r <- r + (t(g) %*% P %*% g) / (t(g) %*% g)
}
rhat <- r / m
rhat

y_tilde <- P %*% y
(t(y_tilde) %*% y_tilde) / (N - 2)

sum(diag(P)) / (N - 2)

