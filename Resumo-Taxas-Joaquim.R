# ============================================================
# DISTRIBUIÇÃO LOG-BILAL REPARAMETRIZADA
# ============================================================

# ------------------------------------------------------------
# 1. FUNÇÃO DE DENSIDADE DE PROBABILIDADE (PDF)
# ------------------------------------------------------------
dLB <- function(y, mu = 0.75) {
  fx1 <- (12 / ((mu / (mu + 24))^(-1/2) - 5)) * 
    y^(4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * 
    (1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
  return(fx1)
}


# ------------------------------------------------------------
# 2. FUNÇÃO DE DISTRIBUIÇÃO ACUMULADA (CDF)
# ------------------------------------------------------------
pLB <- function(y, mu = 0.75) {
  cdf <- 3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 
    2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  return(cdf)
}


# ------------------------------------------------------------
# 3. FUNÇÃO INVERSA (AUXILIAR PARA QUANTIL E GERAÇÃO)
# ------------------------------------------------------------
inverse <- function(f, lower, upper) {
  function(y) {
    uniroot(function(x) f(x) - y, lower = lower, upper = upper, tol = 1e-10)$root
  }
}


# ------------------------------------------------------------
# 4. FUNÇÃO QUANTÍLICA (GERAÇÃO DE VALORES ALEATÓRIOS)
# ------------------------------------------------------------
qLB <- function(u, mu) {
  f <- function(y) {
    3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 
      2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  }
  q <- inverse(f, lower = 0, upper = 1)
  return(q(u))
}


# ------------------------------------------------------------
# 5. TESTE: GERAÇÃO DE NÚMEROS E PLOTAGEM
# ------------------------------------------------------------
set.seed(123)
n <- 1000
mu <- 0.75

u_vals <- runif(n)
y_vals <- sapply(u_vals, qLB, mu = mu)

hist(y_vals, breaks = 30, freq = FALSE, col = "lightblue", 
     main = "Histograma da Log-Bilal com PDF teórica", xlab = "y")

curve(dLB(x, mu = mu), col = "red", lwd = 2, add = TRUE)
legend("topleft", legend = c("PDF teórica"), col = "red", lwd = 1, cex = 0.7)



# ------------------------------------------------------------
# 6. LOG-VEROSSIMILHANÇA
# ------------------------------------------------------------
log_lik_LB <- function(theta, y) {
  mu <- theta[1]
  
  if (mu <= 0 || any(y <= 0 | y >= 1)) {
    return(-Inf)
  }
  
  ll <- log(12 / ((mu / (mu + 24))^(-1/2) - 5)) + 
    (4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * log(y) + 
    log(1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
  
  return(sum(ll))
}


# ------------------------------------------------------------
# Exemplo de estimação via optim
# ------------------------------------------------------------
set.seed(123)
y_vals <- sapply(runif(100), qLB, mu = 0.75)

theta_start <- c(0.5)

optim_result <- optim(theta_start, log_lik_LB, y = y_vals,
                      method = "BFGS",
                      control = list(fnscale = -1))

print(optim_result$par)

# 1. Simulate data Log-Bilal
set.seed(123)
mu_true <- 0.75
y_vals <- sapply(runif(1000), qLB, mu = mu_true)


log_pdf_sum <- sum(log(dLB(y_vals, mu = mu_true)))
log_lik_val <- log_lik_LB(c(mu_true), y_vals)


print(log_pdf_sum)
print(log_lik_val)
print(all.equal(log_pdf_sum, log_lik_val))

