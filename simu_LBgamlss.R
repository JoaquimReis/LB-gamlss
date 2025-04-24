source("Functions_LBgamlss.R")
source("LB_gamlss.R")

#------------------------------------------------------------------------------------------
# Simulação Sem Regressor
#------------------------------------------------------------------------------------------
set.seed(10)
n <- 1000
R <- 100
mu_true <- 0.7
mu_result <- c()

for (i in 1:R) {
  y <- rLB(n, mu_true)
  fit <- gamlss(y ~ 1, family = LB(), trace = F)
  logit_link<-make.link("logit")
  mu_result[i] <- logit_link$linkinv(fit$mu.coefficients)
}

result1 <- matrix(c(mu_true, mean(mu_result)), 1, 2)
colnames(result1) <- c("mu verdadeiro", "mu estimado")
rownames(result1) <- c("summary")
print(round(result1, 2))

#------------------------------------------------------------------------------------------
# Simulação Com Regressor
#------------------------------------------------------------------------------------------
X <- runif(n)
logit_link <- make.link("logit")
log_link <- make.link("identity")
b1 <- .7
b2 <- .5
mu_true <- logit_link$linkinv(b1 + b2 * X)

mu_result <- matrix(NA, R, 2)

for (i in 1:R) {
  y <- rLB(n, mu_true)
  fit1 <- gamlss(y ~ X, family = LB(), trace = F)
  mu_result[i, ] <- fit1$mu.coefficients
}

true_values <- c(b1, b2)
mean_values <- colMeans(mu_result)
bias_values <- (true_values - mean_values) / true_values * 100
eqm_values <- apply(mu_result, 2, var) + (true_values - mean_values)^2

result2 <- cbind(true_values, mean_values, bias_values, eqm_values)
colnames(result2) <- c("true value", "mean", "bias (%)", "eqm")
rownames(result2) <- c("b1", "b2")
print(round(result2, 4))

