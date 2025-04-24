library(gamlss)

LB <- expression(
  
  log(12 / ((mu / (mu + 24))^(-1/2) - 5)) + 
    (4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * log(y) + 
    log(1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
)


mLB <- D(LB,"mu")

LB<-function (mu.link = "logit")
{
  mstats <- checklink("mu.link", "LB", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  structure(list(family = c("LB", "Log-Bilal"),
                 parameters = list(mu = TRUE),
                 nopar = 1,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 
                 mu.linkfun = mstats$linkfun,
                 
                 mu.linkinv = mstats$linkinv,
                 
                 mu.dr = mstats$mu.eta,
                 
                 dldm = function(y, mu) {
                   dldm <- eval(mLB)
                   dldm
                 },
                 d2ldm2 = function(y,mu) {
                   dldm <- eval(mLB)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 
                 G.dev.incr = function(y, mu, w, ...) -2 * log(dLB(y=y, mu=mu)),
                 rqres = expression(
                   rqres(pfun = "pLB", type = "Continuous", y = y, mu = mu)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}

#------------------------------------------------------------------------------------------
# Função de densidade Log-Bilal
#------------------------------------------------------------------------------------------
dLB <- function(y, mu = 0.75, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  
  fy1 <- (12 / ((mu / (mu + 24))^(-1/2) - 5)) * y^(4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * (1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
  
  if (log == FALSE) fy <- fy1 else fy <- log(fy1)
  fy
}

#integrate(dLB, 0, 1) # checando se a densidade integra para 1

#------------------------------------------------------------------------------------------
# Função de distribuição acumulada Log-Bilal
#------------------------------------------------------------------------------------------
pLB <- function(q, mu = 0.75, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  
  cdf1<- 3 * q^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * q^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  
  if (lower.tail == TRUE) cdf <- cdf1 else cdf <- 1 - cdf1
  if (log.p == FALSE) cdf <- cdf else cdf <- log(cdf)
  
  cdf
}

# pLB(0.5)
# integrate(dLB, 0, 0.5) # checando com a densidade

#------------------------------------------------------------------------------------------
# Função inversa (quantílica) Log-Bilal via método da inversão numérica
#------------------------------------------------------------------------------------------
inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol = 1e-10)$root
  }
}

qLB <- function(u, mu){
  if (any(u <= 0 | u >= 1)) stop(paste("u must be between 0 and 1", "\n", ""))
  
  f<- function(y) {3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  }
  q<-inverse(f, lower = 0, upper = 1)
  return(q(u))
}

# y=.7
# mu=.15
# u=pLB(y,mu)
# qLB(u,mu)

#------------------------------------------------------------------------------------------
# Função de geração de números aleatórios - Distribuição Log-Bilal
#------------------------------------------------------------------------------------------
rLB <- function(n, mu) {
  u <- runif(n)
  y <- mapply(function(ui, mui) qLB(ui, mu = mui), u, mu)
  return(y)
}

# Exemplo de uso:
# y <- rLB(1000, mu = 0.6)
# hist(y, breaks = 30, main = "Amostras da distribuição Log-Bilal")

#------------------------------------------------------------------------------------------
# Simulação Com Regressor
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

                 