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