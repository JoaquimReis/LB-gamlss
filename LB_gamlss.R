library(gamlss)
source("Functions_LBgamlss.R")

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
