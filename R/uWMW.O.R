# Calculation of uWMW without housekeeping genes.
# 
# This function is internal and should not be called by the user. The input is the same as in \code{\link{.uWMW}}

uWMW.O <- function(data, groups)
{
  d.m <- data
  g1 <- groups == unique(groups)[1]
  g2 <- !g1
  
  f.info <- t(apply(d.m, 1, wmw.f, g1, g2))
  f.info[,1] <- ifelse(f.info[,1] > 0.999, 0.999, f.info[,1])
  f.info[,1] <- ifelse(f.info[,1] < 0.001, 0.001, f.info[,1])
  fv <- f.info[,1]
  
  odds.f.vec <- f.info[,1]/(1-f.info[,1])
  log.odds.int <- mean(log(odds.f.vec))
  le <- nrow(f.info)
  beta.est <- c(intercept = log.odds.int, log(odds.f.vec[-le]) - log.odds.int)
  
  Outer <- outer.o.f(beta.length = length(beta.est), 
                       weights = fv*(1-fv)*f.info[,2]*f.info[,3])
  
  Middle1 <- outer.o.f(beta.length = length(beta.est), 
                         weights = apply(d.m,1, sandwich1.f, g1, g2) - 
                           fv^2*f.info[,2]*f.info[,3]^2)
  Middle2 <- outer.o.f(beta.length = length(beta.est), 
                         weights = apply(d.m,1, sandwich2.f, g1, g2) - 
                           fv^2*f.info[,2]^2*f.info[,3])   
  Middle3 <- outer.o.f(beta.length = length(beta.est), 
                         weights = f.info[,2]*f.info[,3]*(fv - fv^2) - 
                           apply(d.m,1, sandwich3.f, g1, g2)/4)
  Middle <- Middle1 + Middle2 - Middle3
  
  cov.beta <- solve(Outer)%*%Middle%*%solve(Outer)
  logOR <- c(beta.est[-1], -sum(beta.est[-1]))
  names(logOR) <- rownames(d.m)
  SE <- c(sqrt(diag(cov.beta))[-1], sqrt(sum(cov.beta[-1,-1])))             
  return(list(logOR = logOR, SE = SE, coef = beta.est, vcov = cov.beta))
}
