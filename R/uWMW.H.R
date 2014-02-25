# Calculate the uWMW test with housekeeping genes
# 
# This is an internal function, not to be called by the user.
# input is the same as in the function \code{\link{.uWMW}}.

uWMW.H <- function(data, groups, housekeeping.names)
{
  d.m <- data
  g1 <- groups == unique(groups)[1]
  g2 <- !g1
  names.hkp <- housekeeping.names
  
  # For each feature estimate the PI and obtain sample sizes
  f.info <- t(apply(d.m, 1, wmw.f, g1, g2))
  # To avoid unstable computations, replace estimated PI of 0 by 0.01 and those of 1 by 0.99
  f.info[,1] <- ifelse(f.info[,1] > 0.999, 0.999, f.info[,1])
  f.info[,1] <- ifelse(f.info[,1] < 0.001, 0.001, f.info[,1])
  
  # In the case of multiple housekeeping features, the PI's should be averaged out. 
  boolean.hkp <- rownames(f.info) %in% names.hkp
  PI.hkp.vec <- matrix(f.info[boolean.hkp,], ncol = ncol(f.info))
  PI.hkp <- sum(PI.hkp.vec[,1]*PI.hkp.vec[,2]*PI.hkp.vec[,3])/sum(PI.hkp.vec[,2]*PI.hkp.vec[,3]) 
  
  fv <- f.info[,1]
  fv[boolean.hkp] <- PI.hkp
  
  # Convert the PI's to odds
  # for the housekeeping PI
  odds.hkp <- PI.hkp/(1 - PI.hkp)
  # for the non-housekeeping features
  f.info.nohkp <- f.info[!(rownames(f.info) %in% names.hkp),1]
  odds.f.vec <- f.info.nohkp/(1-f.info.nohkp)
  
  beta.est <- c(intercept = log(odds.hkp), log(odds.f.vec) - log(odds.hkp))
  
  # Calculate sandwich estimator in an effecient manner
  
  Outer <- outer.h.f(beta.length = length(beta.est), 
                       weights = fv*(1-fv)*f.info[,2]*f.info[,3],
                       boolean.hkp = boolean.hkp)
  
  
  Middle1 <- outer.h.f(beta.length = length(beta.est), 
                         weights = apply(d.m,1, sandwich1.f, g1, g2) - 
                           fv^2*f.info[,2]*f.info[,3]^2,
                         boolean.hkp = boolean.hkp)
  
  
  
  Middle2 <- outer.h.f(beta.length = length(beta.est), 
                         weights = apply(d.m,1, sandwich2.f, g1, g2) -
                           fv^2*f.info[,2]^2*f.info[,3],
                         boolean.hkp = boolean.hkp)
  Middle3 <- outer.h.f(beta.length = length(beta.est), 
                         weights = f.info[,2]*f.info[,3]*(fv - fv^2) - 
                           apply(d.m,1, sandwich3.f, g1, g2)/4,
                         boolean.hkp = boolean.hkp)
  Middle <- Middle1 + Middle2 - Middle3
  
  cov.beta <- solve(Outer)%*%Middle%*%solve(Outer)
  return(list(logOR = beta.est[-1], SE = sqrt(diag(cov.beta))[-1], 
              coef = beta.est, vcov = cov.beta))
}
