# Internal function outer.hme.f
# 
# This function is internal and not to be called by the user

outer.h.f <- function(beta.length, weights, boolean.hkp)
{
  Outer <- matrix(nrow = beta.length, ncol = beta.length, 0)
  
  Outer[1,1] <- sum(weights)
  Outer[-1,1] <- Outer[1, -1] <- diag(Outer)[-1] <-  weights[!boolean.hkp]
  return(Outer)
}