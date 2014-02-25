# Internal calculation function outer.ome.f
# 
# Function is internal and not to be called by the user.

outer.o.f <- function(beta.length, weights)
{
  le <- length(weights)
  Outer1 <- matrix(nrow = beta.length, ncol = beta.length, 0)
  
  Outer1[1,1] <- sum(weights[-le])
  Outer1[-1,1] <- Outer1[1, -1] <- diag(Outer1)[-1] <- weights[-le]
  
  Outer2 <- matrix(nrow = beta.length, ncol = beta.length, 1)
  Outer2[-1,1] <- Outer2[1, -1] <- -1
  Outer2 <- weights[le]*Outer2
  
  Outer <- Outer1 + Outer2
  return(Outer)
}

