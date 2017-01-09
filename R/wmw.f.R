# Estimation function based on wilcox.test
# 
# This function is internal and not to be called by the user.

#schattersfunctie op basis van wilcox
wmw.f <- function(x, g1, g2) 
{
  # deals with missing values as well (see wilcox.test.default)
  x <- x[is.finite(x)]
  
  #Calculation of the W statistics (with g2 as group 1 !!!)
  x1 <- x[g2]
  x2 <- x[g1]
  r <- rank(c(x1,x2))
  n1 <- as.double(length(x1))
  n2 <- as.double(length(x2))
  W <- sum(r[seq_along(x1)]) - n1 *(n1+1) / 2
  
  PI.est <- W/(n1*n2)
  return(c(PI.est, n1 = n1, n2 = n2))
}
