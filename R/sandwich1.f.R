# Sandwich estimator 1
# 
# This function is internal and not to be called by the user

sandwich1.f <- function(Y, g1, g2)
{
  Y1 <- na.omit(Y[g1])
  Y2 <- na.omit(Y[g2])
  return(sum(sapply(Y1, function(x) sum(1*(x < Y2) + 0.5*(x == Y2)))^2))  
} 