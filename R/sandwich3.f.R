# Sandwich estimator 3
# 
# This function is internal and not to be called by the user

sandwich3.f <- function(Y, g1, g2)
{
  Y1 <- na.omit(Y[g1])
  Y2 <- na.omit(Y[g2])
  return(sum(sapply(Y2, function(x) sum(1*(Y1 == x) ))))   
}
