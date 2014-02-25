# Estimation function based on wilcox.test
# 
# This function is internal and not to be called by the user.

#schattersfunctie op basis van wilcox
wmw.f <- function(x, g1, g2) 
{
  x <- as.matrix(x)
  res <- wilcox.test(x[g2,], x[g1,], alternative = "two.sided", 
                     paired = FALSE, exact = FALSE, corrrect= FALSE, na.omit = TRUE)
  n1 <- length(na.omit(x[g1,]))
  n2 <- length(na.omit(x[g2,]))
  PI.est <- res[["statistic"]]/(n1*n2)
  return(c(PI.est, n1 = n1, n2 = n2))
}