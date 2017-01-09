#Internal function for calculating the unified Wilcox-Mann-Whitney
# 
# This function is used internally in the S4 method uWMW, and should not be called by the user in general. 
# 
# @param data A matrix with the qPCR data, with the rownames as names for the genes and the columns as replicates.
# @param groups A vector as long as the number of columns in the data matrix. The vector indicates from which group every replicate comes.
# @param housekeeping.names An optional character vector with the names of the housekeeping genes used in the test. If set to NULL, the test is carried out without taking housekeeping genes into account.


.uWMW <- function(data, groups, housekeeping.names = NULL)
{
  if(sum(is.element(housekeeping.names, rownames(data))) != length(housekeeping.names)) 
    stop("One or more housekeeping features names are not valid, please check for typos")
  if(length(groups) != ncol(data))
    stop("Length of 'groups' is not equal to the number of columns of 'data' ")
  
  grouplev <- as.character(unique(groups))
  if(length(grouplev) != 2)
    stop(paste("uWMW needs exactly 2 groups, but",length(grouplev),"were specified."))
  
  if(is.null(housekeeping.names))  
  {
    t.type <- "O"
    housekeeping.names <- character(0)
    results <- uWMW.O(data, groups)
  } else
  {
    t.type <- "H"
    results <- uWMW.H(data, groups, housekeeping.names)
  }
  if(sum(diag(results[[4]]) < 0) > 0) warning("One or more variance estimates are negative")
  z.tmp <- results[[1]]/results[[2]]
  p.tmp <- 2*(1 - pnorm(abs(z.tmp)))

  rownames(results[[4]]) <- colnames(results[[4]]) <- names(results[[3]])

  uwmwRes(
    type         = t.type,
    housekeeping = housekeeping.names,
    names        = names(results[[1]]),
    logOR        = results[[1]],
    se           = results[[2]],
    OR           = exp(results[[1]]),
    z.value      = z.tmp,
    p.value      = p.tmp,
    coef         = results[[3]],
    vcov         = results[[4]],
    groupinfo    = grouplev
    )
}
