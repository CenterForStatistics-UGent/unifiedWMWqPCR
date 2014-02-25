# This file contains S3 generics for plot, coef and vcov
# This way they don't need to be generated automatically.
# Check uwmw_Accessors and plotMethoduwmwRes for 
# the S4 counterparts

# To avoid conflicts.

plot.uwmwRes <- function(x,y,...){
  x.sort <- sort(x,which="p")
  id <- which(x.sort[,"p.value"] < 0.05)
  forestplot(x.sort,order=id,...)
}   