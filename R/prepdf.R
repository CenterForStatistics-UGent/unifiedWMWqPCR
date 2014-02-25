
# Internal Function to prepare a data frame as a matrix
prep.df <- function(x,feat.names,subjects,value){

  su <- unique(x[,subjects])
  fu <- unique(x[,feat.names])
  
  ns <- length(su)
  nf <- length(fu)
  if(nrow(x) != ns*nf)
    stop("Incomplete data currently not supported.")
  
  id <- order(x[,subjects],x[,feat.names])
  
  dmat <- matrix(x[id,value],ncol=ns,nrow=nf)
  rownames(dmat) <- sort(fu)
  colnames(dmat) <- sort(su)
  dmat
}
