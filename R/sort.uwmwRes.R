#' Sort and order method for uWMWRes objects
#' 
#' This functions provide sorting functionality for \code{\link{uwmwRes}} objects. It allows to sort the values in the object in order to get the genes with the highest OR, lowest p value, ... The function uses \code{\link{order}} underneath
#' 
#' The function does not change the internal order, but changes the slot \code{id} in the object. This slot is used by other functions to give the requested values 
#' 
#' @param x a uwmwRes object
#' @param decreasing a logical value indicating whether values should be sorted in increasing or decreasing order. 
#' @param which a character value indicating on which values should be used to sort on. The possible values are: "or" for sorting on the odds ratio, "p" for sorting on the p value, "se" for sorting on the standard error or "name" for sorting on the gene names.
#' @param na.last a logical value indicating whether NA results should be sorted at the end. See \code{\link{order}} for more information.
#' @param ... currently ignored
#' 
#' @return \code{sort} returns a sorted uwmwRes object.
#' 
#' @seealso The functions \code{\link{is.unsorted}}, \code{\link{orderedBy}} and \code{\link{getOrder}} to check on the ordered state of a \code{\link{uwmwRes}} object.
#' 
#' \code{order} returns the order of the uwmwRes object. 
#' 
#' \code{unorder} returns a uwmwRes object with the order removed.
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#' NBsort <- sort(NBtest,which="name")
#' NBsort[1:10]
#' 
#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @name sort.uwmwRes
#' @rdname sort.uwmwRes
#' @aliases unorder
#' @export
#' @docType methods
#' @export
setMethod("sort",signature="uwmwRes",function(x,decreasing=FALSE,which=c("or","p","se","name","z"),na.last=TRUE,...){
  which <- match.arg(which)
  
  which <- match.arg(which)
  getf <- switch(which,"or"=logor,"p"=pval,"se"=se,"name"=names,"z"=zval)
  x@id <- order(getf(x),decreasing=decreasing,na.last=na.last)
  x@orderedBy <- which
  
  return(x)
}  )

#' @rdname sort.uwmwRes
#' @export
setMethod("unorder",
          signature="ANY",
          function(x,...){
            theclass <- class(x)
            stop(paste("unorder is not defined for",theclass,"objects."))
          })

#' @rdname sort.uwmwRes
#' @export
setMethod("unorder",
          signature="uwmwRes",
          function(x,...){
            if(!is.unsorted(x)){
              x@orderedBy <- "none"
              x@id <- seq_len(length(x))
            }
            return(x)
          })
# Just to make sure that they don't try this one.
#' @rdname sort.uwmwRes
setMethod("order",
          signature="uwmwRes",
          function(...,na.last=TRUE,decreasing=FALSE
                   ){
            stop("order gives no sensible output when used on uwmwRes objects.")
            NULL
            
          })
