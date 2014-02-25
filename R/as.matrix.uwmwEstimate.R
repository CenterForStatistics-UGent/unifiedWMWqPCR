#' S3 method as.matrix for uwmwEstimate objects.
#' 
#' For \code{\link{uwmwEstimate}} objects, an \code{as.matrix} method is defined that transforms the object to a numeric matrix with the following columns:
#' \describe{
#'   \item{\code{est}:}{ The estimates}
#'   \item{\code{se}:}{ The standard errors on estimates.}
#'   \item{\code{ll}}{lower limit of the confidence interval}
#'   \item{\code{ul}}{upper limit of the confidence interval}
#' }
#' The row names are the names of the tested genes/features. The matrix takes the ordering in the object into account if necessary.
#' 
#' @param x a uwmwEstimate object
#' @param ... currently ignored
#' 
#' @return A matrix with the columns specified above.
#' 
#' @author Joris Meys
#' 
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#' NBest <- getEstimate(NBtest,"p")
#' as.matrix(NBtest)
#' 
#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @name as.matrix.uwmwEstimate
#' @rdname as.matrix.uwmwEstimate
#' @method as.matrix uwmwEstimate
#' @export

as.matrix.uwmwEstimate <- function(x,...){
  out <- cbind(
    est     = x@est,       
    se      = x@se,
    ul      = x@ul,
    ll      = x@ll
  )
  rownames(out) <- names(x)

#  class(out) <- c("uwmwMatrix",class(out))
  out
}
# NOTE : there's no separate getters for est, se, ...
# as.matrix is used to extract that information. 
# This to avoid an overload of different accessor functions.

# This adds the S4 requested by BioConductor. 
## S4 method calls S3
#' @name as.matrix.uwmwEstimate
#' @rdname as.matrix.uwmwEstimate
#' @aliases as.matrix,uwmwEstimate-method
setMethod("as.matrix", "uwmwEstimate", as.matrix.uwmwEstimate)

## S4 coercion calls S4 method
setAs("uwmwEstimate", "matrix", function(from) as.matrix(from))

