#' S3 method as.matrix for uwmwRes objects.
#' 
#' For \code{\link{uwmwRes}} objects, an \code{as.matrix} method is defined that transforms the object to a numeric matrix with the following columns:
#' \describe{
#'   \item{\code{logor}:}{ The log odds ratio values}
#'   \item{\code{se}:}{ The standard errors on the log OR values.}
#'   \item{\code{or}}{The odds ratio values}
#'   \item{\code{z.value}}{the z values related to the log OR values}
#'   \item{\code{p.value}}{The p values related to the log OR values}
#' }
#' The row names are the names of the tested genes/features. The matrix takes the ordering in the object into account if necessary.
#' 
#' @param x a uwmwRes object
#' @param ... currently ignored
#' 
#' @return A matrix containing the columns specified above.
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#' as.matrix(NBtest)
#' 
#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @name as.matrix.uwmwRes
#' @rdname as.matrix.uwmwRes
#' @method as.matrix uwmwRes
#' @export
#' 

as.matrix.uwmwRes <- function(x,...){
  out <- cbind(
    logor        = logor(x),
    se           = se(x),
    or           = oddsRatio(x),
    z.value      = zval(x),
    p.value      = pval(x)
  )
  rownames(out) <- names(x)
  out
}

# This adds the S4 requested by BioConductor. 
## S4 method calls S3
#' @name as.matrix.uwmwRes
#' @rdname as.matrix.uwmwRes
#' @aliases as.matrix,uwmwRes-method
setMethod("as.matrix", "uwmwRes", as.matrix.uwmwRes)

## S4 coercion calls S4 method
setAs("uwmwRes", "matrix", function(from) as.matrix(from))
