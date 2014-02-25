#' Quick forest plot of significantly up- and downregulated features.
#' 
#' This function plots a forest plot, normally used in meta analysis, to visualize the odds ratios (OR) and confidence intervals resulting from a call to \code{\link{uWMW}}. It is the default plot function for the \code{\link{uwmwRes}} objects. The results are ordered according to significance, and only the significant results are plotted. This function calls \code{\link{forestplot}} directly. 
#' 
#' @param x an uwmwRes object.
#' @param y ignored for uwmwRes objects
#' @param ... arguments passed down to \code{\link{forestplot}}
#' 
#' @return NULL invisibly
#' 
#' @seealso \code{\link{forestplot}}
#' 
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#' plot(NBtest)
#' 
#' @include uwmwRes_Class.R uwmwEstimate_Class.R S3Generics.R
#' @importFrom graphics plot
#' @export
#' @rdname plot
#' @aliases plot,uwmwRes,ANY-method
setMethod("plot",signature=c("uwmwRes","ANY"),plot.uwmwRes)

