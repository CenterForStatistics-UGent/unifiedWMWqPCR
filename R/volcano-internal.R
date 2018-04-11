#' Parameters used for the function volcanoplot
#'
#' The parameters described can all be used in the different \link{volcanoplot} methods.
#'
#' The arguments \code{highlight} and \code{names} are chosen in such a way that the function can mimick the behaviour of the \code{volcanoplot} function in the package \code{limma}.
#'
#' @param highlight an integer, indicating how many of the top features should be highlighted. Defaults to 0.
#' @param names an optional character vector giving the names of the genes to be highlighted. If not specified, the features are numbered from most significant to least significant. Default value is \code{NULL}.
#' @param xlab see \code{\link{par}}. Defaults to \code{"Fit"}
#' @param ylab see \code{\link{par}}. Defaults to \code{"Significance"}
#' @param pch see \code{\link{par}}. Defaults to 16
#' @param cex see \code{\link{par}}. Defaults to 0.35 as in the \code{limma} package.
#' @param transf.x a function used to transpose the values on the X axis. Defaults to NULL, unless otherwise specified in the specific methods.
#' @param transf.y a function used to transpose the values on the Y axis. Defaults to NULL, unless otherwise specified in the specific method.
#' @param add.ref a charater value indicating if reference lines for x, y or both axes should be drawn. It takes the values "none", "x", "y" or "both". Default value is "none". If not specified but ref.x or ref.y is, then the reference line is drawn for the X and/or the Y axis.
#' @param ref.x a numerical vector, indicating the value or values at which the reference line(s) should be drawn on the X axis. Defaults to \code{NULL}.
#' @param ref.y a numerical value, indicating the value at which the reference line should be drawn on the Y axis. Defaults to \code{NULL} unless otherwise specified in the specific methods.
#' @param col.x The color for the reference line on the X axis. The argument is passed to the col argument of \code{\link{lines}}, so check that help page for possible values. Defaults to \code{"darkgrey"}
#' @param col.y see col.x, but for the Y axis. If \code{col.x} is specified and \code{col.y} is not, then the value of \code{col.x} is used for \code{col.y}.
#' @param lwd.x The line width for the reference line on the X axis. The argument is passed to the lwd argument of \code{\link{lines}}, so check that help page for possible values. Defaults to 1.
#' @param lwd.y see \code{lwd.x}, but for the Y axis. If \code{lwd.y} is not specified and \code{lwd.x} is, that value is used for \code{lwd.y} as well.
#' @param lty.x The line type for the reference line on the X axis. The argument is passed to the lty argument of \code{\link{lines}}, so check that help page for possible values. Defaults to 1.
#' @param lty.y see \code{lty.x}, but for the Y axis. If \code{lty.y} is not specified and \code{lty.x} is, that value is used for \code{lty.y} as well.
#'
#' @note These parameters belong to the internal function, eventually called by the different S4 methods. This internal function is not exported, and should not be called directly.
#' @examples
#' # see the help page of volcanoplot
#'
#' @importFrom graphics text
#' @include volcanoplot.R
#' @name volcanopar
#' @rdname volcanopar
#' @aliases volcano-par volcanoplot-par volcanoplotpar
NULL
# The function
.volcano <- function(fit,sig,
         highlight= 0,
         names = NULL,
         xlab = "Fit",
         ylab = "Significance",
         pch = 16,
         cex = 0.35,
         col = "black",      # Needed due to that pesky partial matching idiocy messing up when col is passed as an argument. You get multiple matching with col.x and col.y.
         transf.x = NULL,
         transf.y = NULL,
         add.ref=c("none","x","y","both"),
         ref.x = NULL,
         ref.y = NULL,
         col.x = "darkgrey",
         col.y = col.x,
         lwd.x = 1,
         lwd.y = lwd.x,
         lty.x = 1,
         lty.y = lty.x,
         ...){
  # Check the arguments
  if(length(fit) != length(sig)) stop("x and y should have the same length.")

  # Transform if needed
  if(identical(FALSE,transf.x)) transf.x <- NULL
  if(identical(FALSE,transf.y)) transf.y <- NULL

  if(is.null(transf.x)) transf.x <- identity
  if(is.null(transf.y)) transf.y <- identity

  x <- transf.x(fit)
  y <- transf.y(sig)

  # Set the arguments to determine which references to draw.
  mref <- missing(add.ref)
  add.ref <- match.arg(add.ref)

  if(mref){
    add.x <- if(is.null(ref.x)) FALSE else TRUE
    add.y <- if(is.null(ref.y)) FALSE else TRUE
  } else {
    add.x <- if(add.ref %in% c("x","both")) TRUE else FALSE
    add.y <- if(add.ref %in% c("y","both")) TRUE else FALSE
  }

  if(add.x & is.null(ref.x))
    stop("If add.x is TRUE, ref.x should be specified")
  if(add.y & is.null(ref.y))
    stop("If add.y is TRUE, ref.y should be specified")

  # Construct panel.first
  # Ancient coding, but very useful to get the reference lines
  # plotted underneath the points instead of on top of them.

  # make sure you don't mess up somebody elses use of panel.first
  mc <- match.call()
  panel.first <- mc$panel.first

  if(add.x || add.y){
    if(add.x)
      panel.first <- c(panel.first,
                       expression(
                         abline(v=ref.x,col=col.x,lwd=lwd.x,lty=lty.x)
                         )
      )
    if(add.y)
      panel.first <- c(panel.first,
                       expression(
                         abline(h=ref.y,col=col.y,lwd=lwd.y,lty=lty.y)
                       )
      )
  }

  plot(x,y,xlab = xlab,ylab = ylab,pch = pch,cex = cex,
       panel.first=eval(panel.first), ...)

  if (highlight > 0) {
    if (is.null(names))
      names <- 1:length(x)
    names <- as.character(names)
    o <- order(y, decreasing = TRUE)
    i <- o[1:highlight]
    text(x[i], y[i], labels = substring(names[i], 1, 8),
         cex = 0.8, col = "blue")
  }
  invisible()
}


