#' Make a volcano plot of the outcome of a uWMW test
#' 
#' This function creates a volcano plot of the outcome of a uWMW test. In this plot, the p-value is plotted against the odd ratio or the log-odds ratio, depending on what interests you most. It allows to quickly see what proportion of results is significant, and what proportion of significant results has a biologically significant chance of being upregulated.
#' 
#' The methods described here all use the same internal function to create the plot
#' 
#' @param fit a uwmwRes object, a matrix with 2 columns or a numeric vector with either odds ratios or the log of odds ratios.
#' @param esttype The estimate that should be plotted from the uwmwRes object. See also \code{\link{getEstimate}}.
#' @param transf.y a function used to transpose the values on the Y axis. For uwmwRes objects, defaults to -log10(y). Passed down to internal function.
#' @param ref.y The location on the Y axis for the reference line.
#' @param coef for compatibility with limma package, see \code{\link{volcanopar}}
#' @param highlight for compatibiity with limma package, see \code{\link{volcanopar}}
#' @param ... arguments passed to the internal function. See \code{\link{volcanopar}}
#' @return invisible NULL 
#' 
#' @examples
#' data(NBmat)
#' housekeeping.id <- grep("let",rownames(NBmat),value=TRUE)
#' NB.Htest <- uWMW(NBmat,NBgroups,housekeeping.id)
#' 
#' volcanoplot(NB.Htest)
#' volcanoplot(NB.Htest,"p")
#' 
#' @author Joris Meys
#' @docType methods
#' @name volcanoplot
#' @rdname volcanoplot-methods
#' @aliases volcanoplot,uwmwRes-method
setMethod("volcanoplot",
  signature=c("uwmwRes"),
  function(fit, 
           esttype=c("logor","or","logodds","odds","p"), 
           transf.y=function(i) -log10(i),ref.y = -log10(0.05),...){
    
    esttype <- match.arg(esttype)
    mc <- match.call()
    # set plotting parameters
    if(is.null(mc$xlab)) xlab <- esttype
    if(is.null(mc$ylab)) ylab <- expression(- log[10]~p ) 
    
    
    yval <- pval(fit)        
    xval <- getEstimate(fit,esttype,se.fit=FALSE,ci=FALSE)
    
    volcanoplot(xval,yval,transf.y=transf.y,ref.y=ref.y,
                xlab=xlab,ylab=ylab,...)
  })


# USE A PLOT METHOD FOR MATRIX
#' @rdname volcanoplot-methods
#' @aliases volcanoplot,matrix-method
setMethod("volcanoplot",
  signature=c("matrix"),
  function(fit,...){
    if(ncol(fit) > 2) stop("matrix should have two columns")
    
    volcanoplot(fit[,1],fit[,2],...)
  }
)

# General method
#' @rdname volcanoplot-methods
#' @aliases volcanoplot,numeric-method
setMethod("volcanoplot",
  signature=c("numeric"),
  function(fit,...){
    .volcano(fit,...)
  }
)

# MArrayLM method
# This is re-implemented to use the internal function of 
# the uWMWqPCR package so it always reacts the same on 
# the different arguments. 
# Passing to limma::volcanoplot was another option, but
# this seemed better in the long run. 
#' @rdname volcanoplot-methods
#' @aliases volcanoplot,MArrayLM-method
setMethod("volcanoplot",
  signature=c("MArrayLM"),
  function(fit,coef=1,highlight=0,...){
    if (is.null(fit$lods)) 
      stop("No B-statistics found, perhaps eBayes() not yet run")
    x <- as.matrix(fit$coef)[, coef]
    y <- as.matrix(fit$lods)[, coef]
    
    # construct the call
    mc <- match.call()
    mc$coef <- NULL
    if(is.null(mc$names)) mc$names <- fit$genes$ID
    if(is.null(mc$xlab)) mc$xlab <- "Log Fold Change"
    if(is.null(mc$ylab)) mc$ylab <- "Log Odds"
    if(is.null(mc$pch)) mc$pch <- 16
    if(is.null(mc$cex)) mc$cex <- 0.35
        
    mc[[1]] <- volcanoplot
    mc$fit <- x
    mc$sig <- y
    eval(mc)
  }        )
