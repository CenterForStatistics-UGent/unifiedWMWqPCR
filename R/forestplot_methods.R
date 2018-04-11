# TODO : indication of significance

#' Making a forest plot of the results of uWMW
#'
#' This function creates a forest plot indicating the (log) odds ratios, the (log) odds or the probabilities for the results of the unified Wilcoxon-Mann-Whitney test.
#'
#' The function has methods for uwmwRes and uwmwEstimate objects. When called for an uwmwRes object, the requested estimate is first calculated using \code{\link{getEstimate}} and the result is passed on to the next method.
#'
#' Note that in either case, it is not possible to use the function on a subset of either type of object. The subsetting functions for \code{uwmwRes} and \code{uwmwEstimate} objects return matrices, and hence necessary information on the reference value is lost. To plot a subset of your data, use the \code{order} argument as shown in the examples.
#'
#' Adding a reference value to the plot only makes sense when plotting the log(odds), odds or probabilities. If log(OR) or OR are plotted, \code{addfit} is set to FALSE.
#'
#' The default settings plot a reference line at a location depending on the plotted estimate. For log(OR), the line is plotted at \code{refline = 0}. For OR, the line is plotted at \code{refline = 1}.
#'
#' @param x An object of class uwmwRes or uwmwEstimate.
#' @param ... parameters passed down to the internal functions. These can be any of the following.
#' @param estimate An optional character string defining which measure should be plotted. It can take the values \code{logor}, \code{or}, \code{logodds}, \code{odds} or \code{p} (for displaying the probability of differential expression, not the p-value!). Defaults to logor. Note that this argument is ignored when plotting an uwmwEstimate object.
#' @param annotate A logical value indicating whether the plot needs to be annotated, i.e. whether the values for the chosen measure and confidence interval should be displayed on the right of the plot. Defaults to TRUE
#' @param addfit A logical value indicating whether the reference measure should be plotted. See Details.
#' @param xlim The horizontal limits of the plot region. If unspecified, the function tries to set the horizontal plot limits to some sensible values. Should not be used by the user.
#' @param alim the actual x axis limits. If unspecified, an educated guess is taken by the function.
#' @param ylim The vertical limits of the plot. If unspecified, the function does what it thinks is best. Should not be used by the user.
#' @param at Position of the x axis tick marks and corresponding labels are placed. If unspecified, the function tries to position them at sensible values.
#' @param steps An integer indicating the number of tick marks on the X axis. Ignored when \code{at} is specified. Defaults to 5
#' @param level Numerical value between 0 and 1 to specify the width of the confidence interval. Defaults to 0.95 (95\% confidence interval).
#' @param digits integer specyfing the number of decimal places for tick mark labels and annotations. Can also be a vector of two integers. In that case, the first value specifies the number of decimal places for the annotations, the second for the x axis labels.
#' @param refline numerical value indicating where a reference line should be drawn. An NA value will prevent the line from being drawn. See Details.
#' @param xlab title for the x axis. If unspecified, the function tries to figure out the fitting title.
#' @param slab optional vector with names for the displayed genes.
#' @param mlab optional character string giving a label to the intercept estimate. If unspecified, this is created in the function if necessary.
#' @param ilab optional vector or matrix with character strings providing additional information that can be plotted next to the genes.
#' @param ilab.xpos Vector of numerical values specifying the x axis positions of the character vectors given via ilab. This has to be specified when ilab is specified.
#' @param ilab.pos integer from 1 to 4 specifying the alignment of the character vector(s) given via ilab (2 is right aligned, 4 is left aligned). Default is to center the labels.
#' @param order optional character string, character vector or numerical vector specifying how the genes should be ordered.
#' @param transf optional argument specifying the name of a function that should be used to transform the observed effect sizes, summary estimates, fitted values and confidence interval bounds (e.g., \code{transf=exp}). Defaults to \code{FALSE}, which means that no transformation is used.
#' @param atransf optional argument specifying the name of a function that should be used to transform the x-axis labels and annotations (e.g., transf=exp). Defaults to \code{FALSE}, which means that no transformation is used.
#' @param targs optional arguments needed by the function specified via transf or atransf.
#' @param rows optional vector specifying the horizontal position for the plotted results. If unspecified, the layout happens automatically. See Details and Examples for more information.
#' @param efac numerical value specifying the vertical expansion of the arrows, summary estimate symbols and Ci limits. Normally the default of 1 should work just fine.
#' @param pch plotting symbol used for the observed effect sizes. By default, it's a filled square.
#' @param psize optional vector with the point sizes for the observed effects. If set to \code{NULL}, the point sizes are drawn proportional to the value of the log of the test statistic. Defaults to 1.
#' @param cex optional numerical value for expansion of text and symbols. See also \code{\link{par}}.
#' @param cex.lab Optional numerical value for expansion of the axis title.
#' @param cex.axis Optional numerical value for expansion of the x axis labels.
#' @param col character string specifying the color used for the individual estimates.
#' @param border character string specifying the color used for the border of the individual estimates.
#' @param refcol Character string specyfying the color of the reference line. Defaults to red.
#' @param predcol character string specifying the color of the estimated reference value. Ignored if estimate is \code{logor} or \code{or}.
#'
#' @return NULL invisibly
#'
#' @author This code is adapted by Joris Meys from the function \code{forest.rma} (\code{metafor} package). The original function is written by W. Viechtbauer.
#'
#' @note Thanks to the work of W.Viechtbauer, \code{forestplot} provides many possibilities for tweaking and customizing the plots. Many of the arguments work the same as in the function \code{forest.rma} (\code{metafor} package). You can always check the help file of \code{forest.rma} for more illustrations on the different arguments.
#'
#' @note This function is currently implemented using an internal function that expects an uwmwEstimate object. In a future version, the internal function will be rewritten to be more generic. This will enable the definition of methods for other classes without need to change the internal function itself.
#'
#' @section Warning:
#'  Although the internal function is shown here (merely for illustration of the arguments and defaults), the user shouldn't be calling this one directly. The function is not exported.
#'
#' @examples
#'
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#' sigid <- which(pval(NBtest) < 0.05)
#' forestplot(NBtest,"logodds",order=sigid)
#'
#' nameid <- c("hsa-mir-30a-3p","hsa-mir-30a-5p")
#' forestplot(NBtest,"p",order=nameid,addfit=FALSE,
#'            refline=NA,main="Comparison 30a")
#'
#' forestplot(NBtest,"p",order=nameid,addfit=FALSE,
#'            refline=0.5,main="Comparison 30a",
#'            alim=c(0,1),xlim=c(-1,2),at=c(0,0.5,1))
#'
#' @importFrom graphics par abline strheight segments polygon axis mtext points
#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @name forestplot
#' @rdname forestplot-methods
#' @aliases forestplot,ANY-method
setMethod('forestplot',
          signature="ANY",
          function(x,...){
            stop("Forestplot currently only works on uwmwRes and uwmwEstimate objects. See ?uWMW and ?uwmwRes.")
          })

#' @rdname forestplot-methods
#' @aliases forestplot,uwmwRes-method
setMethod('forestplot',
          signature="uwmwRes",
          function(x,estimate=c("logor","logodds","or","odds","p"),
                   level=0.95,...){

            estimate <- match.arg(estimate)

            est <- getEstimate(x,estimate,ci=level)

            forestplot(est,level=level,...)
          })

#' @rdname forestplot-methods
#' @aliases forestplot,uwmwEstimate-method
setMethod('forestplot',
          signature="uwmwEstimate",
          function(x,...){
             forestplot.internal(x,...)
          })

#' @rdname forestplot-methods
forestplot.internal <- function (x,
                        annotate = TRUE,
                        addfit = TRUE,
                        xlim = NULL,      # Horizontal limit
                        alim = NULL,
                        ylim = NULL,
                        at = NULL,
                        steps = 5,
                        level = 0.95,
                        digits = 2,
                        refline = NULL,
                        xlab = NULL,
                        slab = NULL,
                        mlab = NULL,
                        ilab = NULL,
                        ilab.xpos = NULL,
                        ilab.pos = NULL,
                        order = NULL,
                        transf = FALSE,
                        atransf = FALSE,
                        targs = NULL,
                        rows = NULL,
                        efac = 1,
                        pch = 15,
                        psize = 1,
                        col = "darkgrey",
                        border = "darkgrey",
                        cex = NULL,
                        cex.lab = NULL,
                        cex.axis = NULL,
                        refcol="red",
                        predcol=refcol, ...)
{



  # Control arguments
  estimate <- esttype(x)

  if(estimate !='logor'){
    if(transf | atransf) warning("transf/atransf arguments are ignored when measure is not 'logor'.")
    transf <- atransf <- FALSE
  }

  na.act <- getOption("na.action")
  if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail",
                            "na.pass")))
    stop("Unknown 'na.action' specified under options().")

  transf.char <- deparse(substitute(transf))
  atransf.char <- deparse(substitute(atransf))
  if (transf.char != "FALSE" && atransf.char != "FALSE")
    stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")

  names.x <- names(x)
  # Check order and reorder observations if necessary.
  estmat <- as.matrix(x)
  if(!is.null(order)){
    if(is.character(order)) order <- match(order,names.x,0L)
    yi <- estmat[order,"est"]         # estimate
    ci.lb <- estmat[order,"ll"]       # ll ci
    ci.ub <- estmat[order,"ul"]       # ul ci
    names.x <- names.x[order]
    vi <- (se(x)^2)[order]     # variance of estimates
  } else {
    yi <- estmat[,"est"]                # estimate
    ci.lb <- estmat[,"ll"]              # ll ci
    ci.ub <- estmat[,"ul"]              # ul ci
    vi <- se(x)^2              # variance of estimates
  }

  k <- length(yi)                     # number of estimates

  vi[which(vi<=0)] <- NA   # Make sure all variances are valid

  # Reference value

  if(!estimate %in% c("logor","or") ){

    preds <- ref(x)
    b  <- preds['refest']
    b.ci.lb <- preds['ll']
    b.ci.ub <- preds['ul']

  } else {
    addfit <- FALSE
  }

  if(is.null(refline))
    refline <- switch(estimate,
                      "logor"= 0,
                      "or" = 1,
                      "logodds" =,
                      "odds" = ,
                      "p" = b)

  ## CONSTRUCTION OF THE PLOT

  # Determining the rows

  if (is.null(rows)) {
    rows <- seq_len(k)
  }
  else {
    if (length(rows) == 1L) {
      rows <- seq_len(rows)
    }
  }
  if (length(rows) != length(yi))
    stop("Number of outcomes does not correspond to the length of the rows argument.")

  # Get the labels if necessary
  if (is.null(slab)) {
    slab <- names.x

  }

  if (length(yi) != length(slab))
    stop("Number of outcomes does not correspond to the length of the slab argument.")

  if (is.vector(ilab))
    ilab <- cbind(ilab)

  # Take care of the pch
  if (length(pch) == 1L)
    pch <- rep(pch, k)

  if (length(pch) != length(yi))
    stop("Number of outcomes does not correspond to the length of the pch argument.")
  options(na.action = "na.exclude")

  # Check psize
  if (!is.null(psize)) {
    if (length(psize) == 1L)
      psize <- rep(psize, k)
    if (length(psize) != length(yi))
      stop("Number of outcomes does not correspond to the length of the psize argument.")
  }

  if (is.null(psize)) {
    if (any(is.na(vi))) {
      psize <- rep(1, k)
    }
    else {
      psize <- vi/sum(vi, na.rm = TRUE)
      psize <- (psize - min(psize, na.rm = TRUE)) /
        (max(psize,na.rm = TRUE) - min(psize, na.rm = TRUE))
      psize <- (psize * 1) + 0.5
    }
  }


  ##

  rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
  if (annotate) {
    plot.multp.l <- 1.2
    plot.multp.r <- 1.2
    axis.multp.l <- 0.2
    axis.multp.r <- 0.2
  }
  else {
    plot.multp.l <- 1.2
    plot.multp.r <- 0.4
    axis.multp.l <- 0.2
    axis.multp.r <- 0.2
  }

  if (is.null(xlim)) {
    xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l,
              max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
    xlim <- round(xlim, digits)
  }

  if (is.null(alim)) {
    if (is.null(at)) {
      alim <- c(min(ci.lb, na.rm = TRUE) - rng * axis.multp.l,
                max(ci.ub, na.rm = TRUE) + rng * axis.multp.r)
      alim <- round(alim, digits)
    }
    else {
      alim <- range(at)
    }
  }
  alim <- sort(alim)
  xlim <- sort(xlim)
  if (xlim[1] > min(yi, na.rm = TRUE)) {
    xlim[1] <- min(yi, na.rm = TRUE)
  }
  if (xlim[2] < max(yi, na.rm = TRUE)) {
    xlim[2] <- max(yi, na.rm = TRUE)
  }
  if (alim[1] > min(yi, na.rm = TRUE)) {
    alim[1] <- min(yi, na.rm = TRUE)
  }
  if (alim[2] < max(yi, na.rm = TRUE)) {
    alim[2] <- max(yi, na.rm = TRUE)
  }
  if (alim[1] < xlim[1]) {
    xlim[1] <- alim[1]
  }
  if (alim[2] > xlim[2]) {
    xlim[2] <- alim[2]
  }
  if (is.null(ylim)) {
    if (addfit) {
      ylim <- c(-1.5, k + 3)
    }
    else {
      ylim <- c(0.5, k + 3)
    }
  }
  else {
    ylim <- sort(ylim)
  }
  if (is.null(at)) {
    at <- seq(alim[1], alim[2], length = steps)
  }
  else {
    at[at < alim[1]] <- alim[1]
    at[at > alim[2]] <- alim[2]
    at <- unique(at)
  }

  at.lab <- at

  at.lab <- round(at.lab, digits)

  par.mar <- par("mar")
  par.mar.adj <- par.mar - c(0, 3, 1, 1)
  par.mar.adj[par.mar.adj < 1] <- 1
  par(mar = par.mar.adj)
  on.exit(par(mar = par.mar))
  plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
       yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
  abline(h = ylim[2] - 2, ...)
  par.usr <- par("usr")
  height <- par.usr[4] - par.usr[3]
  lheight <- strheight("O")
  cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 *
                                                          k * lheight), 1)
  if (is.null(cex)) {
    cex <- par("cex") * cex.adj
  }
  else {
    if (is.null(cex.lab))
      cex.lab <- cex
    if (is.null(cex.axis))
      cex.axis <- cex
  }
  if (is.null(cex.lab))
    cex.lab <- par("cex") * cex.adj
  if (is.null(cex.axis))
    cex.axis <- par("cex") * cex.adj

  if (is.numeric(refline))
    segments(refline, ylim[1] - 5, refline, ylim[2] - 2,
             lty = "dotted", col=refcol, ...)

  if (addfit) {

    polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 +
                                                   (height/100) * cex * efac, -1, -1 - (height/100) *
                                                   cex * efac), col = predcol, ...)
    if (is.null(mlab))
      mlab <- paste("Reference",estimate)
    text(xlim[1], -1, mlab, pos = 4, cex = cex, ...)
  }


  axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis,
       ...)

  if(is.null(xlab)) xlab <- paste("Observed Outcome:",estimate)

  if (!is.null(xlab))
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2,
          line = 2.5, cex = cex.lab, ...)
  for (i in 1:k) {
    if (is.na(yi[i]) )
      next
    segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i],
                                                  alim[2]), rows[i], ...)
    if (ci.lb[i] >= alim[1]) {
      segments(ci.lb[i], rows[i] - (height/150) * cex *
                 efac, ci.lb[i], rows[i] + (height/150) * cex *
                 efac, ...)
    }
    else {
      polygon(x = c(alim[1], alim[1] + (1.4/100) * cex *
                      (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex *
                      (xlim[2] - xlim[1]), alim[1]), y = c(rows[i],
                                                           rows[i] + (height/150) * cex * efac, rows[i] -
                                                             (height/150) * cex * efac, rows[i]), col = "black",
              ...)
    }
    if (ci.ub[i] <= alim[2]) {
      segments(ci.ub[i], rows[i] - (height/150) * cex *
                 efac, ci.ub[i], rows[i] + (height/150) * cex *
                 efac, ...)
    }
    else {
      polygon(x = c(alim[2], alim[2] - (1.4/100) * cex *
                      (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex *
                      (xlim[2] - xlim[1]), alim[2]), y = c(rows[i],
                                                           rows[i] + (height/150) * cex * efac, rows[i] -
                                                             (height/150) * cex * efac, rows[i]), col = "black",
              ...)
    }
  }
  text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
  if (!is.null(ilab)) {
    for (l in 1:NCOL(ilab)) {
      text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l],
           cex = cex, ...)
    }
  }

  if (annotate) {


      if (addfit) {
        annotext <- round(cbind(c(yi, b), c(ci.lb, b.ci.lb),
                                c(ci.ub, b.ci.ub)), digits)
      }
      else {
        annotext <- round(cbind(yi, ci.lb, ci.ub), digits)
      }



    annotext <- matrix(apply(annotext, 2, format, nsmall = digits),
                         ncol = 3)
    annotext <- cbind(annotext[, 1], " [ ", annotext[,
                                                       2], " , ", annotext[, 3], " ]")

    annotext <- apply(annotext, 1, paste, collapse = "")
    if (addfit) {
      text(x = xlim[2], c(rows, -1), labels = annotext,
           pos = 2, cex = cex, ...)
    }
    else {
      text(x = xlim[2], rows, labels = annotext, pos = 2,
           cex = cex, ...)
    }
  }
  points(yi, rows, pch = pch, cex = cex * psize, ...)
  if ( addfit)
    abline(h = 0, ...)
  invisible()


}

