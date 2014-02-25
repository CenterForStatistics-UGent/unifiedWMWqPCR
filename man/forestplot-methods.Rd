\docType{methods}
\name{forestplot}
\alias{forestplot}
\alias{forestplot,ANY-method}
\alias{forestplot,uwmwEstimate-method}
\alias{forestplot,uwmwRes-method}
\alias{forestplot.internal}
\title{Making a forest plot of the results of uWMW}
\usage{
\S4method{forestplot}{ANY}(x, ...)

\S4method{forestplot}{uwmwRes}(x, estimate = c("logor", "logodds", "or",
  "odds", "p"), level = 0.95, ...)

\S4method{forestplot}{uwmwEstimate}(x, ...)

forestplot.internal(x, annotate = TRUE, addfit = TRUE, xlim = NULL,
  alim = NULL, ylim = NULL, at = NULL, steps = 5, level = 0.95,
  digits = 2, refline = NULL, xlab = NULL, slab = NULL, mlab = NULL,
  ilab = NULL, ilab.xpos = NULL, ilab.pos = NULL, order = NULL,
  transf = FALSE, atransf = FALSE, targs = NULL, rows = NULL,
  efac = 1, pch = 15, psize = 1, col = "darkgrey",
  border = "darkgrey", cex = NULL, cex.lab = NULL, cex.axis = NULL,
  refcol = "red", predcol = refcol, ...)
}
\arguments{
  \item{x}{An object of class uwmwRes or uwmwEstimate.}

  \item{...}{parameters passed down to the internal
  functions. These can be any of the following.}

  \item{estimate}{An optional character string defining
  which measure should be plotted. It can take the values
  \code{logor}, \code{or}, \code{logodds}, \code{odds} or
  \code{p} (for displaying the probability of differential
  expression, not the p-value!). Defaults to logor. Note
  that this argument is ignored when plotting an
  uwmwEstimate object.}

  \item{annotate}{A logical value indicating whether the
  plot needs to be annotated, i.e. whether the values for
  the chosen measure and confidence interval should be
  displayed on the right of the plot. Defaults to TRUE}

  \item{addfit}{A logical value indicating whether the
  reference measure should be plotted. See Details.}

  \item{xlim}{The horizontal limits of the plot region. If
  unspecified, the function tries to set the horizontal
  plot limits to some sensible values. Should not be used
  by the user.}

  \item{alim}{the actual x axis limits. If unspecified, an
  educated guess is taken by the function.}

  \item{ylim}{The vertical limits of the plot. If
  unspecified, the function does what it thinks is best.
  Should not be used by the user.}

  \item{at}{Position of the x axis tick marks and
  corresponding labels are placed. If unspecified, the
  function tries to position them at sensible values.}

  \item{steps}{An integer indicating the number of tick
  marks on the X axis. Ignored when \code{at} is specified.
  Defaults to 5}

  \item{level}{Numerical value between 0 and 1 to specify
  the width of the confidence interval. Defaults to 0.95
  (95\% confidence interval).}

  \item{digits}{integer specyfing the number of decimal
  places for tick mark labels and annotations. Can also be
  a vector of two integers. In that case, the first value
  specifies the number of decimal places for the
  annotations, the second for the x axis labels.}

  \item{refline}{numerical value indicating where a
  reference line should be drawn. An NA value will prevent
  the line from being drawn. See Details.}

  \item{xlab}{title for the x axis. If unspecified, the
  function tries to figure out the fitting title.}

  \item{slab}{optional vector with names for the displayed
  genes.}

  \item{mlab}{optional character string giving a label to
  the intercept estimate. If unspecified, this is created
  in the function if necessary.}

  \item{ilab}{optional vector or matrix with character
  strings providing additional information that can be
  plotted next to the genes.}

  \item{ilab.xpos}{Vector of numerical values specifying
  the x axis positions of the character vectors given via
  ilab. This has to be specified when ilab is specified.}

  \item{ilab.pos}{integer from 1 to 4 specifying the
  alignment of the character vector(s) given via ilab (2 is
  right aligned, 4 is left aligned). Default is to center
  the labels.}

  \item{order}{optional character string, character vector
  or numerical vector specifying how the genes should be
  ordered.}

  \item{transf}{optional argument specifying the name of a
  function that should be used to transform the observed
  effect sizes, summary estimates, fitted values and
  confidence interval bounds (e.g., \code{transf=exp}).
  Defaults to \code{FALSE}, which means that no
  transformation is used.}

  \item{atransf}{optional argument specifying the name of a
  function that should be used to transform the x-axis
  labels and annotations (e.g., transf=exp). Defaults to
  \code{FALSE}, which means that no transformation is
  used.}

  \item{targs}{optional arguments needed by the function
  specified via transf or atransf.}

  \item{rows}{optional vector specifying the horizontal
  position for the plotted results. If unspecified, the
  layout happens automatically. See Details and Examples
  for more information.}

  \item{efac}{numerical value specifying the vertical
  expansion of the arrows, summary estimate symbols and Ci
  limits. Normally the default of 1 should work just fine.}

  \item{pch}{plotting symbol used for the observed effect
  sizes. By default, it's a filled square.}

  \item{psize}{optional vector with the point sizes for the
  observed effects. If set to \code{NULL}, the point sizes
  are drawn proportional to the value of the log of the
  test statistic. Defaults to 1.}

  \item{cex}{optional numerical value for expansion of text
  and symbols. See also \code{\link{par}}.}

  \item{cex.lab}{Optional numerical value for expansion of
  the axis title.}

  \item{cex.axis}{Optional numerical value for expansion of
  the x axis labels.}

  \item{col}{character string specifying the color used for
  the individual estimates.}

  \item{border}{character string specifying the color used
  for the border of the individual estimates.}

  \item{refcol}{Character string specyfying the color of
  the reference line. Defaults to red.}

  \item{predcol}{character string specifying the color of
  the estimated reference value. Ignored if estimate is
  \code{logor} or \code{or}.}
}
\value{
NULL invisibly
}
\description{
This function creates a forest plot indicating the (log)
odds ratios, the (log) odds or the probabilities for the
results of the unified Wilcoxon-Mann-Whitney test.
}
\details{
The function has methods for uwmwRes and uwmwEstimate
objects. When called for an uwmwRes object, the requested
estimate is first calculated using
\code{\link{getEstimate}} and the result is passed on to
the next method.

Note that in either case, it is not possible to use the
function on a subset of either type of object. The
subsetting functions for \code{uwmwRes} and
\code{uwmwEstimate} objects return matrices, and hence
necessary information on the reference value is lost. To
plot a subset of your data, use the \code{order} argument
as shown in the examples.

Adding a reference value to the plot only makes sense when
plotting the log(odds), odds or probabilities. If log(OR)
or OR are plotted, \code{addfit} is set to FALSE.

The default settings plot a reference line at a location
depending on the plotted estimate. For log(OR), the line is
plotted at \code{refline = 0}. For OR, the line is plotted
at \code{refline = 1}.
}
\note{
Thanks to the work of W.Viechtbauer, \code{forestplot}
provides many possibilities for tweaking and customizing
the plots. Many of the arguments work the same as in the
function \code{forest.rma} (\code{metafor} package). You
can always check the help file of \code{forest.rma} for
more illustrations on the different arguments.

This function is currently implemented using an internal
function that expects an uwmwEstimate object. In a future
version, the internal function will be rewritten to be more
generic. This will enable the definition of methods for
other classes without need to change the internal function
itself.
}
\section{Warning}{
  Although the internal function is shown here (merely for
  illustration of the arguments and defaults), the user
  shouldn't be calling this one directly. The function is
  not exported.
}
\examples{
data(NBmat)
NBtest <- uWMW(NBmat, groups=NBgroups)
sigid <- which(pval(NBtest) < 0.05)
forestplot(NBtest,"logodds",order=sigid)

nameid <- c("hsa-mir-30a-3p","hsa-mir-30a-5p")
forestplot(NBtest,"p",order=nameid,addfit=FALSE,
           refline=NA,main="Comparison 30a")

forestplot(NBtest,"p",order=nameid,addfit=FALSE,
           refline=0.5,main="Comparison 30a",
           alim=c(0,1),xlim=c(-1,2),at=c(0,0.5,1))
}
\author{
This code is adapted by Joris Meys from the function
\code{forest.rma} (\code{metafor} package). The original
function is written by W. Viechtbauer.
}

