\name{volcanopar}
\alias{volcano-par}
\alias{volcanopar}
\alias{volcanoplot-par}
\alias{volcanoplotpar}
\title{Parameters used for the function volcanoplot}
\arguments{
  \item{highlight}{an integer, indicating how many of the
  top features should be highlighted. Defaults to 0.}

  \item{names}{an optional character vector giving the
  names of the genes to be highlighted. If not specified,
  the features are numbered from most significant to least
  significant. Default value is \code{NULL}.}

  \item{xlab}{see \code{\link{par}}. Defaults to
  \code{"Fit"}}

  \item{ylab}{see \code{\link{par}}. Defaults to
  \code{"Significance"}}

  \item{pch}{see \code{\link{par}}. Defaults to 16}

  \item{cex}{see \code{\link{par}}. Defaults to 0.35 as in
  the \code{limma} package.}

  \item{transf.x}{a function used to transpose the values
  on the X axis. Defaults to NULL, unless otherwise
  specified in the specific methods.}

  \item{transf.y}{a function used to transpose the values
  on the Y axis. Defaults to NULL, unless otherwise
  specified in the specific method.}

  \item{add.ref}{a charater value indicating if reference
  lines for x, y or both axes should be drawn. It takes the
  values "none", "x", "y" or "both". Default value is
  "none". If not specified but ref.x or ref.y is, then the
  reference line is drawn for the X and/or the Y axis.}

  \item{ref.x}{a numerical vector, indicating the value or
  values at which the reference line(s) should be drawn on
  the X axis. Defaults to \code{NULL}.}

  \item{ref.y}{a numerical value, indicating the value at
  which the reference line should be drawn on the Y axis.
  Defaults to \code{NULL} unless otherwise specified in the
  specific methods.}

  \item{col.x}{The color for the reference line on the X
  axis. The argument is passed to the col argument of
  \code{\link{lines}}, so check that help page for possible
  values. Defaults to \code{"darkgrey"}}

  \item{col.y}{see col.x, but for the Y axis. If
  \code{col.x} is specified and \code{col.y} is not, then
  the value of \code{col.x} is used for \code{col.y}.}

  \item{lwd.x}{The line width for the reference line on the
  X axis. The argument is passed to the lwd argument of
  \code{\link{lines}}, so check that help page for possible
  values. Defaults to 1.}

  \item{lwd.y}{see \code{lwd.x}, but for the Y axis. If
  \code{lwd.y} is not specified and \code{lwd.x} is, that
  value is used for \code{lwd.y} as well.}

  \item{lty.x}{The line type for the reference line on the
  X axis. The argument is passed to the lty argument of
  \code{\link{lines}}, so check that help page for possible
  values. Defaults to 1.}

  \item{lty.y}{see \code{lty.x}, but for the Y axis. If
  \code{lty.y} is not specified and \code{lty.x} is, that
  value is used for \code{lty.y} as well.}
}
\description{
The parameters described can all be used in the different
\link{volcanoplot} methods.
}
\details{
The arguments \code{highlight} and \code{names} are chosen
in such a way that the function can mimick the behaviour of
the \code{volcanoplot} function in the package
\code{limma}.
}
\note{
These parameters belong to the internal function,
eventually called by the different S4 methods. This
internal function is not exported, and should not be called
directly.
}
\examples{
# see the help page of volcanoplot
}

