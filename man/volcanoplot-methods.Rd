\docType{methods}
\name{volcanoplot}
\alias{volcanoplot}
\alias{volcanoplot,MArrayLM-method}
\alias{volcanoplot,matrix-method}
\alias{volcanoplot,numeric-method}
\alias{volcanoplot,uwmwRes-method}
\title{Make a volcano plot of the outcome of a uWMW test}
\usage{
\S4method{volcanoplot}{uwmwRes}(fit, esttype = c("logor", "or", "logodds",
  "odds", "p"), transf.y = function(i) -log10(i), ref.y = -log10(0.05), ...)

\S4method{volcanoplot}{matrix}(fit, ...)

\S4method{volcanoplot}{numeric}(fit, ...)

\S4method{volcanoplot}{MArrayLM}(fit, coef = 1, highlight = 0, ...)
}
\arguments{
  \item{fit}{a uwmwRes object, a matrix with 2 columns or a
  numeric vector with either odds ratios or the log of odds
  ratios.}

  \item{esttype}{The estimate that should be plotted from
  the uwmwRes object. See also \code{\link{getEstimate}}.}

  \item{transf.y}{a function used to transpose the values
  on the Y axis. For uwmwRes objects, defaults to
  -log10(y). Passed down to internal function.}

  \item{ref.y}{The location on the Y axis for the reference
  line.}

  \item{coef}{for compatibility with limma package, see
  \code{\link{volcanopar}}}

  \item{highlight}{for compatibiity with limma package, see
  \code{\link{volcanopar}}}

  \item{...}{arguments passed to the internal function. See
  \code{\link{volcanopar}}}
}
\value{
invisible NULL
}
\description{
This function creates a volcano plot of the outcome of a
uWMW test. In this plot, the p-value is plotted against the
odd ratio or the log-odds ratio, depending on what
interests you most. It allows to quickly see what
proportion of results is significant, and what proportion
of significant results has a biologically significant
chance of being upregulated.
}
\details{
The methods described here all use the same internal
function to create the plot
}
\examples{
data(NBmat)
housekeeping.id <- grep("let",rownames(NBmat),value=TRUE)
NB.Htest <- uWMW(NBmat,NBgroups,housekeeping.id)

volcanoplot(NB.Htest)
volcanoplot(NB.Htest,"p")
}
\author{
Joris Meys
}

