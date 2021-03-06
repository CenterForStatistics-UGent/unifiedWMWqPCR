% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMethoduwmwRes.R
\docType{methods}
\name{plot,uwmwRes,ANY-method}
\alias{plot,uwmwRes,ANY-method}
\title{Quick forest plot of significantly up- and downregulated features.}
\usage{
\S4method{plot}{uwmwRes,ANY}(x, y, ...)
}
\arguments{
\item{x}{an uwmwRes object.}

\item{y}{ignored for uwmwRes objects}

\item{...}{arguments passed down to \code{\link{forestplot}}}
}
\value{
NULL invisibly
}
\description{
This function plots a forest plot, normally used in meta analysis, to visualize the odds ratios (OR) and confidence intervals resulting from a call to \code{\link{uWMW}}. It is the default plot function for the \code{\link{uwmwRes}} objects. The results are ordered according to significance, and only the significant results are plotted. This function calls \code{\link{forestplot}} directly.
}
\examples{
data(NBmat)
NBtest <- uWMW(NBmat, groups=NBgroups)
plot(NBtest)

}
\seealso{
\code{\link{forestplot}}
}
