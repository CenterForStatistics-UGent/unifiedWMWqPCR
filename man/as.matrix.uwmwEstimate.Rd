\docType{methods}
\name{as.matrix.uwmwEstimate}
\alias{as.matrix,uwmwEstimate-method}
\alias{as.matrix.uwmwEstimate}
\title{S3 method as.matrix for uwmwEstimate objects.}
\usage{
\method{as.matrix}{uwmwEstimate}(x, ...)

\S4method{as.matrix}{uwmwEstimate}(x, ...)
}
\arguments{
  \item{x}{a uwmwEstimate object}

  \item{...}{currently ignored}
}
\value{
A matrix with the columns specified above.
}
\description{
For \code{\link{uwmwEstimate}} objects, an \code{as.matrix}
method is defined that transforms the object to a numeric
matrix with the following columns: \describe{
\item{\code{est}:}{ The estimates} \item{\code{se}:}{ The
standard errors on estimates.} \item{\code{ll}}{lower limit
of the confidence interval} \item{\code{ul}}{upper limit of
the confidence interval} } The row names are the names of
the tested genes/features. The matrix takes the ordering in
the object into account if necessary.
}
\examples{
data(NBmat)
NBtest <- uWMW(NBmat, groups=NBgroups)
NBest <- getEstimate(NBtest,"p")
as.matrix(NBtest)
}
\author{
Joris Meys
}

