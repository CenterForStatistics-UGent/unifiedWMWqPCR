\docType{methods}
\name{as.matrix.uwmwRes}
\alias{as.matrix,uwmwRes-method}
\alias{as.matrix.uwmwRes}
\title{S3 method as.matrix for uwmwRes objects.}
\usage{
\method{as.matrix}{uwmwRes}(x, ...)

\S4method{as.matrix}{uwmwRes}(x, ...)
}
\arguments{
  \item{x}{a uwmwRes object}

  \item{...}{currently ignored}
}
\value{
A matrix containing the columns specified above.
}
\description{
For \code{\link{uwmwRes}} objects, an \code{as.matrix}
method is defined that transforms the object to a numeric
matrix with the following columns: \describe{
\item{\code{logor}:}{ The log odds ratio values}
\item{\code{se}:}{ The standard errors on the log OR
values.} \item{\code{or}}{The odds ratio values}
\item{\code{z.value}}{the z values related to the log OR
values} \item{\code{p.value}}{The p values related to the
log OR values} } The row names are the names of the tested
genes/features. The matrix takes the ordering in the object
into account if necessary.
}
\examples{
data(NBmat)
NBtest <- uWMW(NBmat, groups=NBgroups)
as.matrix(NBtest)
}

