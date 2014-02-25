\docType{methods}
\name{uwmw_Extract}
\alias{Extract}
\alias{[,uwmwEstimate,ANY,ANY-method}
\alias{[,uwmwEstimate,ANY-method}
\alias{[,uwmwEstimate,character,ANY-method}
\alias{[,uwmwEstimate,character-method}
\alias{[,uwmwRes,ANY,ANY-method}
\alias{[,uwmwRes,ANY-method}
\alias{[,uwmwRes,character,ANY-method}
\alias{[,uwmwRes,character-method}
\alias{uwmw_Extract}
\title{Extract data from uwmwRes and uwmwEstimate objects.}
\usage{
\S4method{[}{uwmwRes,character}(x, i, j, drop = TRUE)

\S4method{[}{uwmwEstimate,character}(x, i, j, drop = TRUE)

\S4method{[}{uwmwRes,ANY}(x, i, j, drop = TRUE)

\S4method{[}{uwmwEstimate,ANY}(x, i, j, drop = TRUE)
}
\arguments{
  \item{i}{numeric or character vector}

  \item{j}{optional numeric or character vector}

  \item{drop}{logical value. If set to FALSE and the result
  of the extraction is a single row or column, dimensions
  are dropped. Defaults to \code{TRUE}.}

  \item{x}{ object from which to extract element(s) or in
  which to replace element(s).  }
}
\value{
mostly a matrix, unless \code{drop=TRUE} and a single row
or column is selected. Then a vector.
}
\description{
For both uwmwRes and uwmwEstimate objects, you can use the
square bracket operators to extract information much like
you would do for a matrix.
}
\details{
The features can be extracted using the row index, and the
estimates as
}
\examples{
# With an uwmwRes object
data(NBmat)
NBtest <- uWMW(NBmat,NBgroups)
# These two lines are the same
NBtest["hsa-mir-1"]
NBtest["hsa-mir-1", ]
# These two not
str(NBtest["hsa-mir-1",,drop=FALSE])
str(NBtest["hsa-mir-1",])
# These two give the same data, but in a different way:
se(NBtest) # unnamed
NBtest[,"se"] #

# With an uwmwEstimate object
NBodds <- getEstimate(NBtest,"odds")
gnames <- grep("let",names(NBodds),value=TRUE)
NBodds[gnames]
NBodds[gnames,c("ll","ul")]
}

