% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uWMW_methods.R
\docType{methods}
\name{uWMW,matrix-method}
\alias{uWMW,matrix-method}
\alias{uWMW}
\alias{uWMW,qPCRset-method}
\alias{uWMW,data.frame-method}
\title{The unified Wilcoxon-Mann-Whitney test for qPCR data}
\usage{
\S4method{uWMW}{matrix}(x, groups, housekeeping.names = NULL,
  transpose = FALSE, feat.names = NULL)

\S4method{uWMW}{qPCRset}(x, feat.names = NULL, ...)

\S4method{uWMW}{data.frame}(x, groups, feat.names, subjects, value, ...)
}
\arguments{
\item{x}{An object containing the qPCR measurements. See details.}

\item{groups}{A vector indicating the groups that need comparing, or a single character telling which variable in the data frame contains the groups. Make sure this vector is as long as the number of replicates in the data set.}

\item{housekeeping.names}{an (optional) vector with the names of one or more housekeeping features. Make sure those names are spelled exactly as in the object.}

\item{transpose}{In case a matrix is used, should the matrix be transposed? A matrix needs to be transposed when the columns do not represent the replicates. The function expects the columns to be replicates and the rows to be the different features.}

\item{feat.names}{An (optional) character vector with the names of the features (typically genes or microRNAs) or a single character giving the name of the feature variable. If not specified, the feature names are derived from the row names of the matrix, or from the feature names of the qPCRset object.}

\item{...}{For passing arguments between methods and to internal functions.}

\item{subjects}{An (optional) character string indicating which variable of the data frame contains the subject id's. Ignored if x is not a data frame.}

\item{value}{An (optional) character string indicating which variable of the data frame contains the values. Ignored if x is not a data frame.}
}
\value{
An object of the class \code{\link{uwmwRes}}, containing the results of the unified Wilcoxon-Mann-Whitney test. See the help page of the class \code{\link{uwmwRes}} for more information.
}
\description{
This function carries out the unified Wilcoxon-Mann-Whitney test for qPCR data. See De Neve et al. (2013) for more details.
}
\details{
This function carries out the unified Wilcoxon-Mann-Whitney test using either Overall normalization (O) or Housekeeping normalization (H) as reference (see De Neve et al., 2013). If the argument \code{housekeeping.names} is specified, housekeeping normalization is considered. Otherwise overall normalization is considered.

The function uWMW can deal with data frames, matrices and \code{qPCRset} objects from the package \code{HTqPCR}. When using a data frame, you need to specify the arguments \code{groups}, \code{feat.names}, \code{subjects} and \code{value}; each one should contain the name of the related variable in the data frame.

When using a matrix, each column is assumed to be a subject and each row a feature. The argument \code{groups} should contain as much values as there are columns in the matrix. Note that the method for \code{qPCRset} objects does nothing else but extract the data matrix and do the analysis. You still need to specify where the groups are to be found.
}
\examples{
data(NBmat)
NBtest <- uWMW(NBmat, groups=NBgroups)

data(NBdata)
NBtest <- uWMW(x = NBdata, groups = "group", sub="subject",feat="miRNA",val="Cq")
head(NBtest)
as.matrix(NBtest)

}
\references{
De Neve, J. Thas, O. Ottoy, J.P. and Clement L. (2013) An extension of the Wilcoxon-Mann-Whitney test for analyzing RT-qPCR data. Statistical Applications in genetics and Molecular Biology. 12, 333-346
}
\seealso{
\code{\link{uwmw_Accessors}} and \code{\link{uwmw_Extract}} for accessing the results, and \code{\link{volcanoplot}} and \code{\link{forestplot}} for plotting them.
}
\author{
Wrapper methods are written by Joris Meys. Internal functions are written by Jan De Neve.
}
