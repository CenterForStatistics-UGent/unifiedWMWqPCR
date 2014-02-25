#' Documentation for the dataset NBdata
#' 
#' The example data used in this package, are a subset of the data provided by Mestdagh et al. (2009). The subset contains quantification cycles of 323 microRNAs in 61 neuroblastoma (NB) tumor samples: 22 MYCN amplified (called MNA) and 39 MYCN single copy samples (called MNSC). The subset was selected so that all microRNAs with a least 85% undetermined values in both groups were removed, see De Neve et al. (2013) for details.
#' 
#' The data exists in different formats. \code{NBdata} gives you a data frame with following variables:
#' \describe{
#'  \item{\code{subject}:}{variable of class \code{"factor"}, indicating the subject code.}
#'  \item{\code{miRNA}:}{variable of class \code{"factor"}, indicating the miRNA code.}
#'  \item{\code{Cq}:}{variable of class \code{"numeric"}, containing the cycle information.}
#'  \item{\code{group}:}{variable of class \code{"factor"}, indicating the group code.}
#' }
#' 
#' The data matrix \code{NBmat} contains the same data in matrix format, where the rows are the different miRNA's and the columns the different subjects. The vector \code{NBgroups} specifies to which group every column of \code{NBmat} belongs.
#' 
#' @references Mestdagh, P., P. Van Vlierberghe, A. De Weer, D. Muth, F. Westermann, F. Speleman, and J. Vandesompele (2009) A novel and universal method for microRNA RT-qPCR data normalization. Genome Biology., 10, R64.
#' 
#' De Neve, J. Thas, O. Ottoy, J.P. and Clement L. (2013) An extension of the Wilcoxon-Mann-Whitney test for analyzing RT-qPCR data. Statistical Applications in Genetics and Molecular Biology. 12, 333-346.
#' 
#' @examples
#' # Look at the data frame
#' data(NBdata)
#' str(NBdata)
#' 
#' # Look at the matrix and grouping vector
#' data(NBmat)
#' str(NBmat)
#' str(NBgroups)
#'  
#' @name NBdata
#' @rdname NBdata
#' @aliases NBmat NBgroups
#' @docType data
#' @keywords data
NULL