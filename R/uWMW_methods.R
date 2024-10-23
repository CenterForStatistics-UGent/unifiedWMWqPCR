#' The unified Wilcoxon-Mann-Whitney test for qPCR data
#'
#' This function carries out the unified Wilcoxon-Mann-Whitney test for qPCR data. See De Neve et al. (2013) for more details.
#'
#' This function carries out the unified Wilcoxon-Mann-Whitney test using either Overall normalization (O) or Housekeeping normalization (H) as reference (see De Neve et al., 2013). If the argument \code{housekeeping.names} is specified, housekeeping normalization is considered. Otherwise overall normalization is considered.
#'
#' The function uWMW can deal with data frames and matrices. When using a data frame, you need to specify the arguments \code{groups}, \code{feat.names}, \code{subjects} and \code{value}; each one should contain the name of the related variable in the data frame.
#'
#' When using a matrix, each column is assumed to be a subject and each row a feature. The argument \code{groups} should contain as much values as there are columns in the matrix.
#'
#' @param x An object containing the qPCR measurements. See details.
#' @param groups A vector indicating the groups that need comparing, or a single character telling which variable in the data frame contains the groups. Make sure this vector is as long as the number of replicates in the data set.
#' @param feat.names An (optional) character vector with the names of the features (typically genes or microRNAs) or a single character giving the name of the feature variable. If not specified, the feature names are derived from the row names of the matrix, or from the feature names of the qPCRset object.
#' @param subjects An (optional) character string indicating which variable of the data frame contains the subject id's. Ignored if x is not a data frame.
#' @param value An (optional) character string indicating which variable of the data frame contains the values. Ignored if x is not a data frame.
#' @param housekeeping.names an (optional) vector with the names of one or more housekeeping features. Make sure those names are spelled exactly as in the object.
#' @param transpose In case a matrix is used, should the matrix be transposed? A matrix needs to be transposed when the columns do not represent the replicates. The function expects the columns to be replicates and the rows to be the different features.
#' @param ... For passing arguments between methods and to internal functions.
#'
#' @return An object of the class \code{\link{uwmwRes}}, containing the results of the unified Wilcoxon-Mann-Whitney test. See the help page of the class \code{\link{uwmwRes}} for more information.
#'
#' @author Wrapper methods are written by Joris Meys. Internal functions are written by Jan De Neve.
#'
#' @seealso \code{\link{uwmw_Accessors}} and \code{\link{uwmw_Extract}} for accessing the results, and \code{\link{volcanoplot}} and \code{\link{forestplot}} for plotting them.
#'
#' @references De Neve, J. Thas, O. Ottoy, J.P. and Clement L. (2013) An extension of the Wilcoxon-Mann-Whitney test for analyzing RT-qPCR data. Statistical Applications in genetics and Molecular Biology. 12, 333-346
#'
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#'
#' data(NBdata)
#' NBtest <- uWMW(x = NBdata, groups = "group", sub="subject",feat="miRNA",val="Cq")
#' head(NBtest)
#' as.matrix(NBtest)
#'
#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @rdname uWMW-methods
#' @aliases uWMW
#' @aliases uWMW,matrix-method
setMethod("uWMW", signature="matrix",
  function(x,groups,housekeeping.names=NULL,transpose=FALSE,feat.names=NULL){
    if(transpose) x <- t(x)
    if(!is.null(feat.names)){
      ngn <- length(feat.names)
      if(ngn != nrow(x)) stop("Length feat.names doesn't match the data")
      rownames(x) <- feat.names
    }
    .uWMW(x,groups,housekeeping.names)
  }
  )

#' @rdname uWMW-methods
#' @aliases uWMW,data.frame-method
setMethod("uWMW",signature="data.frame",
  function(x,groups,feat.names,subjects,value,...){
    # check all arguments
    if(any(missing(groups),
           missing(feat.names),
           missing(value),
           missing(subjects)))
      stop("When using uWMW on a data frame, the arguments groups, feat.names, value and subjects have to be specified.")

    mat <- prep.df(x,feat.names,subjects,value)

    id <- match(colnames(mat),x[,subjects])
    uWMW(mat,groups=x[id,groups],...)


  })
