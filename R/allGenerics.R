#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @importFrom BiocGenerics sort
#' 
#' @export
setGeneric("forestplot",function(x,...) standardGeneric("forestplot"))

#' @export
setGeneric("volcanoplot",function(fit,...) standardGeneric("volcanoplot"))

#' @export
setGeneric("getEstimate",function(x,...) standardGeneric("getEstimate"))

#' @export
setGeneric("unorder", function(x,...) standardGeneric("unorder"))

#' @export
setGeneric("uWMW",function(x,...) standardGeneric("uWMW"))

#' @export
setGeneric("type",function(x) standardGeneric("type"))

#' @export
setGeneric('ref',function(x,...) standardGeneric('ref'))

#' @export
setGeneric("logor",function(x,...) standardGeneric("logor"))

#' @export
setGeneric('esttype',function(x,...) standardGeneric('esttype'))

#' @export
setGeneric("housekeeping",function(x) standardGeneric("housekeeping"))

#' @export
setGeneric("se",function(x,...) standardGeneric("se"))

#' @export
setGeneric("pval",function(x,...) standardGeneric("pval"))

#' @export
setGeneric("zval",function(x,...) standardGeneric("zval"))

#' @export
setGeneric("getOrder",function(x,...) standardGeneric("getOrder"))

#' @export
setGeneric("orderedBy",function(x,...)standardGeneric("orderedBy"))

#' @export
setGeneric("oddsRatio",function(x,...) standardGeneric("oddsRatio"))

#' @export
setGeneric("groupinfo",function(x,...) standardGeneric("groupinfo"))

#' @export
setGeneric("plot")
#' @export
setGeneric("coef")
#' @export
setGeneric("vcov")
#' @export
setGeneric("as.matrix")