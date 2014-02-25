#' The class uwmwEstimate
#' 
#' This class represents an estimate object resulting from a call to \code{\link{getEstimate}}. It contains all information about the estimate, including standard errors and confidence intervals if requested. For this class a number of methods is foreseen, including the accessors for the slots. The class is sortable and can be indexed, so you can use this for making custom forest plots using the function \code{\link{forestplot}}.
#' 
#' #' @section Slots:
#' \describe{
#'    \item{\code{esttype}:}{object of class \code{"character"}, containing the estimate type. This can be logor for log odds ratio, or for odds ratio, odds, logodds for the log odds or p for the probability.}
#'    \item{\code{names}:}{object of class \code{"character"}, containing the names of the genes for which the estimates are calculated.}  
#'    \item{\code{est}:}{object of class \code{"numeric"}, containing the estimates itself.}
#'    \item{\code{se}:}{object of class \code{"numeric"}, containing the estimates for the standard error, if applicable.}
#'    \item{\code{ll}:}{object of class \code{"numeric"}, containing the lower limit of the confidence interval.}
#'    \item{\code{ul}:}{object of class \code{"numeric"}, containing the upper limit of the confidence interval.}
#'    \item{\code{refest}:}{object of class \code{"numeric"}, containing the estimate for the reference used in the analysis. Note that this only makes sense for log odds, odds and probabilities.}
#'    \item{\code{refse}:}{object of class \code{"numeric"}, containing the se estimate for the reference if applicable.}
#'    \item{\code{refll}:}{object of class \code{"numeric"}, containing the lower limit for the reference if applicable.}
#'    \item{\code{reful}:}{object of class \code{"numeric"}, containing the upper limit for the reference if applicable.}
#'    \item{\code{type}:}{vector of class \code{"character"}, containing the type of reference used in the original analysis. This can be either "O" or "H" for Overall respectively Housekeeping Expression as reference.} 
#'    \item{\code{confint}:}{vector of class \code{"numeric"}}, indicating the limit used for the confidence interval. 0.95 represents the 95\% confidence interval.
#'    \item{\code{housekeeping}:}{object of class \code{"character"}, containing either NULL or the names of the housekeeping genes used in the H version of \code{\link{uWMW}}.}
#'    \item{\code{groupinfo}:}{character vector of length 2, indicating the groups. This slot is mainly used to show how the probabilistic indices are calculated.} 
#' }
#' 
#' @note For this class, \code{\link{show}} and \code{\link{length}} methods are defined. \code{\link{length}} will give you the number of features.
#' 
#' @name uwmwEstimate-class
#' @rdname uwmwEstimate-class
#' @aliases uwmwEstimate 
#' @exportClass uwmwEstimate
#' @author Joris Meys

setClass("uwmwEstimate",
         representation=list(esttype      = 'character',
                             names        = 'character',
                             est          = 'numeric',
                             se           = 'numeric',
                             ll           = 'numeric',
                             ul           = 'numeric',
                             refest       = 'numeric',
                             refse        = 'numeric',
                             refll        = 'numeric',
                             reful        = 'numeric',
                             type         = 'character',
                             confint      = 'numeric',
                             housekeeping = 'character',
                             groupinfo    = 'character'
         ),
         prototype=list(esttype      = character(0),
                        names        = character(0),
                        est          = numeric(0),
                        se           = numeric(0),
                        ll           = numeric(0),
                        ul           = numeric(0),
                        refest       = numeric(0),
                        refse        = numeric(0),
                        refll        = numeric(0),
                        reful        = numeric(0),
                        type         = character(0),
                        confint      = numeric(0),
                        housekeeping = character(0),
                        groupinfo    = character(0)
         ),
         validity=function(object){
           
           esttype <- object@esttype
           type <- object@type
           refest <- object@refest
           refse <- object@refse
           refll <- object@refll
           reful <- object@reful
           
           if(!esttype %in% c("or","logor","odds","logodds","p"))
             return("esttype should be one of or, logor, odds, logodds or p.")
                              
           if(!type %in% c("H","O")) 
             return("type should be either H or O.")
           else if(type =="H" && is.null(object@housekeeping))
             return("housekeeping genes have to be provided when type is H.")
           
           if(esttype %in% c("or","logor"))
             if(!all(is.na(c(refest,refse,refll,reful))))
               return("Reference value makes no sense for odds ratio or logor.")
           
           if(length(object@confint) > 1)
             return("invalid specification for confint.")
           if(object@confint > 1 | object@confint < 0)
             return("confint should be between 0 and 1.")
           
           if(length(object@groupinfo) != 2)
             return("uwmwEstimate groupinfo must be of length 2.")
           

           ns <- slotNames(object)[2:6]
           lengths <- sapply(ns,function(i)length(slot(object,i)))
           if(length(unique(lengths))==1) TRUE else
             "The dimensions of the slots are incompatible." 
           
         }
)

# Constructor method

uwmwEstimate <- function(esttype,
                    names,
                    est,
                    se,
                    ll,
                    ul,
                    refest=NA,
                    refse=NA,
                    refll=NA,
                    reful=NA,
                    type,
                    confint,
                    housekeeping=character(0),
                    groupinfo
                         ){
  new("uwmwEstimate",
      esttype      = esttype,
      names        = names,
      est          = est,
      se           = as.numeric(se),
      ll           = ll,
      ul           = ul,
      refest       = as.numeric(refest),
      refse        = as.numeric(refse),
      refll        = as.numeric(refll),
      reful        = as.numeric(reful),
      type         = type,
      confint      = as.numeric(confint),
      housekeeping = housekeeping,
      groupinfo    = groupinfo
      )
}

