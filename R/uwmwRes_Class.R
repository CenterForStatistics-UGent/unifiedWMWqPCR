#' Class uwmwRes
#' 
#' This class represents the results of the unified Wilcoxon-Mann-Whitney test. It contains all necessary information for the vulcano and forest plots. For this class a number of methods is foreseen, among which accessors for every slot.
#' 
#' @section Slots:
#' \describe{
#'    \item{\code{type}:}{object of class \code{"character"}, containing the type of analysis (either H or O, see \code{\link{uWMW}} for more details.)}
#'    \item{\code{housekeeping}:}{object of class \code{"character"}, containing either NULL or the names of the housekeeping features used in the HME version of \code{\link{uWMW}}.}
#'    \item{\code{names}:}{object of class \code{"character"}, containing the names of the features that were used in the test.}
#'    \item{\code{logOR}:}{object of class \code{"numeric"}, containing the estimated log odds ratios from the uWMW test.}
#'    \item{\code{se}:}{object of class \code{"numeric"}, containing the standard errors on the estimated log odds ratios.}
#'    \item{\code{OR}:}{object of class \code{"numeric"}, containing the odds ratios estimated by the uWMW test. This slot is accessed using the function \code{oddsRatio()}}
#'    \item{\code{z.value}:}{object of class \code{"numeric"}, containing the Z values related to the odds ratios estimated by the uWMW test. These z-values relate to the chance that a specific feature is up- or downregulated, and are used as the basis for determining the p values.}
#'    \item{\code{p.value}:}{object of class \code{"numeric"}, containing the p values related to the odds ratios estimated by the uWMW test. These p-values relate to the chance that a specific feature is up- or downregulated.}
#'    \item{\code{coef}:}{object of class \code{"numeric"}, containing the estimated coefficient of the PIM model that's used in the uWMW test.}
#'    \item{\code{vcov}:}{matrix of class \code{"numeric"}, containing the variance-covariance matrix related to the estimated coefficients.}  
#'    \item{\code{id}:}{vector of class \code{"numeric"}, containing the sorting order of the features. This slot is set using the function \code{sort}}  
#'    \item{\code{orderedBy}:}{character value, indicating whether the object contains an order and if so, based on which slot. Possible values are "none", "p", "z", "or", "se" or "name". Defaults to "none".}
#'    \item{\code{groupinfo}:}{character vector of length 2, indicating the groups. This slot is mainly used to show how the probabilistic indices are calculated.}  
#' }
#' @note For this class, \code{\link{show}} and \code{\link{length}} methods are defined. \code{\link{length}} will give you the number of features.
#' 
#' @seealso \code{\link{uwmw_Accessors}}
#'  
#' @name uwmwRes-class
#' @rdname uwmwRes-class
#' @aliases uwmwRes 
#' @exportClass uwmwRes
#' @author Joris Meys

setClass("uwmwRes",
         representation=list(type         = 'character',
                             housekeeping = 'character',
                             names        = 'character',
                             logOR        = 'numeric',
                             se           = 'numeric',
                             OR           = 'numeric',
                             z.value      = 'numeric',
                             p.value      = 'numeric',
                             coef         = 'numeric',
                             vcov         = 'matrix',
                             id           = 'numeric',
                             orderedBy      = 'character',
                             groupinfo    = 'character'
                          ),
         prototype=list(type         = character(0),
                        housekeeping = character(0),
                        names        = character(0),
                        logOR        = numeric(0),
                        se           = numeric(0),
                        OR           = numeric(0),
                        z.value      = numeric(0),
                        p.value      = numeric(0),
                        coef         = numeric(0),
                        vcov         = matrix(numeric(0),ncol=0,nrow=0),
                        id           = numeric(0),
                        orderedBy    = "none",
                        groupinfo    = character(0)
                        ),
         validity=function(object){
           
           type <- object@type
           housekeeping <- object@housekeeping
           
           if(!type %in% c("H","O")) return("type should be either H or O")
           
           # Check the lengths
           ns <- slotNames(object)
           lengths <- sapply(ns[c(3:9,11)],function(i)length(slot(object,i)))
           dims <- dim(object@vcov)

           if(type=='H'){
             # You have one coef more than features in the case of HME.
             lengths["coef"] <- lengths["coef"] - 1
             dims <- dims - 1
             if(length(housekeeping)==0)
                return("Housekeeping features are not determined")
           } else if(type=='O'){
             if(length(housekeeping)!=0)
               return("Type O does not allow housekeeping features.")
           }
           
           if(length(object@groupinfo) != 2)
             return("uwmwRes must have 2 groups")
           
           if(length(unique(c(lengths,dims)))==1) TRUE else
                 "The dimensions of the slots are incompatible" 
         }
)

# Constructor method

uwmwRes <- function(type,housekeeping,names,logOR,se,OR,
                    z.value,p.value,coef,vcov,
                    id=seq_len(length(names)),
                    orderedBy=c("none","or","p","se","name","z"),
                    groupinfo){
  orderedBy <- match.arg(orderedBy)
  new("uwmwRes",
      type=type,
      housekeeping=housekeeping,
      names=names,
      logOR=logOR,
      se=se,
      OR=OR,
      z.value=z.value,
      p.value=p.value,
      coef=coef,
      vcov=vcov,
      id=id,
      orderedBy=orderedBy,
      groupinfo=groupinfo)
}



