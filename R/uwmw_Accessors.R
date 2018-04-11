#' Extract info from uwmwRes and uwmwEstimate objects
#'
#' This help file describes different ways to access the slots and values contained in \code{\link{uwmwRes}} objects resulting from calls to \code{\link{uWMW}}, and in \code{\link{uwmwEstimate}} objects resulting from calls to \code{\link{getEstimate}}.
#'
#' @param x a uwmwRes object or a uwmwEstimate object.
#' @param ordered logical value. If set to \code{FALSE}, the extracted values are given in the original order, also after a uwmwRes object has been ordered using \code{\link{sort}}. Defaults to \code{TRUE}
#' @param na.rm for compatibility with the base function \code{\link{is.unsorted}}. Ignored for uwmwRes objects.
#' @param strictly for compatibility with the base function \code{\link{is.unsorted}}. Ignored for uwmwRes objects.
#'
#' @seealso \code{\link{uwmw_Extract}} for matrix like extraction of data.
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat,NBgroups)
#' coef(NBtest)[1:10]
#' type(NBtest)
#'
#' # With a sorted object
#' NBsort <- sort(NBtest, which="p")
#' is.unsorted(NBtest)
#' is.unsorted(NBsort)
#' orderedBy(NBsort)
#'
#' # On an Estimate object
#' NBlogodds <- getEstimate(NBsort,"logodds")
#' se(NBlogodds)
#' # and so on...

#' @name uwmw_Accessors
#' @rdname uwmw_Accessors
#' @aliases type coef vcov ref logor esttype getOrder housekeeping orderedBy pval se zval oddsRatio groupinfo type,uwmwRes-method
#' @docType methods
#' @include uwmwRes_Class.R uwmwEstimate_Class.R allGenerics.R S3Generics.R
#' @importFrom stats coef vcov

#' @return \code{type(x)} returns the type of uWMW carried out (i.e. O or H for using overall respectively housekeeping expression as a reference.)
#' @rdname uwmw_Accessors
#' @name type
setMethod("type",signature(x="uwmwRes"),
          function(x) x@type)

#' @rdname uwmw_Accessors
#' @aliases type,uwmwEstimate-method
setMethod("type",signature(x="uwmwEstimate"),
          function(x) x@type)

# housekeeping accessor
#' @return \code{housekeeping(x)} returns the content of the housekeeping slot from the object. Or, in case of overall normalization, it returns NULL.
#' @rdname uwmw_Accessors
#' @aliases housekeeping,uwmwRes-method
setMethod("housekeeping",signature(x="uwmwRes"),
          function(x) x@housekeeping)
#' @rdname uwmw_Accessors
#' @aliases housekeeping,uwmwEstimate-method
setMethod("housekeeping",signature(x="uwmwEstimate"),
          function(x) x@housekeeping)

# names accessor
#' @export
#' @return \code{names(x)} returns the names of the genes in the object, and in the order defined in the object.
#' @rdname uwmw_Accessors
#' @aliases names,uwmwRes-method
setMethod("names",signature(x="uwmwRes"),
          function(x){
            if(!is.unsorted(x))
              x@names[x@id]
            else
              x@names
          } )
#' @rdname uwmw_Accessors
#' @aliases names,uwmwEstimate-method
setMethod("names",signature(x="uwmwEstimate"),
          function(x){
              x@names
          } )



# logor accessor
#' @return \code{logor(x)} returns a numeric vector with the log OR values.
#' @rdname uwmw_Accessors
setMethod("logor",signature(x="uwmwRes"),
          function(x,ordered=TRUE){
            if(ordered)
              x@logOR[x@id]
            else
              x@logOR
          } )

# se accessor
#' @return \code{se(x)} returns a numeric vector with the standard errors on the logor.

#' @rdname uwmw_Accessors
setMethod("se",signature(x="uwmwRes"),
          function(x,ordered=TRUE){
            if(ordered)
              x@se[x@id]
            else
              x@se
          } )

#' @rdname uwmw_Accessors
setMethod("se",signature(x="uwmwEstimate"),
          function(x){
              x@se
          } )

# OR accessor
#' @return \code{oddsRatio(x)} returns a numeric vector with the odds ratios.
#' @rdname uwmw_Accessors
#' @aliases oddsRatio,uwmwRes-method
setMethod("oddsRatio",signature(x="uwmwRes"),
          function(x,ordered=TRUE){
            if(ordered)
              x@OR[x@id]
            else
              x@se
          } )

# z.value accessor
#' @return \code{zval(x)} returns a numerical vector with the Z values contained in the object.
#' @rdname uwmw_Accessors
#' @aliases zval,uwmwRes-method
setMethod("zval",signature(x="uwmwRes"),
          function(x, ordered=TRUE){
            if(ordered)
              x@z.value[x@id]
            else
              x@z.value
          } )

# p.value accessor
#' @return \code{pval(x)} returns a numerical vector containing the p values in the object.
#' @rdname uwmw_Accessors
#' @aliases pval,uwmwRes-method
setMethod("pval",signature(x="uwmwRes"),
          function(x,ordered=TRUE){
            if(ordered)
              x@p.value[x@id]
            else
              x@p.value
          } )

# coef method\
#' @param object see \code{x}
#' @rdname uwmw_Accessors
#' @aliases coef,uwmwRes-method
#' @export
setMethod("coef",signature(object="uwmwRes"),
          function(object) {
            object@coef
          }
          )

# vcov method
#' @rdname uwmw_Accessors
#' @aliases vcov,uwmwRes-method
#' @export
setMethod("vcov",signature(object="uwmwRes"),
          function(object){
            object@vcov
          }
          )

# ref method
#' @rdname uwmw_Accessors
#' @aliases ref,uwmwEstimate-method
setMethod('ref',signature(x="uwmwEstimate"),
          function(x){
            c(refest=x@refest,
              se = x@refse,
              ll = x@refll,
              ul = x@reful)
          })
# esttype method
#' @export
#' @rdname uwmw_Accessors
#' @aliases esttype,uwmwEstimate-method
setMethod('esttype',signature(x='uwmwEstimate'),
          function(x){
            x@esttype
          })

# method getOrder
#' @rdname uwmw_Accessors
#' @aliases getOrder,uwmwRes-method
setMethod("getOrder",signature(x='uwmwRes'),
          function(x){
            x@id
          })

# method is.unsorted
#' @export
#' @rdname uwmw_Accessors
#' @aliases is.unsorted,uwmwRes-method
setGeneric("is.unsorted")
setMethod("is.unsorted",
          signature="uwmwRes",
          function(x,na.rm=FALSE,strictly=FALSE){
            if(na.rm || strictly){
              warning("The arguments na.rm and strictly are ignored when checking an uwmwRes object.")
            }
            orderedBy(x) == "none"
          })

# method orderedBy
#' @return \code{orderedBy(x)} returns the slot by which the object is ordered (i.e. the value of the slot orderedBy)
#' @rdname uwmw_Accessors
#' @aliases orderedBy,uwmwRes-method
setMethod("orderedBy",
          signature="uwmwRes",
          function(x) x@orderedBy
          )

#' @rdname uwmw_Accessors
#' @aliases orderedBy,uwmwEstimate-method
setMethod("orderedBy",
          signature="uwmwEstimate",
          function(x) x@orderedBy
)

# method orderedBy
#' @return \code{groupinfo(x)} returns the groupinfo slot, i.e. a character vector of length 2 that indicates in which order the groups are compared by uWMW.
#' @rdname uwmw_Accessors
#' @aliases groupinfo,uwmwRes-method
setMethod("groupinfo",
          signature="uwmwRes",
          function(x) x@groupinfo
)

#' @rdname uwmw_Accessors
#' @aliases groupinfo,uwmwEstimate-method
setMethod("groupinfo",
          signature="uwmwEstimate",
          function(x) x@groupinfo
)


# length method
#' @rdname uwmw_Accessors
#' @aliases length,uwmwEstimate-method
setMethod('length','uwmwEstimate',
          function(x){
            length(names(x))
          })

# length method
#' @rdname uwmw_Accessors
#' @aliases length,uwmwRes-method
setMethod('length','uwmwRes',
          function(x){
            length(names(x))
          })
