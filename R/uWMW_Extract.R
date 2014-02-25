#' Extract data from uwmwRes and uwmwEstimate objects.
#' 
#' For both uwmwRes and uwmwEstimate objects, you can use the square bracket operators to extract information much like you would do for a matrix. 
#' 
#' The features can be extracted using the row index, and the estimates as 
#' 
#' @param i numeric or character vector
#' @param j optional numeric or character vector
#' @param drop logical value. If set to FALSE and the result of the extraction is a single row or column, dimensions are dropped. Defaults to \code{TRUE}.
#' 
#' @return mostly a matrix, unless \code{drop=TRUE} and a single row or column is selected. Then a vector.
#' 
#' @examples
#' # With an uwmwRes object
#' data(NBmat)
#' NBtest <- uWMW(NBmat,NBgroups)
#' # These two lines are the same
#' NBtest["hsa-mir-1"]
#' NBtest["hsa-mir-1", ]
#' # These two not
#' str(NBtest["hsa-mir-1",,drop=FALSE])
#' str(NBtest["hsa-mir-1",])
#' # These two give the same data, but in a different way:
#' se(NBtest) # unnamed
#' NBtest[,"se"] # 
#' 
#' # With an uwmwEstimate object
#' NBodds <- getEstimate(NBtest,"odds")
#' gnames <- grep("let",names(NBodds),value=TRUE)
#' NBodds[gnames]
#' NBodds[gnames,c("ll","ul")]
#' 
#' @name uwmw_Extract
#' @rdname uwmw_Extract
#' @aliases Extract


# method [] method
#' @export
#' @rdname uwmw_Extract
#' @aliases [,uwmwRes,character,ANY-method [,uwmwRes,character-method
setMethod("[",signature(x="uwmwRes",i="character",j="ANY"),
          function(x,i,j,drop=TRUE){
            id <- match(i,names(x),0L)
            if(any(id==0)) warning("Some features could not be found in the object.")
            
            out <- x[id,j,drop=drop]
            return(out)
            
          }
)

# method [] method
#' @export
#' @rdname uwmw_Extract
#' @aliases [,uwmwEstimate,character,ANY-method [,uwmwEstimate,character-method
setMethod("[",signature(x="uwmwEstimate",i="character",j="ANY"),
          function(x,i,j,drop=TRUE){
            id <- match(i,names(x),0L)
            if(any(id==0)) warning("Some features could not be found in the object.")
            
            out <- x[id,j,drop=drop]
            return(out)
            
          }
)


# method [] method
#' @export
#' @rdname uwmw_Extract
#' @aliases [,uwmwRes,ANY,ANY-method [,uwmwRes,ANY-method
setMethod("[",signature(x="uwmwRes",i="ANY",j="ANY"),
          function(x,i,j,drop=TRUE){
            Narg <- nargs() - 1
            no.i <- missing(i)
            
            if(Narg==2L){
              out <- as.matrix(x)[i,]
            }else{
              if(no.i) i <- seq_len(length(x))
              out <- as.matrix(x)[i,j,drop=drop]
              
            }
            return(out)
          }
)

# method [] method
#' @export
#' @rdname uwmw_Extract
#' @aliases [,uwmwEstimate,ANY,ANY-method [,uwmwEstimate,ANY-method
setMethod("[",signature(x="uwmwEstimate",i="ANY",j="ANY"),
          function(x,i,j,drop=TRUE){
            Narg <- nargs() - 1
            no.i <- missing(i)
            
            if(Narg==2L){
              out <- cbind(
                y     = x@est[i],
                se    = x@se[i],
                ll    = x@ll[i],
                ul    = x@ul[i]
              )
              rownames(out) <- names(x)[i]
            }else{
              if(no.i) i <- seq_len(length(x))
              xmat <- as.matrix(x)
              out <- xmat[i,j,drop=drop]
              
            }
            return(out)
          }
)
