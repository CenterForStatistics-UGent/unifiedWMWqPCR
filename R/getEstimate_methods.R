#' Extract logor, or, odds and percentages from a uwmwRes object
#'
#' With this function you can extract the modelled (log) odds ratio, odds or
#' percentages that represent the chance on differential expression as
#' estimated by the uWMW function. It also allows to extract either the
#' standard error of or the confidence interval around the estimates.
#' See Details for more explanation.
#'
#' The function can only calculate standard errors for the log OR and
#' the log odds. In all other cases, \code{se.fit} is ignored. The function
#' takes into account a possible ordering in the object
#' (see also \code{\link[unifiedWMWqPCR]{sort}}). So take into account that you
#' get the estimates in the specified order. In case you want this different,
#' either use the function \code{\link{unorder}} on the object first,
#' or check if any of the \code{\link{uwmw_Accessors}} can help you out.
#'
#' The argument \code{se.fit} is mainly to be used to save calculation time.
#' Normally there's no need to set it to \code{FALSE}.
#'
#' @param x an object of the clas uwmwRes
#' @param esttype a character string indicating the measure you want to extract.
#' It can take the values \code{logor} for the log odds ratio, \code{or} for
#' the odds ratios, \code{logodds} for the log odds, \code{odds} for the odds
#' or \code{p} for the percentages.
#' @param se.fit logical value indicating whether the standard errors of the
#' logor or the log odds should be returned as well. Ignored when type has a
#' value different from \code{logor} or \code{logodds}.
#' Note that you can also use the accessor \code{\link{se}} to
#' get only the standard errors.
#' @param ci numerical value indicating the confidence interval
#' (0.95 is 95\% confidence interval). If set to TRUE, the 95% confidence
#' interval is returned. If set to NULL, no confidence interval is returned.
#' @param drop a logical value. If set to \code{TRUE} and neither \code{se} nor
#' \code{ci} is calculated, the function returns a vector instead of
#' an uwmwEstimate object.
#' @param ... passes on arguments to the next method
#'
#' @return In general, a \code{\link{uwmwEstimate}} object with the
#' requested estimate. See \code{\link{uwmwEstimate}} for details.
#' In case \code{drop=TRUE} and neither the standard error nor
#' the confidence interval is calculated, a numeric named vector.
#'
#' @examples
#' data(NBmat)
#' NBtest <- uWMW(NBmat, groups=NBgroups)
#' getEstimate(NBtest,'logodds')
#' getEstimate(NBtest,'odds',ci=0.9)
#'
#' @include uwmwRes_Class.R uwmwEstimate_Class.R
#' @name getEstimate
#' @rdname getEstimate-methods
#' @aliases getEstimate,uwmwRes-method
setMethod("getEstimate",signature(x="uwmwRes"),
  function(x,esttype=c("logor","or","logodds","odds","p"),se.fit=(esttype %in% c("logor","logodds")),ci=TRUE,drop=TRUE){


    esttype <- match.arg(esttype)
    if(isTRUE(ci)) ci <- 0.95
    if(identical(FALSE,ci)) ci <- NULL

    xlogor <- logor(x,ordered=FALSE)
    xse    <- se(x,ordered=FALSE)
    nx     <- length(x)
    names.x <- names(x)
    xvcov <- vcov(x)

    est <- if(esttype=='logor')
              xlogor
           else if(esttype=='or')
              exp(xlogor)
           else if(esttype=='logodds')
              xlogor + coef(x)[1]
           else if(esttype=='odds')
              exp(xlogor + coef(x)[1])
            else
              1/(1 + exp(-xlogor - coef(x)[1]))


    se <- if(se.fit & esttype=="logor")
            xse
          else if(se.fit & esttype=="logodds"){

            covar <- if(type(x)=="O")
              c(xvcov[1,-1], -sum(xvcov[1,-1]))
            else
              xvcov[1,-1]

            sqrt(xse^2 + xvcov[1,1] + 2* covar)
          }
          else rep(NA,nx)

# Calculate the CI. mp is also used further on (for reference values)
    mp <- qnorm(ci+(1-ci)/2)

    if(!is.null(ci)){
      ul.logor <- (xlogor + mp * xse)
      ll.logor <- (xlogor - mp * xse)

      ul <- if(esttype=='logor')
        ul.logor
      else if(esttype=='or')
        exp(ul.logor)
      else if(esttype=='odds')
        exp(ul.logor + coef(x)[1])
      else if(esttype=="logodds")
        ul.logor + coef(x)[1]
      else if(esttype=="p")
        1/(1 + exp(-ul.logor - coef(x)[1]))
      else NA

      ll <- if(esttype=='logor')
        ll.logor
      else if(esttype=='or')
        exp(ll.logor)
      else if(esttype=='odds')
        exp(ll.logor + coef(x)[1])
      else if(esttype=="logodds")
        ll.logor + coef(x)[1]
      else if(esttype=="p")
        1/(1 + exp(-ll.logor - coef(x)[1]))
      else NA
    } else {
      ul <- ll <- as.numeric(rep(NA,nx))
    }

# Calculate the reference value for log odds and odds.
    if(!esttype %in% c("logor","or")){

      ref.logodds <- coef(x)[1]
      ref.se.logodds <- sqrt(xvcov[1,1])


      ref.ci.ub.logodds <- ref.logodds + mp*ref.se.logodds
      ref.ci.lb.logodds <- ref.logodds - mp*ref.se.logodds

      if(esttype=="odds"){
        ref <- exp(ref.logodds)
        ref.se <- NA
        ref.ci.lb <- exp(ref.ci.lb.logodds)
        ref.ci.ub <- exp(ref.ci.ub.logodds)
      } else if (esttype=="p"){
        ref <- plogis(ref.logodds)
        ref.se <- NA
        ref.ci.lb <- plogis(ref.ci.lb.logodds)
        ref.ci.ub <- plogis(ref.ci.ub.logodds)
      } else if (esttype=="logodds"){
        ref <- ref.logodds
        ref.se <- ref.se.logodds
        ref.ci.lb <- ref.ci.lb.logodds
        ref.ci.ub <- ref.ci.ub.logodds
      }
    } else {
      ref <- ref.se <- ref.ci.lb <- ref.ci.ub <- NA
    }

    id <- getOrder(x)

    if(!se.fit & is.null(ci) & drop){
      out <- est[id]
      names(out) <- names.x
      out
    } else {
    uwmwEstimate(
      esttype=esttype,
      names=names.x,
      est = est[id],
      se = se[id],
      ll = ll[id],
      ul = ul[id],
      refest = ref,
      refse = ref.se,
      reful = ref.ci.ub,
      refll = ref.ci.lb,
      type = type(x),
      confint = ci,
      housekeeping = housekeeping(x),
      groupinfo = groupinfo(x)
      )
    }

  })
