

#' Class for \code{terms} objects which contain information about \code{nauf} contrasts.
#'
#' When \code{\link{nauf_model.frame}} is called, a
#' \code{nauf.frame} is returned, and this object's \code{terms}
#' attribute has the (S3) class \code{nauf.terms}.  The \code{nauf.terms} object
#' has an attribute \code{nauf.info} which contains all of the information
#' necessary to implement \code{\link{nauf_contrasts}} in regression.
#'
#' The \code{nauf.info} attribute is a list with the following elements:
#' \describe{
#'   \item{resp}{The name of the response variable.}
#'   \item{groups}{A named list of random effects grouping factor levels.}
#'   \item{uf}{A named list with an elment for each unordered factor. Each of
#'     these elements is a list of character vectors indicating the names of the
#'     levels that correspond to the elment's number's set of contrasts; i.e.
#'     the first element represents the main effect contrasts for the factor,
#'     the second element (if present) represents \code{.c2.} contrasts, and so
#'     on (see \code{\link{nauf_contrasts}}).}
#'   \item{of}{A named list with an element for each ordered factor containing
#'     its levels and contrasts.}
#'   \item{num}{A named list with an element for each numeric vector variable
#'     containing the variables' means.}
#'   \item{mat}{A named list with an element for each matrix variable containing
#'     the variables' colmun means.}
#'   \item{extras}{A character vector giving the names of offsets, weights,
#'     mustart, etc.}
#'   \item{cc}{A list of contrast changes required as described in
#'     \code{\link{nauf_contrasts}}.  The first element is a list of the changes
#'     required for the fixed effects.  If the model has random effects, then
#'     there is an additional element for each element of the list returned by
#'     \code{\link[lme4]{findbars}}.  Each element of \code{cc} is a named list
#'     indicating which contrasts in the \code{uf} element of \code{nauf.info}
#'     should be used.}
#'   \item{hasna}{A named logical vector with an entry for each variable
#'     indicating whether nor not the variable contains \code{NA} values.}
#'   \item{ncs_scale}{The \code{ncs_scale} argument from the call to
#'     \code{\link{nauf_model.frame}}.}
#' }
#'
#' @seealso \code{\link{nauf_contrasts}} and \code{\link{nauf_model.frame}}.
#'
#' @name nauf.terms
NULL


#' Class for fitted mixed effects models with \code{nauf} contrasts.
#'
#' Models fit with \code{\link{nauf_lmer}} have class \code{nauf.lmerMod}
#' (inheriting from \code{\linkS4class{lmerMod}}) and models fit with
#' \code{\link{nauf_glmer}} and \code{\link{nauf_glmer.nb}} have class
#' \code{nauf.glmerMod} (inheriting from \code{\linkS4class{glmerMod}}).
#'
#' @slot rsp,Gp,call,frame,flist,cnms,lower,theta,beta,u,devcomp,pp,optinfo See
#'   \code{\linkS4class{merMod}}.
#'
#' @section Generic Methods:
#' There are S3 methods specific to the
#' \code{nauf.lmerMod} and \code{nauf.glmerMod} classes for the generic
#' functions \code{\link[=predict.nauf.merMod]{predict}} and
#' \code{\link[=anova.nauf.merMod]{anova}}, which behave differently from the
#' methods for \code{\linkS4class{merMod}} objects. Other methods for generic
#' functions from the \code{stats}, \code{MASS}, and \code{lme4} packages
#' should work as they would for \code{\linkS4class{merMod}} objects, with the
#' exception of the method for the \code{\link[stats]{simulate}} function,
#' for which the \code{nauf} method has more limited options than
#' \code{\link[lme4]{simulate.merMod}}, and is not meant for end user use
#' at this time (the method is necessary for the \code{PB} method to work
#' for \code{\link[=anova.nauf.merMod]{anova}}). If you encounter a generic
#' function in these packages which does not function properly, please report
#' the issue at \url{https://github.com/CDEager/nauf/issues}.
#'
#' @seealso \code{\link{nauf_glmer}}, \code{\link{nauf_contrasts}}, and
#'   \code{\linkS4class{merMod}}.
#'
#' @name nauf.merMod
NULL


#' @rdname nauf.merMod
nauf.lmerMod <- setClass("nauf.lmerMod", contains = "lmerMod")


#' @rdname nauf.merMod
nauf.glmerMod <- setClass("nauf.glmerMod", contains = "glmerMod")


#' Class for fitted Bayesian models with \code{nauf} contrasts.
#'
#' ADD DESCRIPTION
#'
#' @seealso \code{\link{nauf_stan_glm}}, \code{\link{nauf_stan_glmer}},
#'   \code{\link{nauf_contrasts}}, and \code{\link[rstanarm]{stanreg-objects}}.
#'
#' @name nauf.stanreg
NULL


# a is array[iteration, chain, coefficient]
# x is vector[feature] or matrix[obs, feature]
# returns array[iteration, chain, prediction]
`%*a%` <- function(x, a) {
  x <- t(x)
  
  d <- dim(a)
  d[3] <- ncol(x)
  
  ans <- array(dim = d)
  for (chain in 1:d[2]) {
    ans[, chain, ] <- a[, chain, ] %*% x
  }
  
  class(ans) <- c("nauf.mcmc", "array")
  
  return(ans)
}

