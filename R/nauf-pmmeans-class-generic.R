

#' Predicted marginal means for \code{nauf} models.
#'
#' Create a reference grid for a \code{nauf} model with \code{nauf_ref.grid},
#' and use the resultnig \code{nauf.ref.grid} as the \code{object} argument to
#' \code{nauf_pmmeans} to obtain predicted marginal means and pairwise
#' comparisons, optionally conditioning these predictions on certain subsets of
#' the data via the \code{subset} argument.
#'
#' A reference grid creates a data frame which contains all possible
#' combinations of the factors in a regression model, holding all covariates
#' at their mean values.  There are many options for
#' \code{\link[lsmeans]{ref.grid}} which are not currently supported for
#' \code{nauf} models.  The main functionality which is not currently supported
#' is that the reference grid cannot be created specifying certain levels for
#' variables (i.e. the \code{at} argument; this is handled through the
#' \code{subset} argument to \code{nauf_pmmeans}).  A direct call to
#' \code{\link[lsmeans]{ref.grid}} will result in warnings (or possibly errors),
#' and inference made with the resulting object will be misleading and/or
#' incorrect.  Only \code{nauf_ref.grid} should be used. The \code{nauf.ref.grid}
#' returned by \code{nauf_ref.grid} can then be used as the \code{object}
#' argument to \code{nauf_pmmeans} to obtain predicted marginal means and
#' pairwise comparisons with p-values that adjust for familywise error rate.
#'
#' The \code{specs} and \code{pairwise} arguments to \code{nauf_pmmeans}
#' indicate what variables marginal means should be calculated for and wheter
#' pairwise comparisons of these means should be made.  If \code{specs} is a
#' character vector, then \code{pairwise} is used; if \code{specs} is a formula,
#' then the full iteraction of the terms on the right hand side of the formula
#' is considered, and the left hand side is used to indicate pairwise
#' comparisons.  For example (where \code{rg} is a \code{nauf.ref.grid}):
#'
#' \preformatted{
#' # all of these calculate pmm's for each combination of the factors f1 and f2
#' # but not pairwise comparisons
#' nauf_pmmeans(rg, c("f1", "f2"))
#' nauf_pmmeans(rg, ~ f1 + f2)
#' nauf_pmmeans(rg, ~ f1 * f2)
#' nauf_pmmeans(rg, ~ f1:f2)
#'
#' # all of these calculate the same pmm's, and additionally pairwise comparions
#' nauf_pmmeans(rg, c("f1", "f2"), pairwise = TRUE)
#' nauf_pmmeans(rg, pairwise ~ f1 + f2)
#' nauf_pmmeans(rg, pairwise ~ f1 * f2)
#' nauf_pmmeans(rg, pairwise ~ f1:f2)
#' }
#'
#' If \code{specs} indicates a single covariate, the effect of an increase of
#' \code{1} in the covariate is computed.  If \code{specs} indicates multiple
#' covariates, the effect of a simultaneous increase of \code{1} in all of the
#' covariates is computed.  If \code{specs} indicates a combination of factors
#' and covariate(s), the the effect of an increase of \code{1} for the
#' covariates is calcualted for each level of the full interaction of the
#' factors.
#'
#' If \code{by} is specified, then \code{pairwise} is forced to 
#' \code{TRUE}, and pairwise comparisons are performed within each level of the 
#' full interaction of the factors listed in \code{by}, rather than performing 
#' all possible pairwise comparisons.  For example, if there are two factors
#' \code{f1} with levels \code{A}, \code{B}, and \code{C}, and a factor
#' \code{f2} with levels \code{D} and \code{E}:
#'
#' \preformatted{
#' # this will produce six pmmeans (A:D, A:E, B:D, B:E, C:D, C:E) and
#' # all 15 pairwise comparisons
#' nauf_pmmeans(rg, c("f1", "f2"), pairwise = TRUE)
#'
#' # this would produce the same six pmmeans, but only three pairwise
#' # comparisons (A:D - A:E, B:D - B:E, C:D - C:E)
#' nauf_pmmeans(rg, c("f1", "f2"), by = "f1")
#'
#' # this would produce the same six pmmeans, but only six pairwise comparisons
#' # (A:D - B:D, A:D - C:D, B:D - C:D, A:E - B:E, A:E - C:E, B:E - C:E)
#' nauf_pmmeans(rg, c("f1", "f2"), by = "f2")
#' }
#'
#' The reference grid returned by \code{nauf_ref.grid} contains combinations of
#' factors which are not actually possible in the data set.  For example,
#' if factor \code{f1} has levels \code{A} and \code{B}, and factor \code{f2}
#' is \code{NA} when \code{f1 = A}, and takes values \code{C} and \code{D}
#' when \code{f1 = B}, the reference grid will still contain the combinations
#' \code{f1 = A, f2 = C}; \code{f1 = A, f2 = D}; and \code{f1 = B, f2 = NA},
#' even though these combinations are not possible.  This is because it is
#' impossible to know without the user's knowledge which combinations
#' make sense.  In many cases, this is inconsequential for the computation
#' of predicted marginal means, since the coding of unordered factors in
#' \code{nauf} regressions will average over the effects.  In cases where
#' these rows in the reference grid will cause invalid estimates and pairwise
#' comparisons, the \code{subset} argument can be used in the call to
#' \code{nauf_pmmeans} to ensure only the correct subsets are considered.
#' The default for the \code{subset} argument is \code{NULL}, indicating that
#' the the entire reference grid should be considered.  If not \code{NULL}, then
#' \code{subset} must be a list which defines the valid subsets as lists of
#' named character vectors, where the name of the character vector is an
#' unordered factor in the model, and the vector itself contains the levels
#' which define the subset (including \code{NA} in the case of factors which
#' have \code{NA} values; when \code{NA} is specified as a level, there should
#' be no quotes around it).  Any row in the reference grid which matches the
#' definition of at least one of the groups defined in \code{subset} is kept,
#' and all others are dropped.  So, continuing with the \code{f1} and \code{f2}
#' example, if \code{f2 = NA} corresponds to \code{f2 = D} in meaning, and is
#' coded as \code{NA} because all \code{f1 = A} observations are by necessity
#' \code{f2 = D}, then to analyze the effect of \code{f1}, we want to compare
#' the groups \code{f1 = A, f2 = NA} and \code{f1 = B, f2 = D}, which we could
#' do with the following call:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f1", subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = "D")))
#' }
#'
#' This would produce an estimate for \code{f1 = A} and \code{f1 = B}, but
#' conditioning on the subset where \code{f1} is truly contrastive based on
#' \code{f2}.  If, on the other hand, \code{f2 = NA} does not correspond in
#' interpretation to either \code{f2 = C} or \code{f2 = D}, but rather indicates
#' that \code{f2} is simply not meaningful when \code{f1 = A}, we would want to
#' average over the effect of \code{f2} within \code{f1 = B}, and compare this
#' result to \code{f1 = A, f2 = NA}, which we could do with the following call:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f1", subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D"))))
#' }
#'
#' In this case, the second sub-list in the \code{subset} list indicates that if
#' \code{f1 = B} and either \code{f2 = C} or \code{f2 = D}, then it belongs to
#' the second subset.  In this case, the \code{subset} argument is actually not
#' necessary, since for \code{f1 = A}, we want to \emph{not consider} the effect
#' \code{f2}, and for \code{f2 = B}, we want to \emph{average over all possible
#' levels} of \code{f2}, and these are actually the same thing computationally
#' for unordered factors in \code{nauf} models.  That is, we would get the same
#' result with:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f1")
#' }
#'
#' Generally speaking, if all of the factors in \code{specs} do \emph{not}
#' contain \code{NA} values, then the \code{subset} argument is unnecessary.
#' If any of the factors in \code{specs} \emph{do} contain \code{NA} values,
#' then you will almost always want to use the \code{subset} argument.  Now
#' consider that we are interested in \code{f2}.  Because \code{f2} is only
#' contrastive when \code{f1 = B}, we probably want to call:
#'
#' \preformatted{
#' # note that because there are not multiple subsets being specified, you
#' # don't have to specify subset = list(list(f1 = "B")); nauf_pmmeans will
#' # assume list(f1 = "B") means list(list(f1 = "B"))
#' nauf_pmmeans(rg, "f2", subset = list(f1 = "B"))
#' }
#'
#' This call will produce two estimates, one for \code{f2 = C} and one for
#' \code{f2 = D}, conditioning on \code{f1 = B}.  There will be no estimate for
#' \code{f2 = NA} because, by default, no estiamtes are produced for
#' combinations of factors where one factor is \code{NA}.  If we wanted to
#' compare the three possible groups (i.e. \code{f1 = A, f2 = NA};
#' \code{f1 = B, f2 = C}; and \code{f1 = B, f2 = D}), then we could additionally
#' use the \code{na_as_level} argument and change our \code{subset}:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f2", subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D"))),
#'   na_as_level = "f2")
#'
#' # this gives the same estimates, but the output will also show the
#' # corresponding level of f1, which is more transparent
#' nauf_pmmeans(rg, c("f1", "f2"), subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D"))),
#'   na_as_level = "f2")
#' }
#'
#' The easiest way to use the \code{subset} argument is to create a list that
#' defines valid subsets for different regression terms of interest outside of
#' \code{nauf_pmmeans}, and then using the relevant element of the list in the
#' \code{nauf_pmmeans} call.  For example:
#'
#' \preformatted{
#' pmmsubs <- list()
#' pmmsubs$f1 <- list(list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D")))
#' pmmsubs$f2 <- list(f1 = "B")
#'
#' nauf_pmmeans(rg, "f1", subset = pmmsubs$f1)
#'
#' nauf_pmmeans(rg, "f2", subset = pmmsubs$f2)
#'
#' nauf_pmmeans(rg, c("f1", "f2"), subset = pmmsubs$f1, na_as_level = "f2")
#' }
#'
#' This way you can just define the different subsets of the data once and not
#' have to think about it at every \code{nauf_pmmeans} call.
#'
#' @param mod A regression model fit with \code{nauf} contrasts.
#' @param KR Only applies when \code{mod} is a \code{\linkS4class{nauf.lmerMod}}
#'   fit with \code{REML = TRUE}. If \code{KR = TRUE}, then the Kenward-Roger
#'   approximation is used to calculate degrees of freedom.  If
#'   \code{KR = FALSE} (the default), then the Satterthwaite approximation is
#'   used.  When \code{mod} is a
#'   \code{\linkS4class{nauf.lmerMod}} fit with \code{REML = FALSE}, then the
#'   Satterthwaite approximation is always used.  The Kenward-Roger method
#'   is implemented with \code{\link[pbkrtest]{Lb_ddf}} and the Satterthwaite
#'   method is implemented with \code{\link[lmerTest]{calcSatterth}}.
#' @param object A \code{nauf.ref.grid} object created with
#'   \code{nauf_ref.grid}.
#' @param specs The fixed effects for which the full interaction term should be
#'   considered in the calculation of predicted marginal means. The preferred
#'   method is to specify the variables as a character vector.  However, they
#'   can also be specified on the right hand side of a formula, optionally with
#'   the keyword \code{pairwise} on the left hand side to indicate that pairwise
#'   comparions should be performed.
#' @param pairwise A logical (default \code{FALSE}) indicating whether pairwise
#'   comparisons of the predicted marginal means should be performed.  If
#'   \code{specs} is a formula, then the \code{pairwise} argument is ignored and
#'   the left hand side of the formula is used to determined whether pairwise
#'   comparsions should be made.  If \code{by} is not \code{NULL} then
#'   \code{pairwise} is forced to \code{TRUE}.
#' @param subset A list indicating which subsets of the reference grid should be
#'   considered in the calculation of the predicted marginal means. See
#'   'Details'.
#' @param na_as_level A character vector of unordered factors in \code{specs}
#'   that have \code{NA} values that should be considered as levels. The default
#'   \code{NULL} indicates that \code{NA} should not be considered as a level
#'   for any unordered factors in \code{specs}. See 'Details'.
#' @param by An optional character vector specifying unordered factors in
#'   \code{specs}.  If specified, then pairwise comparisons are performed
#'   within each level of the full interaction term of the factors, rather than
#'   for all possible combinations.  If an unordered factor listed in \code{by}
#'   is not included in \code{specs}, it is added to \code{specs} automatically.
#' @param ... Additional arguments are ignored with a warning.
#'
#' @return
#' \code{nauf_ref.grid} returns a \code{nauf.ref.grid} object, which is just a
#' list with one element \code{ref.grid} of class
#' \code{\link[lsmeans]{ref.grid-class}}.  This reference grid should not be
#' used directly with \code{\link[lsmeans]{lsmeans}}, but rather only with
#' \code{nauf_pmmeans}. \code{nauf_pmmeans} returns a 
#' \code{\link{nauf.pmm.list}} object.
#'
#' @seealso \code{\link{nauf_contrasts}}, \code{\link{nauf_glm}},
#'   \code{\link{nauf_glmer}}, \code{\link{nauf_stan_glm}}, and
#'   \code{\link{nauf_stan_glmer}}.
#'
#' @name nauf-pmmeans
NULL


#' List of predicted marginal means objects for \code{nauf} models.
#'
#' The \code{\link{nauf_pmmeans}} function returns an object of class
#' \code{nauf.pmm.list}.
#'
#' The \code{nauf.pmm.list} object contains a first element \code{pmmeans}
#' which contains the predicted marginal means, and possibly additional
#' \code{contrasts} elements depending on whether the call to 
#' \code{\link{nauf_pmmeans}} indicated that pairwise comparisons should be made.  
#' If the \code{by} argument was specified, there will be a numbered 
#' \code{contrasts} element for each level of the full interaction of the 
#' factors listed in the \code{by} argument.  The object also has a \code{specs}
#' attribute which contains information about the variables and subsets from the 
#' call to the function.
#'
#' The object has \code{summary} and \code{print} methods which first print 
#' information from the \code{specs} attribute, and then pass addtional function
#' arguments along to each element in the list.  For frequentist regressions, 
#' the elements in the \code{nauf.pmm.list} are 
#' \code{\link[lsmeans]{lsmobj-class}} objects. In this case, the \code{summary} 
#' and \code{print} methods for the \code{nauf.pmm.list} take arguments listed in 
#' \code{\link[lsmeans]{summary.ref.grid}} (e.g. \code{infer}, \code{type}, 
#' \code{adjust}, etc.). For Bayesian regressions, the elements in the 
#' \code{nauf.pmm.list} are \code{\link{nauf.pmm.stan}} objects. In this case, 
#' the \code{summary} and \code{print} methods for the \code{nauf.pmm.list} take 
#' arguments which are listed in \code{\link{nauf.pmm.stan}}.
#'
#' For posterior marginal means from a Bayesian \code{nauf} model, there is
#' a method for \code{\link[=as.shinystan,nauf.pmm.list-method]{as.shinystan}} 
#' which allows \code{\link[shinystan]{launch_shinystan}} to be used.
#'
#' @name nauf.pmm.list
setOldClass("nauf.pmm.list")


#' Posterior samples of marginal means from Bayesian \code{nauf} models.
#'
#' When \code{\link{nauf_pmmeans}} is used with Bayesian regressions, the 
#' elements of the resulting \code{\link{nauf.pmm.list}} have class 
#' \code{nauf.pmm.stan}.
#'
#' The \code{nauf.pmm.stan} object is a list with the following elements.
#'
#' \describe{
#'   \item{names}{A data frame with the levels of each factor (or 'inc_1' for
#'     covariates), with one row for each element in the third dimension of 
#'     \code{samples}.}
#'   \item{contrasts}{The fixed effects model matrix showing the contrasts 
#'     applied to the regression coefficients (rows correspond to the 
#'     \code{names} element).}
#'   \item{samples}{An array with three dimensions. The first corresponds to 
#'     iterations, the second to chains, and the third to parameters (the same 
#'     structure returned by \code{\link[rstan]{as.array.stanfit}}).}
#'   \item{family}{The regression family.}
#'   \item{inv.lbl}{The label of the inverse link (e.g. probability, rate, etc.).}
#'   \item{misc}{A list with additional information based on the type of model.}
#' }
#'
#' @param object A \code{nauf.pmm.stan} object.
#' @param probs A vector of quantiles to calculate for the estimates.
#'   The default is a 95% credible interval (also called an uncertainty interval). 
#'   The median (0.5) is always added regardless of whether it is specified.
#' @param type If \code{"link"} (the default), then the estimates are not 
#'   transformed prior to summary; if \code{"response"}, then the inverse link 
#'   function in the object's \code{family} element is applied to the 
#'   \code{samples} element prior to summary.
#' @param x Either a \code{nauf.pmm.stan} object, or the \code{summ.nauf.pmm.stan} 
#'   object returned by calling \code{summary} on a \code{nauf.pmm.stan} object.
#' @param row.names,optional Changes from the defaults are ignored.
#' @param mcse A logical indicating whether the Monte Carlo standard errors
#'   should be printed (default \code{FALSE} since the column is always recoverable
#'   by dividing the \code{SD} column by the square root of the \code{ESS} column).
#' @param rhat An optional logical indicating whether or not to print the 
#'   Gelman-Rubin R-hat statistic.  If \code{rhat = NULL} (the default), then 
#'   the \code{Rhat} column of the summary is only printed if any of the 
#'   statistics are greater that \code{1.1}. Regardless of the value of the 
#'   \code{rhat} argument, a warning is issued if any of the statistics are 
#'   greater than \code{1.1}. 
#' @param ... See the \code{Methods} section.
#'
#' @return The returned object depends on the function.  See the 'Methods'
#'   section.
#'
#' @section Methods:
#' \describe{
#'   \item{summary}{Returns a data frame with class
#'     \code{summ.nauf.pmm.stan}, with means, standard deviations, medians, mean 
#'     absolute differences, the posterior probability that the estimate is 
#'     greater than zero, quantiles specified in \code{probs}, the effective 
#'     sample size, Monte Carlo standard error, and Gelman-Rubin R-hat statistic.
#'     Additional arguments in \code{...} are ignored.  Quantiles, effective
#'     samples size, and Monte Carlo standard error, and R-hat statistics are
#'     computed using \code{\link[rstan]{monitor}}.}
#'   \item{print}{For \code{summ.nauf.pmm.stan} objects, the data frame is
#'     printed, omitting the Monte Carlo standard error and 
#'     R-hat statistic based on \code{mcse} and \code{rhat}.  If the inverse
#'     link function was used, then a message indicating the response type
#'     is also printed. Additional arguments in \code{...} are passed to 
#'     \code{\link[base]{print.data.frame}}. For \code{nauf.pmm.stan} objects,
#'     first \code{summary} is called, passing along \code{...},
#'     and then print is called on the resulting \code{summ.nauf.pmm.stan}, also
#'     passing along \code{...} arguments.}
#'   \item{as.data.frame}{For \code{summ.nauf.pmm.stan} objects, removes the 
#'     \code{data.frame} contained in a \code{summ.nauf.pmm.stan}.  Additional
#'     arguments in \code{...} are ignored. For \code{nauf.pmm.stan} objects,
#'     first \code{summary} is called, passing along \code{...}, and then
#'     \code{as.data.frame} is called on the resulting \code{summ.nauf.pmm.stan}
#'     object, ignoring \code{...} arugments.}
#'   \item{as.array}{Returns the \code{samples} element of a 
#'     \code{nauf.pmm.stan} object. Additional arguments in \code{...} are 
#'     ignored.}
#'   \item{as.matrix}{Returns the \code{samples} element of a 
#'     \code{namf.pmm.stan} object, flattening the array into a matrix such that 
#'     the rows represent iterations ordered by chain and the columns represent 
#'     parameters.}
#' }
#'
#' @name nauf.pmm.stan
setOldClass("nauf.pmm.stan")


###### nauf.ref.grid ######

#' @export
terms.nauf.ref.grid <- function(x, ...) {
  return(x$ref.grid@model.info$terms)
}


#' @export
print.nauf.ref.grid <- function(x, ...) {
  show(x$ref.grid)
}



###### nauf.pmm.stan ######

#' @rdname nauf.pmm.stan
#'
#' @method summary nauf.pmm.stan
#'
#' @importFrom rstan monitor
#' 
#' @export
summary.nauf.pmm.stan <- function(object, probs = c(0.025, 0.975),
                                  type = c("link", "response"), ...) {
  type <- match.arg(type)
  lbl <- "link"
  
  if (type == "response") {
    if (!is.null(object$family) && is.function(inv <- object$family$linkinv)) {
      object$samples <- inv(object$samples)
      if (is.null(lbl <- object$inv.lbl)) lbl <- "response"
    } else {
      warning("'type' is 'response' but there is no inverse link function.",
        "  No transformation performed.")
    }
  }
  
  probs <- sort(union(probs, 0.5))
  summ <- as.data.frame(rstan::monitor(object$samples, warmup = 0,
    probs = probs, print = FALSE))
  
  summ$MAD <- apply(object$samples, 3, mad)
  summ[["P(>0)"]] <- apply(object$samples > 0, 3, mean)
  nc <- ncol(summ)
  med <- match("50%", colnames(summ))
  
  # mean se_mean sd ... 50%    ... n_eff   Rhat   MAD     P(>0)
  # 1    2       3      med        nc-3    nc-2   nc-1    nc
  # Mean MCSE    SD ... Median ... ESS     Rhat   MAD     P(>0)
  colnames(summ)[c(1:3, med, nc - 3)] <- c("Mean", "MCSE", "SD", "Median", "ESS")
  
  # Mean   SD   Median MAD   P(>0)  ...        ESS    MCSE   Rhat
  # 1      3    med    nc-1  nc     grep("%")  nc-3   2      nc-2
  ord <- c(1, 3, med, nc - 1, nc, grep("%", colnames(summ)), nc - 3, 2, nc - 2)
  summ <- summ[ord]
  
  if (is.data.frame(object$names)) {
    summ <- cbind(object$names, summ)
  } else if (is.character(object$names)) {
    summ <- cbind(data.frame(Name = object$names), summ)
  } else {
    summ <- cbind(data.frame(Name = dimnames(object$samples)[[3]]), summ)
  }
  rownames(summ) <- NULL
  
  return(structure(summ, class = "summ.nauf.pmm.stan", type = lbl,
    misc = object$misc))
}


#' @rdname nauf.pmm.stan
#'
#' @export
print.summ.nauf.pmm.stan <- function(x, row.names = FALSE, mcse = FALSE,
                                     rhat = NULL, ...) {
  type <- attr(x, "type")
  attr(x, "type") <- NULL
  misc <- attr(x, "misc")
  attr(x, "misc") <- NULL
  class(x) <- "data.frame"
  
  if (!mcse) x <- x[-which(colnames(x) == "MCSE")]
  if (w <- any(x$Rhat > 1.1)) {
    warning("Some Rhat's are greater than 1.1; chains have not converged.")
  }
  if (isFALSE(rhat) || (is.null(rhat) && !w)) {
    x <- x[-which(colnames(x) == "Rhat")]
  }
  
  print(x, row.names = FALSE, ...)

  if (type != "link") {
    cat("\nEstimates are given on the", type, "scale.\n")
  }
}


#' @rdname nauf.pmm.stan
#'
#' @export
as.array.nauf.pmm.stan <- function(x, ...) {
  return(x$samples)
}


#' @rdname nauf.pmm.stan
#'
#' @export
as.matrix.nauf.pmm.stan <- function(x, ...) {
  m <- apply(x$samples, 3, function(y) y)
  if (!is.matrix(m)) {
    m <- t(m)
  }
  dimnames(m) <- list(iterations = NULL, parameters = dimnames(x$samples)[[3]])
  return(m)
}


#' @rdname nauf.pmm.stan
#'
#' @export
print.nauf.pmm.stan <- function(x, ...) {
  print(summary(x, ...), ...)
}


#' @rdname nauf.pmm.stan
#'
#' @export
as.data.frame.nauf.pmm.stan <- function(x, row.names = NULL, optional = FALSE,
                                        ...) {
  return(as.data.frame(summary(x, ...)))
}


#' @rdname nauf.pmm.stan
#'
#' @export
as.data.frame.summ.nauf.pmm.stan <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  class(x) <- "data.frame"
  attr(x, "misc") <- NULL
  attr(x, "type") <- NULL
  return(x)
}


###### nauf.pmm.list ######

#' Create a shinystan object for posterior marginal means from nauf models.
#'
#' Joins the \code{samples} elements from the \code{\link{nauf.pmm.stan}} objects
#' in a \code{\link{nauf.pmm.list}}, adds the \code{log-posterior} to them,
#' and passes the resulting array to the \code{array} method for
#' \code{\link[shinystan]{as.shinystan}}.
#'
#' @param X A \code{\link{nauf.pmm.list}}.
#' @param ... Not used.
#'
#' @return A \code{\link[shinystan]{shinystan}} object based on the posterior
#' marginal means (and possibly pairwise comparisons) in \code{X}.  This allows
#' the results to be viewed in the \code{shinystan} GUI as if they had been
#' included in the generated quantities block of the model.
#'
#' @seealso \code{\link[shinystan]{launch_shinystan}}
#'
#' @export
setMethod("as.shinystan", signature(X = "nauf.pmm.list"), function(X, ...) {
  if (!all(sapply(X, is.nauf.pmm.stan))) {
    stop("'X' is not from a Bayesian model")
  }
  
  a <- X[[1]]$misc$shiny
  samples <- cbind(do.call(cbind, lapply(X, as.matrix)), a$lp)
  samples <- make_stan_mcmc_array(samples, nchain(X))
  
  out <- shinystan::as.shinystan(samples,
    model_name = "Posterior Marginal Means (nauf)",
    sampler_params = a$params, algorithm = "NUTS", max_treedepth = a$mtd)
})


#' @method summary nauf.pmm.list
#'
#' @export
summary.nauf.pmm.list <- function(object, ...) {
  specs <- attr(object, "specs")
  attr(object, "specs") <- NULL
  
  if (!specs$bayes) {
    class(object) <- c("lsm.list", "list")
    summ <- summary(object, ...)
  } else {
    class(object) <- "list"
    summ <- mlapply(object = object, same = list(...), fun = summary)
  }
  
  first_class(summ) <- "summ.nauf.pmm.list"
  attr(summ, "specs") <- specs
  return(summ)
}


#' @export
print.summ.nauf.pmm.list <- function(x, ...) {
  specs <- attr(x, "specs")
  attr(x, "specs") <- NULL
  drop_class(x) <- "summ.nauf.pmm.list"
  do.call(print_pmm_list, c(list(x = x), specs, list(...)))
}


print_pmm_list <- function(x, variables, pairwise, averaged_over, held_at_mean,
                           conditioned_on, keep_NA, drop_NA, subset, by, bayes,
                           note, ...) {
  cat("\n", if (bayes) "Posterior" else "Predicted", " marginal means for '",
    paste(variables, collapse = ":"), "'", sep = "")
  
  if (length(keep_NA)) {
    cat("\nNA considered a level for:", add_quotes(keep_NA))
  }
  if (length(drop_NA)) {
    cat("\nNA not considered a level for:", add_quotes(drop_NA))
  }
  if (length(by)) {
    cat("\nPairwise comparisons by '", paste(by, collapse = ":"), "'", sep = "")
  }
  if (length(note)) {
    cat("\nNote:", note)
  }
  cat("\n")

  if (length(averaged_over)) {
    cat("\nFactors averaged over:", add_quotes(averaged_over))
  }
  if (length(held_at_mean)) {
    cat("\nCovariates held at their means:", add_quotes(held_at_mean))
  }
  if (length(conditioned_on)) {
    cat("\nFactors conditioned on:", add_quotes(conditioned_on),
      "\n\nSee the 'subset' element of the 'specs' attribute for subsetted",
      "groups")
  }

  cat("\n")
  if (sum(lengths(list(averaged_over, held_at_mean, conditioned_on)))) {
    cat("\n")
  }

  if (bayes) {
    for (i in 1:length(x)) {
      cat("$", names(x)[i], "\n", sep = "")
      print(x[[i]], ...)
      cat("\n")
    }
  } else {
    drop_class(x) <- "summ.nauf.pmm.list"
    print(x, ...)
  }
}


#' @export
print.nauf.pmm.list <- function(x, ...) {
  print(summary(x, ...), ...)
}

