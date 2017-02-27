
#' nauf: An R package for fitting regressions with not applicable values in unordered factors.
#'
#' The \code{nauf} package provides tools for fitting regressions where
#' unordered factors may have \code{NA} values in addition to their
#' regular levels.  In this package, \code{NA} is used to encode when an
#' unordered factor is truly \emph{not applicable}.  This is different
#' than "not available" or "missing at random".  The concept applies only to
#' unordered factors, and indicates that the factor is simply not meaningful
#' for an observation, or that while the observation may technically be
#' definable by one of the factor levels, the interpretation of its belonging
#' to that group isn't the same.  For imbalanced observational data, coding
#' as \code{NA} may also be used to control for a factor which is only
#' contrastive within a subset of the data due to the sampling scheme.  The
#' \code{nauf} package provides functions to implement these values
#' automatically; the user simply needs to code the relevant observations as
#' \code{NA}.
#'
#' The \code{nauf} package works by implementing new methods for the generics
#' \code{model.frame} and \code{model.matrix} from the \code{stats} package,
#' ensuring that these generics are used by already existing regression fitting
#' functions. All unordered factors are coded with sum contrasts using
#' \code{\link{named_contr_sum}}, which uses the contrasts returned by
#' \code{\link[stats]{contr.sum}}, but names the dummy variables rather than
#' numbering them so the output is more easily interpreted.  These contrasts
#' are assigned ignoring \code{NA} values, and then in the model matrix
#' all \code{NA} values are set to zero.  Interaction terms
#' which involve more than one unordered factor, with at least one of these
#' unordered factors having \code{NA} values, may require different contrasts
#' to be used in the interaction term than in the main effects, as determined
#' with \code{\link{nauf_interaction}}.  A list of the contrast matrices used
#' in a regression can be obtained with \code{\link{nauf_contrasts}}.
#'
#' The function \code{\link{nauf_reg}}, is a wrapper function that assigns
#' the formula for a regression an extra (first) class attribute of
#' \code{\link[=nauf-class]{nauf}}, and then calls an appropriate regression
#' fitting function from the \code{stats} or \code{MASS} packages, causing these
#' functions to use the \code{nauf} generics for \code{model.frame} and
#' \code{model.matrix} (see \code{\link{nauf_model_frame}} and
#' \code{\link{nauf_model_matrix}}).  Currently, any regression which can be
#' fit with the functions \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#' and \code{\link[MASS]{glm.nb}} can be fit using \code{\link{nauf_reg}} to
#' allow \code{NA} values in the unordered factors.  The functions
#' \code{\link{nauf_grid}} and \code{\link{nauf_pmmeans}} interface with
#' \code{\link[lsmeans]{ref.grid}} and \code{\link[lsmeans]{lsmeans}} to
#' obtain predicted marginal means (also called least-squares means) and
#' pairwise comparisons for factors, interactions between factors, and
#' interactions between factors and covariates, allowing for the estimates
#' to be generated for only the relevant subsets of the data.
#'
#' @docType package
#' @name nauf-package
#'
#' @import stats
#' @import MASS
#' @import methods
#' @import lsmeans
NULL

