

#' nauf: Add new title.
#'
#' Add Description
#'
#' Add Details
#'
#' Documentation plan. Explain all of the nauf contrast methods in the doc for
#' nauf_contrasts and do not export nauf_interaction. Explain the basic fixef
#' functionality in the documentation for nauf.glm class.  For the nauf.glm
#' functions, link to the nauf.glm class doc and to the functions they draw
#' their code from.  Document nauf.lmerMod and nauf.glmerMod in one doc and
#' reference same treatment as for fixef in nauf.glm, and add note about
#' grouping factor NAs.  Documentation pages should be nauf-package,
#' nauf-contrasts, nauf.glm-class, nauf.merMod-class, and a page for each
#' regression fitting function.  Modular nauf lme4 functions should not be
#' exported.  For generics, minimal documentation except for predict for
#' nauf.merMod's since the functionality is restricted.
#'
#' @section TODO:
#' \describe{
#'   \item{lme4}{Finish implementation of lme4 functions with nauf methods and
#'     add unit tests.}
#'   \item{afex}{Implement nauf methods at least for the mixed function and
#'     possibly some others.}
#'   \item{lsmeans}{Implement nauf methods for ref.grid and pmmeans.  Only
#'     relatively basic functionality will be supported.}
#'   \item{documentation}{Update all documentation.}
#'   \item{rstanarm}{Initial release on CRAN without rstanarm compatibility?}
#' }
#'
#' @docType package
#' @name nauf-package
#'
#' @import stats
#' @import MASS
#' @import methods
#' @import lme4
#' @import Matrix
#' @import standardize
NULL

