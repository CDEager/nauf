

#' Regression with NA values in unordered factors.
#'
#' It is often the case that a factor only makes sense in a subset of a dataset
#' (i.e. for some observations a factor may simply not be meaningful), or that
#' with observational datasets there are no observations in some levels of an
#' interaction term.  There are also cases where a random effects grouping
#' factor is only applicable in a subset the data, and it is desireable to model
#' the noise introduced by the repeated measures on the group members within the
#' subset of the data where the repeated measures exist. The \code{nauf} package
#' allows unordered factors and random effects grouping factors to be coded as
#' \code{NA} in the subsets of the data where they are not applicable or
#' otherwise not contrastive.  It is highly recommended that variables be put on
#' the same scale with \code{\link[standardize]{standardize}} prior to using
#' \code{nauf} functions (though this is not required).
#'
#' @section Contrasts:
#' A detailed description of how \code{NA} values are treated is given in
#' \code{\link{nauf_contrasts}}. These contrasts are implemented automatically
#' through \code{\link{nauf_model.frame}}, which stores the information required
#' to make fixed effects and random effects model matrices with \code{nauf}
#' contrasts in its \code{terms} attribute.  For details on the \code{terms}
#' attribute, see \code{\link{nauf.terms}}.  For fixed effects and random
#' effects model matrices, see \code{\link{nauf_model.matrix}} and
#' \code{\link{nauf_glFormula}}.
#'
#' @section Regressions:
#' \code{nauf} contrasts have been implemented for fixed effects regressions
#' that would normally be fit with \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, and \code{\link[MASS]{glm.nb}}; and for mixed
#' effects regressions that would normally be fit with \code{\link[lme4]{lmer}},
#' \code{\link[lme4]{glmer}}, and \code{\link[lme4]{glmer.nb}}.  For fixed
#' effects \code{nauf} regressions, see \code{\link{nauf_glm}}, and for mixed
#' effects \code{nauf} regressions, see \code{\link{nauf_glmer}}.
#'
#' @section ANOVAs:
#' For fixed effects \code{nauf} models, the \code{anova} function uses the
#' methods for corresponding non-\code{nauf} models (i.e.
#' \code{\link[stats]{anova.lm}}, \code{\link[stats]{anova.glm}}, or
#' \code{\link[MASS]{anova.negbin}} depending on the regression family).  For
#' mixed effects \code{nauf} models, the \code{anova} function uses the
#' \code{\link{anova.nauf.merMod}} method, which with default arguments is the
#' same as \code{\link[lme4]{anova.merMod}}.  The \code{\link{anova.nauf.merMod}}
#' method also allows Type III tests to be made via likelihood ratio tests,
#' parametric bootstrapping, and, for linear models, the Satterthwaite and
#' Kenward-Roger approximations of denominator degrees of freedom (similar to
#' the \code{\link[afex]{mixed}} function in the \code{afex} package).
#'
#' @section Predicted Marginal Means:
#' The \code{\link{nauf_ref.grid}} function can be used to construct reference
#' grids for \code{nauf} models (constructing a
#' \code{\link[lsmeans]{ref.grid-class}} object).  Predicted marginal means
#' (often called least-squares means) can be calculated with this reference grid
#' using the \code{\link{nauf_pmmeans}} function.  The function also allows the
#' user to flexibly specify a subset of the reference grid to use when
#' calculating the marginal means, so that the effect of a factor can be tested
#' with regards to the subset of the data where it is contrastive.
#'
#' @section Datasets and Vignette:
#' For detailed examples of how to use the \code{nauf} package, see the
#' "Using the nauf package" vignette, which makes use of the
#' \code{\link{plosives}} and \code{\link{fricatives}} datasets included in the
#' package.
#'
#' @import stats
#' @import MASS
#' @import methods
#' @import lme4
#' @import standardize
#' @import stringr
#' @importClassesFrom lme4 merMod lmerMod glmerMod
"_PACKAGE"

