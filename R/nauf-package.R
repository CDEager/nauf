

#' nauf: Regression with Not Applicable Unordered Factors in R.
#'
#' It is often the case that a factor only makes sense in a subset of a dataset
#' (i.e. for some observations a factor may simply not be meaningful), or that
#' with observational datasets there are no observations in some levels of an
#' interaction term.  The \code{nauf} package allows unordered factors to be
#' coded as \code{NA} in the subsets of the data where they are not applicable
#' or otherwise not contrastive.  Sum contrasts are used for all unordered
#' factors, setting \code{NA} values to zero.  This allows all of the data to be
#' modeled together without creating collinearity or making the output difficult
#' to interpret.  It is highly recommended that regression variables be put
#' on the same scale with \code{\link[standardize]{standardize}} prior to
#' using \code{nauf} functions (though this is not required for the functions
#' to work).
#'
#' @section Contrasts:
#' A detailed description of how \code{NA} values are treated is given in 
#' \code{\link{nauf_contrasts}}. These contrasts are implemented automatically 
#' through \code{\link{nauf_model.frame}}, which stores the information required 
#' to make fixed effects and random effects model matrices with \code{nauf} 
#' contrasts in its \code{terms} attribute.  For details on the \code{terms} 
#' attribute, see \code{\linkS4class{nauf.terms}}.  For fixed effects and random 
#' effects model matrices, see \code{\link{nauf_model.matrix}} and 
#' \code{\link{nauf_glFormula}}.
#'
#' @section Fixed Effects Regressions:
#' There are currently three fixed effects model fitting functions for which 
#' \code{nauf} contrasts have been implemented.  The \code{\link{nauf_lm}} 
#' function fits linear models (based on \code{\link[stats]{lm}} in the 
#' \code{stats} package).  The \code{\link{nauf_glm}} function fits generalized 
#' linear models other than negative binomial models (e.g. \code{binomial}, 
#' \code{poisson}, etc.; based on \code{\link[stats]{glm}} in the \code{stats} 
#' package).  And the \code{\link{nauf_glm.nb}} function fits negative binomial 
#' models (based on \code{\link[MASS]{glm.nb}} in the \code{MASS} package).  
#' These regression functions return objects with the same elements as the 
#' functions they are based on, and have class \code{\linkS4class{nauf.glm}}.  
#' Standard generic methods defined for classes \code{negbin}, \code{glm}, and 
#' \code{lm} work as they normally would.
#' 
#' @section Mixed Effects Regressions:
#' There are currently three mixed effects model fitting functions for which 
#' \code{nauf} contrasts have been implemented.  They are 
#' \code{\link{nauf_lmer}}, \code{\link{nauf_glmer}}, and 
#' \code{\link{nauf_glmer.nb}} (based on \code{\link[lme4]{lmer}}, 
#' \code{\link[lme4]{glmer}}, and \code{\link[lme4]{glmer.nb}} in the 
#' \code{lme4} package), and fit the same types of regressions as 
#' \code{\link{nauf_lm}}, \code{\link{nauf_glm}}, and \code{\link{nauf_glm.nb}}, 
#' but with mixed effects rather than just fixed effects.  
#' These regression functions return objects of class 
#' \code{\linkS4class{nauf.lmerMod}} (for \code{nauf_lmer}) and 
#' \code{\linkS4class{nauf.glmerMod}} (for \code{nauf_glmer} and 
#' \code{nauf_glmer.nb}), which inherit from \code{\linkS4class{lmerMod}} and 
#' \code{\linkS4class{glmerMod}}, respectively.  Most generic methods defined for 
#' these classes function as they normally would, with the notable exception 
#' of \code{\link[=predict.nauf.merMod]{predict}}, where the \code{nauf.merMod} 
#' method does not incoporate the full range of functionality of 
#' \code{\link[lme4]{predict.merMod}}.
#'
#' @section p-values:
#' For \code{\linkS4class{nauf.glm}} models, the standard 
#' \code{\link[stats]{anova}} function can be used to obtain p-values for each 
#' predictor.  For \code{\linkS4class{nauf.lmerMod}} and 
#' \code{\linkS4class{nauf.glmerMod}} models, there is an 
#' \code{\link[=anova.nauf.merMod]{anova method}} that performs Type III tests 
#' with various options for p-value calculation, such as likelihood ratio tests, 
#' parameteric bootstrapping, and the Kenward-Roger and Satterthwaite 
#' approximations for F-test denominator degrees of freedom.
#'
#' For both fixed effects models and mixed effects models fit with \code{nauf} 
#' contrasts, the \code{\link{nauf_ref.grid}} and \code{\link{nauf_pmmeans}} 
#' functions can be used to calculate predicted marginal means (also often 
#' called least-squares means or lsmeans), and further allow these estimates to 
#' be generated over only the subset of the reference grid where the factors in 
#' question are contrastive.
#'
#' @docType package
#' @name nauf-package
#' @aliases nauf
#'
#' @import stats
#' @import MASS
#' @import methods
#' @import lme4
#' @import standardize
#' @import stringr
#' @importClassesFrom lme4 merMod lmerMod glmerMod
NULL

