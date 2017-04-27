

#' Class for \code{nauf} regression formulas.
#'
#' The \code{formula} attribute of a \code{\linkS4class{nauf.frame}} has the
#' (S3) class \code{nauf.formula} to ensure that calls to
#' \code{\link[stats]{model.frame}} which use the formula call
#' \code{\link{nauf_model.frame}}.  There are currently no generic methods for
#' \code{nauf.formula} objects meant for end user use.
#'
#' @seealso \code{\link{nauf_contrasts}}, \code{\link{nauf_model.frame}}, and
#'   \code{\link{is.nauf.formula}}.
#'
#' @aliases nauf.formula
#'
#' @name nauf.formula-class
NULL


#' Class for \code{terms} objects which contain information about \code{nauf} contrasts.
#'
#' When \code{\link{nauf_model.frame}} is called, a
#' \code{\linkS4class{nauf.frame}} is returned, and this object's \code{terms}
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
#' There are currently no generic methods for \code{nauf.terms} objects meant
#' for end user use.
#'
#' @seealso \code{\link{nauf_contrasts}}, \code{\link{nauf_model.frame}}, and
#'   \code{\link{is.nauf.terms}}.
#'
#' @aliases nauf.terms
#'
#' @name nauf.terms-class
NULL


#' Class for model frames created with \code{nauf_model.frame}.
#'
#' A model frame created by \code{\link{nauf_model.frame}}.  It has a \code{formula}
#' attribute, which is a \code{\linkS4class{nauf.formula}}, and a \code{terms}
#' attribute, which is a \code{\linkS4class{nauf.terms}}.  The (S3) class
#' \code{nauf.frame} comes \emph{after} the \code{data.frame} class attribute so
#' that it can be included as the \code{frame} slot in a
#' \code{\linkS4class{nauf.lmerMod}} or \code{\linkS4class{nauf.glmerMod}}.
#' There are currently no generic methods for \code{nauf.frame} objects meant
#' for end user use.
#'
#' @seealso \code{\link{nauf_contrasts}}, \code{\link{nauf_model.frame}}, and
#'   \code{\link{is.nauf.frame}}.
#'
#' @aliases nauf.frame
#'
#' @name nauf.frame-class
NULL


#' Class of fitted fixed effects models with \code{nauf} contrasts.
#'
#' A fitted fixed effects regression model returned by \code{\link{nauf_lm}},
#' \code{\link{nauf_glm}}, or \code{\link{nauf_glm.nb}} with an additional first
#' (S3) class attribute \code{nauf.glm}.  There are currently no generic methods
#' for class \code{nauf.glm}, and generics from the \code{base}, \code{stats},
#' and \code{MASS} packages such as \code{summary}, \code{predict}, etc. work
#' just as they would for models fit with \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, or \code{\link[MASS]{glm.nb}}.  If you encounter a
#' generic function in these packages which does not function properly, please
#' report the issue at \url{https://github.com/CDEager/nauf/issues}.
#'
#' @seealso \code{\link{nauf_lm}}, \code{\link{nauf_glm}},
#'   \code{\link{nauf_glm.nb}}, and \code{\link{nauf_contrasts}}
#'
#' @aliases nauf.glm
#'
#' @name nauf.glm-class
NULL


#' Class for fitted linear mixed effects models with \code{nauf} contrasts.
#'
#' A fitted mixed effects regression model returned by \code{\link{nauf_lmer}}
#' of S4 class \code{nauf.lmerMod}, which inherits from
#' \code{\linkS4class{lmerMod}}.
#'
#' There are S3 methods for
#' the \code{\link[=predict.nauf.merMod]{predict}} and 
#' \code{\link[=anova.nauf.merMod]{anova}} functions meant for end user use.
#' There is also an S3 method for \code{\link[stats]{simulate}}, but with more
#' limited functionality than \code{\link[lme4]{simulate.merMod}}; the
#' \code{nauf.lmerMod} method for this function is not meant for end users, and
#' is necessary for the \code{PB} method to work for
#' \code{\link[=anova.nauf.merMod]{anova}}.
#' Other generics from the \code{base}, \code{stats}, and
#' \code{lme4} packages such as \code{summary} etc. work just as
#' they would for models fit with \code{\link[lme4]{lmer}}.  If you encounter a
#' generic function in these packages which does not function properly, please
#' report the issue at \url{https://github.com/CDEager/nauf/issues}.
#'
#' @seealso \code{\link{nauf_lmer}} and \code{\link{nauf_contrasts}}.
#'
#' @slot rsp,Gp,call,frame,flist,cnms,lower,theta,beta,u,devcomp,pp,optinfo See
#'   \code{\linkS4class{merMod}}.
nauf.lmerMod <- setClass("nauf.lmerMod", contains = "lmerMod")


#' Class for fitted generalized linear mixed effects models with \code{nauf} contrasts.
#'
#' A fitted mixed effects regression model returned by \code{\link{nauf_glmer}}
#' or \code{\link{nauf_glmer.nb}} of S4 class \code{nauf.glmerMod}, which
#' inherits from \code{\linkS4class{glmerMod}}.
#'
#' There are S3 methods for
#' the \code{\link[=predict.nauf.merMod]{predict}} and 
#' \code{\link[=anova.nauf.merMod]{anova}} functions meant for end user use.
#' There is also an S3 method for \code{\link[stats]{simulate}}, but with more
#' limited functionality than \code{\link[lme4]{simulate.merMod}}; the
#' \code{nauf.glmerMod} method for this function is not meant for end users, and
#' is necessary for the \code{PB} method to work for
#' \code{\link[=anova.nauf.merMod]{anova}}.
#' Other generics from the \code{base}, \code{stats}, and
#' \code{lme4} packages such as \code{summary} etc. work just as
#' they would for models fit with \code{\link[lme4]{glmer}} and
#' \code{\link[lme4]{glmer.nb}}.  If you encounter a generic function in these
#' packages which does not function properly, please report the issue at
#' \url{https://github.com/CDEager/nauf/issues}.
#'
#' @seealso \code{\link{nauf_glmer}}, \code{\link{nauf_glmer.nb}}, and
#'   \code{\link{nauf_contrasts}}.
#'
#' @slot rsp,Gp,call,frame,flist,cnms,lower,theta,beta,u,devcomp,pp,optinfo See
#'   \code{\linkS4class{merMod}}.
nauf.glmerMod <- setClass("nauf.glmerMod", contains = "glmerMod")


#' Class for reference grids for nauf models
#'
#' A list returned by \code{\link{nauf_ref.grid}} with one element
#' \code{ref.grid} which is a  \code{\linkS4class{ref.grid}} object. The 
#' reference grid the \code{nauf.ref.grid} contains should not be manipulated 
#' directly, or used as an argument to \code{\link[lsmeans]{lsmeans}}.  It 
#' should only be used with \code{\link{nauf_pmmeans}}.
#'
#' @aliases nauf.ref.grid
#' 
#' @name nauf.ref.grid-class
NULL


#' Class for \code{nauf} model anovas with nested models.
#'
#' A list returned by \code{\link{anova.nauf.merMod}} when \code{method} is one
#' of \code{nested-KR}, \code{LRT}, or \code{PB}. Similar to the object returned
#' by \code{\link[afex]{mixed}}, it has three elements: \code{anova_table} is
#' an \code{\link[stats]{anova}} table with the significance tests comparing
#' the full model to each of the nested models; \code{full_model} is the model
#' with all the fixed effects; and \code{restricted_models} is a named list
#' of nested models, each lacking one of the fixed effects.
#'
#' @seealso \code{\linkS4class{nauf.lmerMod}}, 
#'   \code{\linkS4class{nauf.glmerMod}}, and \code{\link{anova.nauf.merMod}}.
#'
#' @aliases nauf.nested.anova
#'
#' @name nauf.nested.anova-class
NULL

