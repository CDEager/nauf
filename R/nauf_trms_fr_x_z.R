

#' Not applicable unordered factor contrasts.
#'
#' The \code{nauf_contrasts} function returns a list of contrasts applied to
#' factors in an object created using a function in the \code{nauf} package.
#' See 'Details'.
#'
#' In the \code{nauf} package, \code{NA} values are used to encode when an
#' unordered factor is truly \emph{not applicable}.  This is different than
#' "not available" or "missing at random".  The concept applies only to
#' unordered factors, and indicates that the factor is simply not meaningful
#' for an observation, or that while the observation may technically be
#' definable by one of the factor levels, the interpretation of its belonging to
#' that level isn't the same as for other observations.  For imbalanced
#' observational data, coding unordered factors as \code{NA} may also be used to
#' control for a factor that is only contrastive within a subset of the data due
#' to the sampling scheme.  To understand the output of the
#' \code{nauf_contrasts} function, the treatment of unordered factor contrasts
#' in the \code{nauf} package will first be discussed, using the
#' \code{\link{plosives}} dataset included in the package as an example.
#'
#' In the \code{\link{plosives}} dataset, the factor \code{ling} is coded as
#' either \code{Monolingual}, indicating the observation is from a monolingual
#' speaker of Spanish, or \code{Bilingual}, indicating the observation is from a
#' Spanish-Quechua bilingual speaker.  The \code{dialect} factor indicates the
#' city the speaker is from (one of \code{Cuzco}, \code{Lima}, or
#' \code{Valladolid}).  The Cuzco dialect has both monolingual and bilingual
#' speakers, but the Lima and Valladolid dialects have only monolingual
#' speakers.  In the case of Valladolid, the dialect is not in contact with
#' Quechua, and so being monolingual in Valladolid does not mean the same
#' thing as it does in Cuzco, where it indicates
#' \emph{monolingual as opposed to bilingual}.  Lima has Spanish-Quechua
#' bilingual speakers, but the research questions the dataset serves to answer
#' are specific to monolingual speakers of Spanish in Lima.  If we leave the
#' \code{ling} factor coded as is in the dataset and use
#' \code{\link[standardize]{named_contr_sum}} to create the contrasts, we obtain
#' the following:
#'
#' \tabular{llrrr}{
#' dialect    \tab ling        \tab dialectCuzco \tab dialectLima \tab lingBilingual \cr
#' Cuzco      \tab Bilingual   \tab            1 \tab           0 \tab             1 \cr
#' Cuzco      \tab Monolingual \tab            1 \tab           0 \tab            -1 \cr
#' Lima       \tab Monolingual \tab            0 \tab           1 \tab            -1 \cr
#' Valladolid \tab Monolingual \tab           -1 \tab          -1 \tab            -1
#' }
#'
#' With these contrasts, the regression coefficient \code{dialectLima} would not
#' represent the difference between the intercept and the mean of the Lima
#' dialect; the mean of the Lima dialect would be the
#' \code{(Intercept) + dialectLima - lingBilingual}.  The interpretation of the
#' \code{lingBilingual} coefficient is similarly affected, and the intercept
#' term averages over the predicted value for the non-existent groups of Lima
#' bilingual speakers and Valladolid bilingual speakers, losing the
#' interpretation as the corrected mean (insofar as there can be a corrected
#' mean in this type of imbalanced data).  With the \code{nauf} package, we can
#' instead code non-Cuzco speakers' observations as \code{NA} for the
#' \code{ling} factor (i.e. execute
#' \code{plosives$ling[plosives$dialect != "Cuzco"] <- NA}).  These \code{NA}
#' values are allowed to pass into the regression's model matrix, and are then
#' set to \code{0}, effectively creating the following contrasts:
#'
#' \tabular{llrrr}{
#' dialect    \tab ling        \tab dialectCuzco \tab dialectLima \tab lingBilingual \cr
#' Cuzco      \tab Bilingual   \tab            1 \tab           0 \tab             1 \cr
#' Cuzco      \tab Monolingual \tab            1 \tab           0 \tab            -1 \cr
#' Lima       \tab NA          \tab            0 \tab           1 \tab             0 \cr
#' Valladolid \tab NA          \tab           -1 \tab          -1 \tab             0
#' }
#'
#' Because sum contrasts are used, a value of \code{0} for a dummy variable
#' averages over the effect of the factor, and the coefficient
#' \code{lingBilingual} only affects the predicted value for observations where
#' \code{dialect = Cuzco}.  In a regression fit with these contrasts, the
#' coefficient \code{dialectLima} represents what it should, namely the
#' difference between the intercept and the mean of the Lima dialect, and the
#' intercept is again the corrected mean.  The \code{lingBilingual} coefficient
#' is now the difference between Cuzco bilingual speakers and the corrected mean
#' \emph{of the Cuzco dialect}, which is \code{(Intercept) + dialectCuzco}.
#' These \code{nauf} contrasts thus allow us to model all of the data in a
#' single model without sacrificing the interpretability of the results. In
#' sociolinguistics, this method is called \emph{slashing} due to the use of a
#' forward slash in GoldVarb to indicate that a factor is not applicable.
#'
#' This same methodology can be applied to other parts of the
#' \code{\link{plosives}} dataset where a factor's interpretation is the same
#' for all observations, but is only contrastive within a subset of the data due
#' to the sampling scheme.  The \code{age} and \code{ed} factors (speaker age
#' group and education level, respectively) are factors which can apply to
#' speakers regardless of their dialect, but in the dataset they are only
#' contrastive within the Cuzco dialect; all the Lima and Valladolid speakers
#' are 40 years old or younger with a university education (in the case of
#' Valladolid, the data come from an already-existing corpus; and in the case of
#' Lima, the data were collected as part of the same dataset as the Cuzco data,
#' but as a smaller control group).  These factors can be treated just as the
#' \code{ling} factor by setting them to \code{NA} for observations from Lima
#' and Valladolid speakers.  Similarly, there is no read speech data for the
#' Valladolid speakers, and so \code{spont} could be coded as \code{NA} for
#' observations from Valladolid speakers.
#'
#' Using \code{NA} values can also allow the inclusion of a random effects
#' structure which only applies to a subset of the data.  The
#' \code{\link{plosives}} dataset has data from both read (\code{spont = FALSE};
#' only Cuzco and Lima) and spontaneous (\code{spont = TRUE}; all three
#' dialects) speech.  For the read speech, there are exactly repeated measures
#' on 54 items, as indicated by the \code{item} factor.  For the
#' spontaneous speech, there are not exactly repeated measures, and so in this
#' subset, \code{item} is coded as \code{NA}.  In a regression fit using
#' \code{nauf_lmer}, \code{nauf_glmer}, or \code{nauf_glmer.nb} with \code{item}
#' as a grouping factor, the random effects model matrix is created for the read
#' speech just as it normally is, and for spontaneous speech observations all of
#' the columns are set to \code{0} so that the \code{item} effects only affect
#' the fitted values for read speech observations.  In this way, the noise
#' introduced by the read speech items can be accounted for while still
#' including all of the data in one model, and the same random effects for
#' \code{speaker} can apply to all observations (both read and spontaneous),
#' which will lead to a more accurate estimation of the fixed, speaker, and item
#' effects since more information is available than if the read and spontaneous
#' speech were analyzed in separate models.
#'
#' There are two situations in which unordered factors will need more than one set
#' of contrasts: (1) when an unordered factor with \code{NA} values interacts
#' with another unordered factor, and some levels are collinear with \code{NA};
#' and (2) when an unordered factor is included as a slope for a random effects
#' grouping factor that has \code{NA} values, but only a subset of the levels
#' for the slope factor occur when the grouping factor is not \code{NA}.  As an
#' example of an interaction requiring new contrasts, consider the interaction
#' \code{dialect * spont} (that is, suppose we are interested in whether the
#' effect of \code{spont} is different for Cuzco and Lima).  We code
#' \code{spont} as \code{NA} when \code{dialect = Valladolid}, as mentioned
#' above.  This gives the following contrasts for the main effects:
#'
#' \tabular{llrrr}{
#' dialect    \tab spont \tab dialectCuzco \tab dialectLima \tab spontTRUE \cr
#' Cuzco      \tab TRUE  \tab            1 \tab           0 \tab         1 \cr
#' Cuzco      \tab FALSE \tab            1 \tab           0 \tab        -1 \cr
#' Lima       \tab TRUE  \tab            0 \tab           1 \tab         1 \cr
#' Lima       \tab FALSE \tab            0 \tab           1 \tab        -1 \cr
#' Valladolid \tab NA    \tab           -1 \tab          -1 \tab         0
#' }
#'
#' If we simply multiply these \code{dialect} and \code{spont} main effect
#' contrasts together to obtain the contrasts for the interaction (which is what
#' is done in the default \code{\link[stats]{model.matrix}} method), we get
#' following contrasts:
#'
#' \tabular{llrr}{
#' dialect    \tab spont \tab dialectCuzco:spontTRUE \tab dialectLima:spontTRUE \cr
#' Cuzco      \tab TRUE  \tab                      1 \tab                     0 \cr
#' Cuzco      \tab FALSE \tab                     -1 \tab                     0 \cr
#' Lima       \tab TRUE  \tab                      0 \tab                     1 \cr
#' Lima       \tab FALSE \tab                      0 \tab                    -1 \cr
#' Valladolid \tab NA    \tab                      0 \tab                     0
#' }
#'
#' However, these contrasts introduce an unnecessary parameter to the model
#' which causes collinearity with the main effects since
#' \code{spontTRUE = dialectCuzco:spontTRUE + dialectLima:spontTRUE} in all
#' cases.  The functions in the \code{nauf} package automatically recognize when
#' this occurs, and create a second set of contrasts for \code{dialect} in which
#' the \code{Valladolid} level is treated as if it were \code{NA} (through and
#' additional call to \code{\link[standardize]{named_contr_sum}}):
#'
#' \tabular{lr}{
#' dialect    \tab dialect.c2.Cuzco \cr
#' Cuzco      \tab                1 \cr
#' Lima       \tab               -1 \cr
#' Valladolid \tab                0
#' }
#'
#' This second set of \code{dialect} contrasts is only used when it needs to be.
#' That is, in this case, these contrasts would be used in the creation of the
#' model matrix columns for the interaction term \code{dialect:spont} term,
#' but not in the creation of the model matrix columns for the main effect terms
#' \code{dialect} and \code{spont}, and when the second set of contrasts is
#' used, \code{.c2.} will appear between the name of the factor and the level so
#' it can be easily identified:
#'
#' \tabular{llrrrr}{
#' dialect    \tab spont \tab dialectCuzco \tab dialectLima \tab spontTRUE \tab dialect.c2.Cuzco:spontTRUE \cr
#' Cuzco      \tab TRUE  \tab            1 \tab           0 \tab         1 \tab                          1 \cr
#' Cuzco      \tab FALSE \tab            1 \tab           0 \tab        -1 \tab                         -1 \cr
#' Lima       \tab TRUE  \tab            0 \tab           1 \tab         1 \tab                         -1 \cr
#' Lima       \tab FALSE \tab            0 \tab           1 \tab        -1 \tab                          1 \cr
#' Valladolid \tab NA    \tab           -1 \tab          -1 \tab         0 \tab                          0
#' }
#'
#' Turning now to an example of when a random slope requires new contrasts,
#' consider a random \code{item} slope for \code{dialect}.  Because
#' \code{dialect = Valladolid} only when \code{item} is \code{NA}, using the
#' main effect contrasts for \code{dialect} for the \code{item} slope would
#' result in collinearity with the \code{item} intercept in the random effects
#' model matrix:
#'
#' \tabular{llrrr}{
#' dialect    \tab item \tab i01:(Intercept) \tab i01:dialectCuzco \tab i01:dialectLima \cr
#' Cuzco      \tab i01  \tab               1 \tab                1 \tab               0 \cr
#' Cuzco      \tab i02  \tab               0 \tab                0 \tab               0 \cr
#' Cuzco      \tab NA   \tab               0 \tab                0 \tab               0 \cr
#' Lima       \tab i01  \tab               1 \tab                0 \tab               1 \cr
#' Lima       \tab i02  \tab               0 \tab                0 \tab               0 \cr
#' Lima       \tab NA   \tab               0 \tab                0 \tab               0 \cr
#' Valladolid \tab NA   \tab               0 \tab                0 \tab               0
#' }
#'
#' This table shows the random effects model matrix for \code{item i01} for all
#' possible scenarios, with the rows corresponding to (in order): a Cuzco
#' speaker producing the read speech plosive in \code{item i01}, a Cuzco speaker
#' producing a read speech plosive in another \code{item}, a Cuzco speaker
#' producing a spontaneous speech plosive, a Lima speaker producing the read
#' speech plosive in \code{item i01}, a Lima speaker producing a read speech
#' plosive in another \code{item}, a Lima speaker producing a spontaneous speech
#' plosive, and a Valladolid speaker producing a spontaneous speech plosive.
#' With the main effect contrasts for \code{dialect},
#' \code{i01:(Intercept) = i01:dialectCuzco + i01:dialectLima} in all cases,
#' causing collinearity.  Because this collinearity exists for all read speech
#' item random effects model matrices, the model is unidentifiable.  The
#' functions in the \code{nauf} package automatically detect that this is the
#' case, and remedy the situation by creating a new set of contrasts used for
#' the \code{item} slope for \code{dialect}:
#'
#' \tabular{llrr}{
#' dialect    \tab item \tab i01:(Intercept) \tab i01:dialect.c2.Cuzco \cr
#' Cuzco      \tab i01  \tab               1 \tab                    1 \cr
#' Cuzco      \tab i02  \tab               0 \tab                    0 \cr
#' Cuzco      \tab NA   \tab               0 \tab                    0 \cr
#' Lima       \tab i01  \tab               1 \tab                   -1 \cr
#' Lima       \tab i02  \tab               0 \tab                    0 \cr
#' Lima       \tab NA   \tab               0 \tab                    0 \cr
#' Valladolid \tab NA   \tab               0 \tab                    0
#' }
#'
#' If we were to, say, fit the model
#' \code{intdiff ~ dialect * spont + (1 + dialect | item)}, then \code{nauf} would
#' additionally recognize that the same set of altered contrasts for
#' \code{dialect} are required in the fixed effects interaction term
#' \code{dialect:spont} and the \code{item} slope for \code{dialect}, and both
#' would be labeled with \code{.c2.}.  In other (rare) cases, more than two sets
#' of contrasts may be required for a factor, in which case they would have
#' \code{.c3.}, \code{.c4.} and so on.
#'
#' In this way, users only need to code unordered factors as \code{NA} in the
#' subsets of the data where they are not contrastive, and \code{nauf} handles
#' the rest.  Having described in detail what \code{nauf} contrasts are, we now
#' return to the \code{nauf_contrasts} function.  The function can be used on
#' objects of any \code{nauf} model, a \code{\link{nauf.terms}} object, or a
#' model frame made by \code{\link{nauf_model.frame}}. It returns a named list
#' with a matrix for each
#' unordered factor in \code{object} which contains all contrasts associated the
#' factor.  For the model \code{intdiff ~ dialect * spont + (1 + dialect | item)},
#' the result would be a list with elements \code{dialect} and \code{spont} that
#' contain the following matrices (see the 'Examples' section for code to
#' generate this list):
#'
#' \tabular{lrrr}{
#' dialect    \tab Cuzco \tab Lima \tab .c2.Cuzco \cr
#' Cuzco      \tab     1 \tab    0 \tab         1 \cr
#' Lima       \tab     0 \tab    1 \tab        -1 \cr
#' Valladolid \tab    -1 \tab   -1 \tab         0
#' }
#'
#' \tabular{lr}{
#' spont \tab TRUE \cr
#' TRUE  \tab    1 \cr
#' FALSE \tab   -1 \cr
#' NA    \tab    0
#' }
#'
#' The default is for the list of contrasts to only contain information about
#' unordered factors.  If \code{inc_ordered = TRUE}, then the contrast matrices
#' for any ordered factors in \code{object} are also included.
#'
#' @param object A \code{\link{nauf.terms}} object, a model frame made with
#'   \code{\link{nauf_model.frame}}, a \code{nauf.glm} model (see
#'   \code{\link{nauf_glm}}), or a \code{\linkS4class{nauf.lmerMod}} or
#'   \code{\linkS4class{nauf.glmerMod}} model.
#' @param inc_ordered A logical indicating whether or not ordered factor
#'   contrasts should also be returned (default \code{FALSE}).
#'
#' @return A named list of contrasts for all unordered factors in \code{object},
#'   and also optionally contrasts for ordered factors in \code{object}.  See
#'   'Details'.
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#'
#' mf <- nauf_model.frame(intdiff ~ dialect * spont + (1 + dialect | item), dat)
#' nauf_contrasts(mf)
#'
#' mf <- nauf_model.frame(intdiff ~ dialect * spont + (1 + dialect | item),
#'   dat, ncs_scale = 0.5)
#' nauf_contrasts(mf)
#'
#' @section Note: The argument \code{ncs_scale} changes what value is used for
#'   the sum contrast deviations.  The default value of \code{1} would give the
#'   contrast matrices in 'Details'.  A value of \code{ncs_scale = 0.5}, for example,
#'   would result in replacing \code{1} with \code{0.5} and \code{-1} with
#'   \code{-0.5} in all of the contrast matrices.
#'
#' @seealso \code{\link{nauf_model.frame}}, \code{\link{nauf_model.matrix}},
#'   \code{\link{nauf_glFormula}}, \code{\link{nauf_glm}}, and
#'   \code{\link{nauf_glmer}}.
#'
#' @export
nauf_contrasts <- function(object, inc_ordered = FALSE) {
  info <- nauf.info(object)

  contr <- mlapply(levs = info$uf, hasna = info$hasna[names(info$uf)],
    same = list(ncs = info$ncs_scale), fun = expand_contr)
  
  if (inc_ordered && length(info$of)) {
    contr <- c(contr, lapply(info$of, `[[`, "contrasts"))
  }

  return(contr)
}


expand_contr <- function(levs, ncs, hasna) {
  contr <- lapply(levs, standardize::named_contr_sum, scale = ncs)
  
  if ((n <- length(contr)) > 1) {
    contr[2:n] <- mlapply(mat = contr[2:n], cj = 2:n,
      same = list(rn = rownames(contr[[1]])), fun = add_contr_zeros)
  }
  
  contr <- do.call(cbind, contr)
  if (hasna) {
    contr <- rbind(contr, 0)
    rownames(contr)[nrow(contr)] <- NA
  }

  return(contr)
}


add_contr_zeros <- function(mat, cj, rn) {
  cn <- paste0(".c", cj, ".", colnames(mat))
  expanded <- matrix(0, length(rn), length(cn), dimnames = list(rn, cn))
  expanded[rownames(mat), ] <- mat
  return(expanded)
}


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


#' @export
model.frame.nauf.formula <- function(formula, data = NULL, subset = NULL,
                                   na.action = na.pass,
                                   drop.unused.levels = TRUE, xlev = NULL,
                                   ncs_scale = attr(formula, "standardized.scale"),
                                   ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf::nauf_model.frame)
  return(eval(mc, parent.frame()))
}


#' @export
model.frame.nauf.terms <- function(formula, data = NULL, subset = NULL,
                                   na.action = na.pass,
                                   drop.unused.levels = TRUE, xlev = NULL,
                                   ncs_scale = attr(formula, "standardized.scale"),
                                   ...) {
  mc <- match.call()
  
  info <- nauf.info(formula)
  drop_class(formula) <- "nauf.terms"
  ncs <- info$ncs_scale
  
  mc[[1]] <- quote(stats::model.frame)
  mc$formula <- formula
  mc$na.action <- na.pass
  mc$drop.unused.levels <- TRUE
  mc$xlev <- NULL
  mc["ncs_scale"] <- NULL
  mf <- eval(mc, parent.frame())
  
  attr(mf, "formula") <- stats::formula(attr(mf, "terms"))
  nauf.info(mf) <- info
  first_class(mf, "terms") <- "nauf.terms"
  first_class(mf, "formula") <- "nauf.formula"
  last_class(mf) <- "nauf.frame"
  
  v <- colnames(mf)
  uf <- intersect(v, names(info$uf))
  of <- intersect(v, names(info$of))
  groups <- intersect(v, names(info$groups))
  
  lvs <- lapply(info$uf[uf], `[[`, 1)
  contr <- lapply(lvs, standardize::named_contr_sum, scale = ncs)
  mf[uf] <- mlapply(x = mf[uf], levels = lvs, contrasts = contr,
    fun = standardize::fac_and_contr)
  
  mf[of] <- mlapply(x = mf[of], levels = lapply(info$of[of], `[[`, "levels"),
    contrasts = lapply(info$of[of], `[[`, "contrasts"),
    same = list(ordered = TRUE), fun = standardize::fac_and_contr)
  
  if (isTRUE(info$allow.new.levels)) {
    mf[groups] <- mlapply(fac = mf[groups], levs = info$groups[groups],
      hasna = info$hasna[groups], fun = function(fac, levs, hasna) {
        wnew <- which(!(fac %in% c(levs, if (hasna) NA)))
        fac <- factor(fac, levels = c(levs, "_NEW_"))
        fac[wnew] <- "_NEW_"
        return(fac)
      }
    )
  } else {
    mf[groups] <- mlapply(x = mf[groups], levels = info$groups[groups],
      same = list(ordered = FALSE), fun = factor)
  }
  
  if (any(sapply(mf, anyNA) & !info$hasna[v])) {
    warning("Some variables which did not have NA values when the model was ",
      "fit have NA values in the new model frame.")
  }

  return(mf)
}


#' @export
model.matrix.nauf.terms <- function(object, data = environment(object),
                                    contrasts.arg = NULL, xlev = NULL, ...) {
  if (!is.nauf.frame(data)) {
    data <- model.frame(object, data)
  }
  return(nauf_mm(data))
}


#' Create a model frame using \code{nauf} contrasts.
#'
#' \code{nauf_model.frame} creates a model frame which employs
#' \code{\link{nauf_contrasts}} for unordered factors.
#'
#' First, the default method for \code{\link[stats]{model.frame}} is called.
#' Then any variable in the resulting model frame that is not an unordered
#' factor but has only two unique non-\code{NA} values is coerced to an
#' unordered factor.  Unordered factors are then assigned contrasts with
#' \code{\link[standardize]{named_contr_sum}}, passing \code{ncs_scale} as the
#' function's \code{scale} argument.  Then, necessary contrast changes in
#' interaction terms and random effects slopes are determined as described in
#' \code{\link{nauf_contrasts}}.
#'
#' The recommended usage is to first \code{\link[standardize]{standardize}} the
#' regression variables, and then use the \code{formula} and \code{data}
#' elements in the resulting \code{standardized} object as arguments to
#' \code{nauf_model.frame}.  When this is done, \code{ncs_scale} is obtained
#' from the \code{standardized.scale} attribute of the \code{formula}, unless
#' \code{ncs_scale} is specified as a value which does not match the
#' \code{standardized} scale, in which case the explicitly specified
#' \code{ncs_scale} argument is used with a warning.  If
#' \code{\link[standardize]{standardize}} is not used prior to calling
#' \code{nauf_model.frame}, then \code{ncs_scale} defaults to \code{1} unless
#' explicitly specified in the function call, in which case the specified value
#' is used.
#'
#' Changes from the following default values are ignored with a warning:
#' \describe{
#'   \item{na.action = na.pass}{This default value is required in order for
#'     \code{NA} values to be treated as defined in
#'     \code{\link{nauf_contrasts}}.}
#'   \item{drop.unused.levels = TRUE}{This default value is set because
#'     \code{nauf_model.frame} assumes that \code{data} is not new data.  To
#'     create a \code{nauf.frame} with new data, the \code{terms}
#'     attribute of an already existing \code{nauf.frame} (which
#'     has class \code{\link{nauf.terms}}) can be used as the
#'     \code{formula} argument to \code{\link[stats]{model.frame}}.}
#'   \item{xlev = NULL}{This default is necessary for the same reasons as the
#'     default value for \code{drop.unused.levels}.}
#'   \item{contrasts = NULL}{For unordered factors, contrasts are automatically
#'     created with \code{\link[standardize]{named_contr_sum}}, as sum contrasts
#'     are necessary to implement \code{\link{nauf_contrasts}}.  To specify
#'     custom contrasts for ordered factors, the custom contrasts should be
#'     explicitly assigned to the ordered factor in \code{data} (this is
#'     automatically done if \code{\link[standardize]{standardize}} is used
#'     first as recommended).}
#' }
#'
#' @param formula,data,subset,... See \code{\link[stats]{model.frame}}.
#' @param na.action,drop.unused.levels,xlev,contrasts Changes from default
#'   values for these arguments are ignored with a warning.
#' @param ncs_scale A positive number passed as the \code{scale} argument to
#'   \code{\link[standardize]{named_contr_sum}} for all unordered factor
#'   contrasts.  The default is to first check whether \code{formula} comes from
#'   a \code{standardized} object returned by
#'   \code{\link[standardize]{standardize}}. If it is, then the \code{scale}
#'   argument from the \code{\link[standardize]{standardize}} call is used.  If
#'   it is not, then \code{ncs_scale} is set to \code{1}.  The value for
#'   \code{ncs_scale} can also be set explicitly.  If it is set explicitly and
#'   \code{formula} is from a \code{standardized} object with a different scale
#'   than the explicitly set value, then the explicitly set value
#'   is used and a warning is issued.
#'
#' @return A model frame with second class attribute \code{nauf.frame}.  Its
#'   \code{formula} attribute has class \code{nauf.formula} and its \code{terms}
#'   attribute has class \code{\link{nauf.terms}}.
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' form <- intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item)
#'
#' ## default behavior when standardize is not used
#' # defaults to ncs_scale = 1
#' mf <- nauf_model.frame(form, dat)
#'
#' # uses specified ncs_scale = 0.5
#' mf_0.5 <- nauf_model.frame(form, dat, ncs_scale = 0.5)
#'
#' ## standardize first (recommended use)
#' sdat <- standardize(form, dat)
#' sdat_0.5 <- standardize(form, dat, scale = 0.5)
#'
#' # uses ncs_scale = 1 from attr(sdat$formula, "standardized.scale")
#' mf_sdat <- nauf_model.frame(sdat$formula, sdat$data)
#'
#' # uses ncs_scale = 0.5 from attr(sdat_0.5$formula, "standardized.scale")
#' mf_sdat_0.5 <- nauf_model.frame(sdat_0.5$formula, sdat_0.5$data)
#'
#' \dontrun{
#' ## not recommended
#' # uses specified ncs_scale = 0.5 and issues a warning since
#' # attr(sdat$formula, "standardized.scale") = 1
#' mf_warning <- nauf_model.frame(sdat$formula, sdat$data, ncs_scale = 0.5)
#' }
#'
#' @seealso \code{\link{nauf_contrasts}} for a description of the contrasts
#'   applied to unordered factors, \code{\link{nauf_model.matrix}} for obtaining
#'   a fixed effects model matrix, and \code{\link{nauf_glFormula}} for
#'   obtaining both fixed effects and random effects model matrices.
#'
#' @export
nauf_model.frame <- function(formula, data = NULL, subset = NULL,
                             na.action = na.pass, drop.unused.levels = TRUE,
                             xlev = NULL, contrasts = NULL,
                             ncs_scale = attr(formula, "standardized.scale"),
                             ...) {
  # Ignore na.action, contrasts, drop.unused.levels, and xlev
  mc <- match.call()
  mc$na.action <- na.pass
  mc$drop.unused.levels <- TRUE
  mc$xlev <- NULL
  mc$contrasts <- NULL
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action', must be na.pass")
  }
  if (!isTRUE(drop.unused.levels)) {
    warning("Ignoring 'drop.unused.levels', must be TRUE")
  }
  if (!is.null(xlev)) {
    warning("Ignoring 'xlev', must be NULL")
  }
  if (!is.null(contrasts)) {
    warning("Ignoring 'contrasts', must be NULL")
  }

  standardized_scale <- attr(formula, "standardized.scale")
  if (is.null(ncs <- ncs_scale)) ncs <- 1
  if (!is.scalar(ncs, 1)) {
    stop("The scale for sum contrasts must be a single positive number")
  }
  if (!is.null(standardized_scale) && standardized_scale != ncs) {
    warning("'formula' is from a standardized object with scale ",
      standardized_scale, " but ncs_scale was specified as ", ncs)
  }

  formula <- stats::formula(formula)
  class(formula) <- "formula"
  bars <- lme4::findbars(formula)
  sb_form <- lme4::subbars(formula)
  if (stats::is.empty.model(sb_form)) stop("There are no predictors in 'formula'")
  fe_form <- stats::terms(lme4::nobars(formula))
  if (!attr(fe_form, "response")) stop("'formula' must have a response")
  if (!attr(fe_form, "intercept")) {
    stop("There must be a fixed effects intercept")
  }
  fe_form <- stats::formula(stats::delete.response(fe_form))
  groups <- check_groups(formula)

  if (!is.data.frame(data)) stop("'data' must be a data.frame")

  mc$formula <- sb_form
  mc["ncs_scale"] <- NULL
  mc[[1]] <- quote(stats::model.frame)
  mf <- eval(mc, parent.frame())
  mt <- attr(mf, "terms")

  fmat <- attr(mt, "factors")
  fmat <- fmat[setdiff(rownames(fmat), groups), , drop = FALSE]
  if (any(fmat > 1)) {
    warning("'nauf' has not been tested for models that violate the ",
      "interaction hierarchy")
  }

  cnms <- colnames(mf)
  extras <- find_extras(mf)
  rgrp <- setNames(cnms %in% groups, cnms)
  mf[rgrp] <- lapply(mf[rgrp], factor, ordered = FALSE)
  check <- which(!extras & !rgrp)
  mf[check] <- lapply(mf[check], charlogbin_to_uf)

  uf <- sapply(mf, is.uf) & !extras & !rgrp
  of <- sapply(mf, is.ordered) & !extras
  mat <- sapply(mf, is.matrix) & !extras
  num <- sapply(mf, is.numeric) & !extras & !mat
  hasna <- sapply(mf, anyNA)
  uf[1] <- of[1] <- mat[1] <- num[1] <- FALSE
  
  if (any(hasna & !(uf | rgrp))) {
    stop("Only unordered factor predictors and random effects grouping factors",
      " can have NA values")
  }

  mf[uf] <- lapply(mf[uf], standardize::named_contr_sum, scale = ncs,
    return_contr = FALSE)
  attr(mt, "dataClasses")[cnms[uf | rgrp]] <- "factor"
  
  changes <- contrast_changes(fe_form, bars, mf, uf)
    
  mt <- nauf.terms(mt, resp = cnms[1], groups = lapply(mf[groups], levels),
    uf = changes$uf, of = lapply(mf[of], levs_and_contr),
    num = lapply(mf[num], mean), mat = lapply(mf[mat], colMeans),
    extras = cnms[extras], cc = changes$cc, hasna = hasna, ncs_scale = ncs)
  first_class(formula) <- "nauf.formula"
  last_class(mf) <- "nauf.frame"
  attr(mf, "terms") <- mt
  attr(mf, "formula") <- formula

  return(mf)
}


nauf.terms <- function(terms, ...) {
  first_class(terms) <- "nauf.terms"
  attr(terms, "nauf.info") <- list(...)
  return(terms)
}


contrast_changes <- function(fixed, bars, mf, uf) {
  ufn <- names(uf)[uf]
  uf <- lapply(mf[ufn], function(x) list(levels(x)))
  changes <- lapply(c(list(fixed), bars), .contrast_changes, mf = mf)
  main <- lapply(changes, `[[`, "lvs")
  inter <- lapply(changes, `[[`, "cc")
  asgn <- lapply(changes, `[[`, "asgn")
  
  # do.call(main) and do.call(c, do.call(c, inter)) are lists of charvecs
  # named by uf; combine all into all_levs
  all_levs <- nsplit(c(do.call(c, main), do.call(c, do.call(c, inter))))[ufn]
  
  # uf is now a named list of unique contrast set levels for each factor
  uf <- mlapply(mec = uf, changelevs = all_levs, fun = function(mec, changelevs)
    unique(c(list(mec), changelevs)))
  
  # convert to named numeric vectors of contrast references
  main <- lapply(main, contr_nums, levlist = uf)
  inter <- lapply(inter, function(x) mlapply(levs = x, same = list(levlist = uf),
    fun = contr_nums))
  
  # join so one element per form (fixed + bars)
  cc <- mlapply(m = main, i = inter, a = asgn, fun = function(m, i, a)
    c(m, mlapply(factors = i, assign = a, fun = list)))
  
  return(rsa_list(uf, cc))
}


# when contrast_changes is called on a ranef bar, we only care if the contrasts
# change in interactions from the main effect contrasts *for the bar* because
# if the main effect contrasts are not .c1., they will have been changed in
# mf prior to calling ccmat
.contrast_changes <- function(form, mf, uf, hasna) {
  lvs <- cc <- asgn <- list()
  
  if (re <- !inherits(form, "formula")) {
    group <- varnms(barform(form, 3))
    if (is.null(varnms(form <- barform(form, 2)))) {
      return(rsa_list(lvs, cc, asgn))
    }
    if (re <- any(hasna[group])) {
      mf <- droplevels(mf[!rowna(mf[group]), , drop = FALSE])
    }
  }
  
  fmat <- attr(stats::terms(form), "factors") > 0
  rn <- rownames(fmat)
  uf <- uf[rn]
  hasna <- hasna[rn]
  nauf <- uf & hasna
  ufmat <- fmat[uf, , drop = FALSE]
  naufmat <- fmat[nauf, , drop = FALSE]
  check_inter <- length(inter <- which(colSums(ufmat) > 1 & colSums(naufmat)))
  check_main <- re & length(main <- intersect(rn[uf], colnames(fmat)))

  if (check_main) {
    lvs <- lapply(mf[main], levels)
  }
  
  if (check_inter) {
    cc <- unique(lapply(inter, function(x) sort(rownames(ufmat)[ufmat[, x]])))
    cc <- mlapply(cols = cc, same = list(x = mf), fun = nauf_interaction)
    changed <- sapply(cc, `[[`, "changed")
    cc <- lapply(cc[changed], `[[`, "levels")
    facs <- lapply(cc, names)
    asgn <- mlapply(facs, same = list(m = fmat), fun = function(f, m)
      which(sapply(m[f, , drop = FALSE], 2, all)))
    names(asgn) <- names(cc) <- sapply(facs, paste, collapse = ":")
  }
  
  return(rsa_list(lvs, cc, asgn))
}


contr_nums <- function(levs, levlist) {
  nums <- mapply(in_list, x = levs, lst = levlist[names(levs)])
  if (length(nums) && (!is.numeric(nums) || any(nums == 0))) {
    stop("failed to create named vector contrast reference numbers")
  }
  return(nums)
}


nauf_interaction <- function(x, cols) {
  nm <- paste(cols, collapse = ":")
  x <- x[cols]
  mlvs <- lapply(x, function(n) reorder_ft(sort(levels(n))))

  # remove any unique combination which involves NAs
  x <- unique(x)
  x <- droplevels(x[!rowna(x), , drop = FALSE])
  if (!nrow(x)) stop("No unique applicable combinations in ", nm)

  # remove levels which are redundant
  # e.g. f1 in [a, b, c], f2 in [d, e, f, g]
  #    when f1 = a, f2 in [d, e, f]
  #    when f1 = b, f2 in    [e, f, g]
  #    when f1 = c, f2 = NA
  #    then at this point [c] has been dropped form f1,
  #    but we still need to drop [d, g] from f2
  if (empty_cells(x)) {
    torm <- mlapply(i = lapply(x, levels), j = 1:ncol(x), same = list(mat = x),
      fun = function(j, i, mat) {
        do.call(c, mlapply(lev = i, same = list(n = j, m = mat),
          fun = function(lev, n, m) {
            check <- droplevels(m[m[[n]] %in% lev, -n, drop = FALSE])
            if (any(sapply(check, nlevels) == 1)) return(lev)
            return(NULL)
          }
        ))
      }
    )
    
    if (length(torm <- torm[lengths(torm) > 0])) {
      f <- names(torm)
      x[f] <- mlapply(fac = x[f], levs = torm, fun = function(fac, levs) {
        fac[fac %in% levs] <- NA
        return(fac)
      })
    }

    x <- unique(droplevels(x[!rowna(x), , drop = FALSE]))
    if (!nrow(x)) stop("No unique applicable combinations in ", nm)
  }

  ilvs <- lapply(x, function(n) reorder_ft(sort(levels(n))))
  changed <- !isTRUE(all.equal(mlvs, ilvs))
  if (any(unlist(lapply(x, nlevels)) < 2)) {
    stop("At least one factor in ", nm,
      " has only one level when NAs are removed")
  }

  # there can still be collinearity, but in this case it isn't structural
  # so warning rather than error (i.e. a column will be dropped from
  # the model matrix if the interaction is included in a regression)
  if (empty_cells(x)) warning("Collinearity in the interaction ", nm)

  return(list(levels = ilvs, changed = changed))
}


#' Create a fixed effects model matrix using \code{nauf} contrasts.
#'
#' \code{nauf_model.matrix} creates a model matrix which employs
#' \code{\link{nauf_contrasts}} for unordered factors.
#'
#' Exactly what happens depends on the values of \code{object} and \code{data}.
#' The following possibilities are evaluated in the order listed:
#' \describe{
#'   \item{object is a nauf.frame}{All arguments besides \code{object} are
#'     ignored, and the information in \code{object} is used to create the model
#'     matrix.}
#'   \item{data is a nauf.frame}{All arguments besides \code{data} are ignored,
#'     and the information in \code{data} is used to create the model matrix.}
#'   \item{object is a formula and data is a data.frame}{
#'     \code{\link{nauf_model.frame}} is called with \code{formula = object}
#'     and \code{data = data}, passing along any additional arguments in
#'     \code{...} (including \code{ncs_scale}).  Then the model matrix is
#'     created using the information in the resulting
#'     \code{nauf.frame}.}
#'   \item{any other argument values}{An error is returned.}
#' }
#'
#' @param object A \code{nauf.frame} or a regression formula.
#'   See 'Details'.
#' @param data A \code{nauf.frame} or a \code{data.frame}
#'   containing the variables in \code{object} if \code{object} is a regression
#'   formula. See 'Details'.
#' @param ... Further arguments to be passed to \code{\link{nauf_model.frame}}
#'   when \code{object} is a regression formula and \code{data} is a
#'   \code{data.frame}. See 'Details'.
#'
#' @return A fixed effects model matrix that implements
#'   \code{\link{nauf_contrasts}}.  Unlike the default
#'   \code{\link[stats]{model.matrix}} method, the model matrix does not have a
#'   \code{contrasts} attribute, since multiple sets of contrasts may be
#'   required for some unordered factors.
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' form <- intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item)
#' sdat <- standardize(form, dat)
#' mf <- nauf_model.frame(sdat$formula, sdat$data)
#'
#' ## the following all result in the same model matrix
#' mm1 <- nauf_model.matrix(mf)
#' mm2 <- nauf_model.matrix(form, mf)  # 'form' ignored
#' mm3 <- nauf_model.matrix(sdat$formula, sdat$data)
#'
#' @seealso \code{\link{nauf_contrasts}} for a description of the contrasts
#'   applied to unordered factors, \code{\link{nauf_model.frame}} for creating a
#'   model frame with \code{nauf} contrasts, and \code{\link{nauf_glFormula}}
#'   for obtaining both fixed effects and random effects model matrices.
#'
#' @export
nauf_model.matrix <- function(object = NULL, data = NULL, ...) {
  mc <- match.call()

  if (is.nauf.frame(object)) return(nauf_mm(object))
  if (is.nauf.frame(data)) return(nauf_mm(data))
  if (inherits(object, "formula")) {
    names(mc)[2] <- "formula"
    mc[[1]] <- quote(nauf::nauf_model.frame)
    mf <- eval(mc, parent.frame())
    return(nauf_mm(mf))
  }

  stop("If 'object' is not a formula, then either 'object' or 'data' must be\n",
    "  a model frame created by nauf_model.frame")
}


nauf_mm <- function(mf, ccn = 1) {
  attr(mf, "na.action") <- "na.pass"
  formula <- attr(mf, "formula")
  mt <- attr(mf, "terms")
  info <- attr(mt, "nauf.info")
  ncs <- info$ncs_scale
  ufc <- info$uf
  cc <- info$cc[[ccn]]

  if (ccn == 1) {
    formula <- stats::delete.response(stats::terms(lme4::nobars(formula)))
  } else {
    formula <- stats::terms(barform(lme4::findbars(formula)[[ccn - 1]], 2))
    ccmain <- intersect(names(ufc), names(cc))
    if (length(ccmain)) {
      mf[ccmain] <- mlapply(fac = mf[ccmain], levs = ufc[ccmain], cj = cc[ccmain],
        same = list(ncs = ncs), fun = apply_contrast_changes)
      cc <- cc[-which(names(cc) %in% ccmain)]
    }
  }

  mm <- stats::model.matrix(formula, mf)

  if (length(cc)) {
    mmlist <- list()
    cnms <- character()
    asgn <- attr(mm, "assign")
    asgn_cc <- sort(unique(unlist(lapply(cc, `[[`, "assign"))))
    mmrm <- which(asgn %in% asgn_cc)
    asgn <- asgn[-mmrm]
    if (length(asgn)) {
      mmlist[[1]] <- mm[, -mmrm, drop = FALSE]
      cnms <- colnames(mmlist[[1]])
    }
    fmat <- attr(formula, "factors")
    
    ccmms <- lapply(cc, ccmat, mf = mf, ufc = ufc, ncs = ncs, fmat = fmat)
    mm <- do.call(cbind, c(mmlist, lapply(ccmms, `[[`, "matrix")))
    asgn <- c(asgn, unlist(lapply(ccmms, `[[`, "assign")))
    names(asgn) <- c(cnms, colnames(mm))
    asgn <- sort(asgn)
    mm <- mm[, names(asgn), drop = FALSE]
    names(asgn) <- NULL
    attr(mm, "assign") <- asgn
  }

  attr(mm, "contrasts") <- NULL
  mm[is.na(mm)] <- 0

  return(mm)
}


apply_contrast_changes <- function(fac, levs, cj, ncs) {
  if (cj == 1) return(fac)
  levs <- levs[[cj]]
  fac <- factor(fac, ordered = FALSE, levels = levs)
  contr <- standardize::named_contr_sum(levs, ncs)
  colnames(contr) <- paste0(".c", cj, ".", colnames(contr))
  contrasts(fac) <- contr
  return(fac)
}


ccmat <- function(cc, mf, ufc, ncs, fmat) {
  uf <- names(cc$factors)
  
  mf[uf] <- mlapply(fac = mf[uf], levs = ufc[uf], cj = cc$factors,
    same = list(ncs = ncs), fun = apply_contrast_changes)
  
  fmat <- fmat[, cc$assign, drop = FALSE] > 0
  main <- rownames(fmat)[rowSums(fmat) > 0]
  fmat <- fmat[main, , drop = FALSE]
  mt <- stats::terms(stats::formula(do.call(paste, c(
    list(paste("~", paste(main, collapse = "+"))),
    lapply(list_mat_cols(fmat), function(x) paste(main[x], collapse = "*")),
    list(sep = "+")))))
  mm <- stats::model.matrix(mt, mf)
  
  asgn_mt <- which(colnames(attr(mt, "factors")) %in% colnames(fmat))
  asgn_mm <- attr(mm, "assign")
  keep <- which(asgn_mm %in% asgn_mt)
  mm <- mm[, keep, drop = FALSE]
  asgn_mm <- cc$assign[as.numeric(factor(asgn_mm[keep]))]
  attr(mm, "assign") <- NULL
  
  return(list(matrix = mm, assign = asgn_mm))
}


#' @importFrom Matrix rBind t sparseMatrix drop0 diag KhatriRao
nauf_mkReTrms <- function(fr, lvs = NULL) {
  # based on lme4::mkReTrms
  if (!is.nauf.frame(fr)) {
    stop("'fr' was not created with nauf_model.frame")
  }
  bars <- lme4::findbars(attr(fr, "formula"))
  if (!length(bars)) {
    stop("No random effects terms specified in formula", call. = FALSE)
  }
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr,
    "data.frame"))

  names(bars) <- lme4_barnames(bars)
  term.names <- vapply(bars, lme4_safeDeparse, "")
  blist <- mlapply(bar = bars, ccn = 1 + 1:length(bars), same = list(fr = fr,
    lvs = lvs), fun = nauf_mkBlist)
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (any(diff(nl) > 0)) {
    ord <- rev(order(nl))
    blist <- blist[ord]
    nl <- nl[ord]
    term.names <- term.names[ord]
  }

  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(Matrix::rBind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1)) / 2)
  nb <- nc * nl
  if (sum(nb) != q) {
    stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
      sum(nb), q))
  }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))

  Lambdat <- Matrix::t(do.call(Matrix::sparseMatrix, do.call(Matrix::rBind,
    lapply(seq_along(blist), function(i) {
      mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
      dd <- diag(nc[i])
      ltri <- lower.tri(dd, diag = TRUE)
      ii <- row(dd)[ltri]
      jj <- col(dd)[ltri]
      data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[, jj]) +
        boff[i], x = as.double(rep.int(seq_along(ii), rep.int(nl[i],
        length(ii))) + thoff[i]))
    }))))

  thet <- numeric(sum(nth))
  ll <- list(Zt = Matrix::drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x),
    Gp = unname(c(0L, cumsum(nb))))
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(Matrix::diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower)
  Lambdat@x[] <- ll$theta[ll$Lind]
  ll$Lambdat <- Lambdat
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else {
    asgn <- seq_along(fl)
  }
  names(fl) <- ufn
  fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist

  return(ll)
}


nauf_mkBlist <- function(bar, ccn, fr, lvs) {
  gvars <- varnms(barform(bar, 3))
  ff <- interaction(fr[gvars])

  if (is.null(lvs)) {
    ff <- droplevels(ff)
    if (all(is.na(ff))) {
      stop("Invalid grouping factor specification, ", deparse(bar[[b]][[3]]),
        call. = FALSE)
    }

  } else {  # implies predict method with new data
    ff <- factor(ff, levels = lvs[[paste(gvars, collapse = ":")]],
      ordered = FALSE)
  }

  mm <- nauf_mm(fr, ccn)
  sm <- Matrix::KhatriRao(Matrix::fac2sparse(ff, to = "d",
    drop.unused.levels = FALSE), t(mm))
  dimnames(sm) <- list(rep(levels(ff), each = ncol(mm)), rownames(mm))

  return(list(ff = ff, sm = sm, nl = nlevels(ff), cnms = colnames(mm)))
}


#' Create a model frame and fixed and random effects model matrices using \code{nauf} contrasts.
#'
#' The same as the \code{lme4} \code{\link[lme4]{modular}} functions
#' \code{glFormula} and \code{lFormula}, but implementing
#' \code{\link{nauf_contrasts}}.  \code{nauf_lFormula} is used for linear mixed
#' effects regressions (i.e. those that would be fit with
#' \code{\link{nauf_lmer}}) and \code{nauf_glFormula} is used for genarlized
#' linear mixed effects regressions (i.e. those that would be fit with
#' \code{\link{nauf_glmer}} or \code{\link{nauf_glmer.nb}}).  Both of the
#' functions contain a call to \code{nauf_mkReTrms}, which serves the same
#' purpose as the \code{lme4} function \code{\link[lme4]{mkReTrms}}, but with
#' \code{\link{nauf_contrasts}}, and, while \code{\link[lme4]{mkReTrms}} is
#' exported by \code{lme4}, \code{nauf_mkReTrms} is an internal function in the
#' \code{nauf} package.
#'
#' @param formula,data,family,REML,subset,weights,offset,control,mustart,etastart,...
#'   See \code{\link[lme4]{glFormula}}.
#' @param na.action,contrasts Changes from default values are ignored.  See
#'   \code{\link{nauf_model.frame}}.
#' @param ncs_scale A positive number to be passed as the \code{scale} argument
#'   to \code{\link[standardize]{named_contr_sum}} for all unordered factors.
#'   See \code{\link{nauf_model.frame}}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{fr}{The model frame (with class \code{nauf.frame}).
#'     See \code{\link{nauf_model.frame}}.}
#'   \item{X}{The fixed effects model matrix with \code{\link{nauf_contrasts}}
#'     applied. See \code{\link{nauf_model.matrix}}.}
#'   \item{reTrms}{A list containing the random effects model matrix and other
#'     information about the random effects structure.  The elements of the list
#'     have the same structure as that returned by \code{\link[lme4]{mkReTrms}},
#'     but incorportating \code{\link{nauf_contrasts}}.}
#'   \item{REML}{(\code{nauf_lFormula} only): A logical indicating if restricted
#'     maximum likelihood was used (copy of argument).}
#'   \item{family}{(\code{nauf_glFormula} only): The regression family (copy
#'     of argument).}
#'   \item{formula}{The \code{formula} argument, but with any double vertical
#'     bars expanded (e.g. \code{(1 + x || subj)} becomes
#'     \code{(1 | subj) + (0 + x | subj)}).}
#'   \item{wmsgs}{Warning messages (if any).}
#' }
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' dat_form <- intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item)
#' sdat <- standardize(dat_form, dat)
#' lmod <- nauf_lFormula(sdat$formula, sdat$data)
#'
#' vless <- droplevels(subset(dat, voicing == "Voiceless"))
#' vless$fully_voiced <- vless$vdur == 0
#' vless_form <- fully_voiced ~ dialect * spont +
#'   (1 + spont | speaker) + (1 + dialect | item)
#' svless <- standardize(vless_form, vless, family = binomial)
#' glmod <- nauf_glFormula(svless$formula, svless$data, family = binomial)
#'
#' @seealso \code{\link{nauf_contrasts}} for a description of the contrasts
#'   applied to unordered factors; \code{\link{nauf_model.frame}} and
#'   \code{\link{nauf_model.matrix}} for the creation of the \code{fr} and
#'   \code{X} elements of the returned list, respectively; and
#'   \code{\link{nauf_lmer}}, \code{\link{nauf_glmer.nb}}, and
#'   \code{\link{nauf_glmer}} for fitting mixed effects regressions with gaussian,
#'   negative binomial, and all other families, respectively.
#'
#' @export
nauf_glFormula <- function(formula, data = NULL, family = gaussian, subset,
                           weights, na.action = na.pass, offset,
                           contrasts = NULL, mustart, etastart,
                           control = lme4::glmerControl(),
                           ncs_scale = attr(formula, "standardized.scale"),
                           ...) {
  # based on lme4::glFormula
  control <- control$checkControl
  mf <- mc <- match.call()

  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  if (is.linear(family <- get_family(family))) {
    mc[[1]] <- quote(nauf::nauf_lFormula)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  } else if (!(is.character(family) && family == "negbin")) {
    if (family$family %in% c("quasibinomial", "quasipoisson", "quasi")) {
      stop("\"quasi\" families cannot be used in glmer")
    }
  }

  ignoreArgs <- c("start", "verbose", "devFunOnly", "optimizer", "control",
    "nAGQ")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4_checkArgs, c(list("glmer"), l...))

  cstr <- "check.formula.LHS"
  lme4_checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4_checkFormulaData(formula, data,
    checkLHS = control$check.formula.LHS == "stop")
  mc$formula <- formula <- stats::as.formula(formula, env = denv)
  m <- match(c("data", "subset", "weights", "offset", "ncs_scale",
    "mustart", "etastart"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(nauf::nauf_model.frame)

  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- formula
  fr <- eval(mf, parent.frame())
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)

  reTrms <- nauf_mkReTrms(fr)
  wmsgNlev <- lme4_checkNlevels(reTrms$flist, n = n, control, allow.n = TRUE)
  wmsgZdims <- lme4_checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
  wmsgZrank <- lme4_checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06,
    allow.n = TRUE)

  mf[[1L]] <- quote(stats::model.frame)
  fixedform <- formula
  lme4_RHSForm(fixedform) <- lme4::nobars(lme4_RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr,
    "terms"), "predvars")

  ranform <- formula
  lme4_RHSForm(ranform) <- lme4::subbars(lme4_RHSForm(lme4_reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), "predvars")

  X <- nauf_model.matrix(fr)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4_chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4_checkScaleX(X, kind = scaleX.chk)

  return(list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}


#' @rdname nauf_glFormula
#' @export
nauf_lFormula <- function(formula, data = NULL, REML = TRUE, subset, weights,
                          na.action = na.pass, offset, contrasts = NULL,
                          control = lme4::lmerControl(),
                          ncs_scale = attr(formula, "standardized.scale"),
                          ...) {
  # based on lme4::lFormula
  control <- control$checkControl
  mf <- mc <- match.call()

  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  ignoreArgs <- c("start", "verbose", "devFunOnly", "control")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4_checkArgs, c(list("lmer"), l...))

  if (!is.null(list(...)[["family"]])) {
    mc[[1]] <- quote(nauf::nauf_glFormula)
    if (missing(control)) mc[["control"]] <- lme4::glmerControl()
    return(eval(mc, parent.frame()))
  }

  cstr <- "check.formula.LHS"
  lme4_checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4_checkFormulaData(formula, data,
    checkLHS = control$check.formula.LHS == "stop")
  formula <- stats::as.formula(formula, env = denv)
  lme4_RHSForm(formula) <- lme4::expandDoubleVerts(lme4_RHSForm(formula))
  mc$formula <- formula
  m <- match(c("data", "subset", "weights", "offset", "ncs_scale"),
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(nauf::nauf_model.frame)

  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- formula
  fr <- eval(mf, parent.frame())
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)

  reTrms <- nauf_mkReTrms(fr)
  wmsgNlev <- lme4_checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4_checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ", "please use ",
      shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- lme4_checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)

  mf[[1L]] <- quote(stats::model.frame)
  fixedform <- formula
  lme4_RHSForm(fixedform) <- lme4::nobars(lme4_RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr,
    "terms"), "predvars")

  ranform <- formula
  lme4_RHSForm(ranform) <- lme4::subbars(lme4_RHSForm(lme4_reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr),
    "predvars")

  X <- nauf_model.matrix(fr)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4_chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4_checkScaleX(X, kind = scaleX.chk)

  return(list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula,
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}

