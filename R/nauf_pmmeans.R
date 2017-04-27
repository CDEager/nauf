

#' Create a reference grid for a \code{nauf} model.
#'
#' \code{nauf_ref.grid} creates a reference grid for a model fit with
#' \code{\link{nauf_contrasts}} which can be used to obtain predicted marginal
#' means and pairwise comparisons for factors with \code{\link{nauf_pmmeans}}.
#' See \code{\link[lsmeans]{ref.grid}} for details of reference grid properties.
#'
#' A reference grid creates a data frame which contains all possible
#' combinations of the factors in a regression model, holding all covariates
#' at their mean values.  There are many options for
#' \code{\link[lsmeans]{ref.grid}} which are not currently supported for
#' \code{nauf} models.  The main functionality which is not currently supported
#' is that the reference grid cannot be created specifying certain levels for
#' variables (i.e. the \code{at} argument; but see the various subsetting options
#' in \code{\link{nauf_pmmeans}}).  For models with polynomial covariates,
#' the results may be misleading, as is usually the case.  A direct call to
#' \code{\link[lsmeans]{ref.grid}} will result in warnings (or possibly errors),
#' and inference made with the resulting object will be misleading and/or
#' incorrect.  Only \code{nauf_ref.grid} should be used.
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
#' comparisons, subsetting arguments can be used in the call to
#' \code{\link{nauf_pmmeans}} to ensure only the correct subsets are considered.
#'
#' @param mod Any \code{nauf} regression model.
#' @param KR Only applies when \code{mod} is a \code{\linkS4class{nauf.lmerMod}}
#'   fit with \code{REML = TRUE}. If \code{KR = TRUE}, then the Kenward-Roger
#'   approximation is used to calculate degrees of freedom.  If 
#'   \code{KR = FALSE}, then the Satterthwaite approximation is used.  If
#'   \code{KR = NULL} (the default), then the Kenward-Roger approximation is
#'   used when there are \code{3000} observations or fewer in the model, and the
#'   Satterthwaite approximation is used when there are more than \code{3000}
#'   observations (since in this case the the Kenward-Roger approximation
#'   requires a lot of computation time).  When \code{mod} is a
#'   \code{\linkS4class{nauf.lmerMod}} fit with \code{REML = FALSE}, then the
#'   Satterthwaite approximation is always used.  The Kenward-Roger method
#'   is implemented with \code{\link[pbkrtest]{Lb_ddf}} and the Satterthwaite
#'   method is implemented with \code{\link[lmerTest]{calcSatterth}}.
#'
#' @return A \code{\linkS4class{nauf.ref.grid}}.
#'
#' @examples
#' \dontrun{
#' mod <- nauf_lm(y ~ f1 + f2 + x)
#' rg <- nauf_ref.grid(mod)
#' }
#'
#' @seealso \code{\link{nauf_pmmeans}}
#'
#' @importFrom pbkrtest Lb_ddf vcovAdj
#' @importFrom lsmeans ref.grid
#' @importFrom lmerTest calcSatterth
#'
#' @export
nauf_ref.grid <- function(mod, KR = NULL) {
  if (!is.nauf.model(mod)) stop("Must supply a nauf model")

  fenr <- stats::delete.response(terms(mod))
  info <- attr(fenr, "nauf.info")
  vars <- varnms(fenr)
  mf <- model.frame(mod)
  
  xlev <- lvs <- mlvs <- list()
  for (v in vars[vars %in% names(info$uf)]) {
    lvs[[v]] <- info$uf[[v]][[1]]
    if (info$hasna[v]) {
      lvs[[v]][length(lvs[[v]]) + 1] <- NA
      mf[[v]] <- addNA(mf[[v]])
    }
    xlev[[v]] <- lvs[[v]]
  }
  for (v in vars[vars %in% names(info$of)]) {
    xlev[[v]] <- lvs[[v]] <- info$of[[v]]$levels
  }
  for (v in vars[vars %in% names(info$num)]) {
    lvs[[v]] <- info$num[[v]]
  }
  for (v in vars[vars %in% names(info$mat)]) {
    mlvs[[v]] <- info$mat[[v]]
  }
  
  g <- expand.grid(lvs)
  g[names(mlvs)] <- lapply(mlvs, function(m) {
    matrix(m, nrow(g), length(m), byrow = TRUE)
  })
  
  mm <- model.matrix(fenr, g)
  asgn <- attr(mm, "assign")
  
  g[[".wgt."]] <- data.frame(xtabs(~ ., mf[names(xlev)]))$Freq
  
  summ <- summary(mod)
  
  rg <- list()
  class(rg) <- c("nauf.ref.grid", "list")
  temp_dat <- data.frame(y = 1:5, x = c(2, 5, 3, 6, 1))
  temp_mod <- lm(y ~ x, temp_dat)
  rg$ref.grid <- lsmeans::ref.grid(temp_mod, data = temp_dat)

  rg$ref.grid@model.info <- list(
    call = summ$call,
    terms = fenr,
    xlev = xlev)
  rg$ref.grid@roles <- list(
    predictors = vars,
    responses = character(),
    multresp = character())
  rg$ref.grid@grid <- g
  rg$ref.grid@levels <- lvs
  rg$ref.grid@matlevs <- mlvs
  rg$ref.grid@linfct <- mm
  rg$ref.grid@bhat <- summ$coefficients[, 1]
  rg$ref.grid@nbasis <- matrix()
  if (is.nauf.merMod(mod)) {
    rg$ref.grid@V <- as.matrix(summ$vcov)
  } else {
    rg$ref.grid@V <- vcov(mod)
  }
  if (is.nauf.lmerMod(mod)) {
    if (!lme4::isREML(mod)) KR <- FALSE
    if (is.null(KR)) KR <- nobs(mod) <= 3000
    if (KR) {
      rg$ref.grid@dffun <- function(k, dfargs) {
        pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
      }
      rg$ref.grid@dfargs <- list(
        unadjV = vcov(mod),
        adjV = pbkrtest::vcovAdj(mod))
    } else {
      rg$ref.grid@dffun <- function(k, dfargs) {
        lmerTest::calcSatterth(dfargs$object, k)$denom
      }
      rg$ref.grid@dfargs <- list(object = mod)
    }
  } else if (is.nauf.glmerMod(mod)) {
    rg$ref.grid@dffun <- function(k, dfargs) NA
    rg$ref.grid@dfargs <- list()
  } else {
    rg$ref.grid@dffun <- function(k, dfargs) dfargs$df
    rg$ref.grid@dfargs <- list(df = mod$df)
  }
  rg$ref.grid@misc <- list(
    estName = "prediction",
    estType = "prediction",
    infer = c(FALSE, FALSE),
    level = 0.95,
    adjust = "none",
    famSize = ncol(g) - 1,
    avgd.over = character(),
    assign = asgn)
  rg$ref.grid@post.beta <- matrix()
  if (!is.linear(family <- get_family(mod))) {
    rg$ref.grid@misc$tran <- family$link
    family <- family$family
    rate <- length(grep("poisson", family)) + length(grep("Negative", family))
    prob <- length(grep("binomial", family))
    if (rate) {
      rg$ref.grid@misc$inv.lbl <- "rate"
    } else if (prob) {
      rg$ref.grid@misc$inv.lbl <- "prob"
    } else {
      rg$ref.grid@misc$inv.lbl <- "response"
    }
  }
  
  return(rg)
}


#' Predicted marginal means for \code{nauf} models.
#' 
#' \code{nauf_pmmeans} computes predicted marginal means (also called
#' least-squares means, a term which is less generally applicable) for a factor
#' or an interaction involving factors in a \code{nauf} regression model,
#' and optionally pairwise comparisons of these group means.  There are several
#' subsetting arguments which allow for the predictions to be generated over
#' only the relevant subset of the reference grid created by
#' \code{\link{nauf_ref.grid}}.  The function also allows covariates to be
#' included. See 'Subsetting', 'Factors', 'Covariates', and
#' 'Factor-Covariate Interactions' for details.
#'
#' @section Subsetting:
#' Because, by design, certain factors in a regression fit with
#' \code{\link{nauf_contrasts}} are applicable only within subsets of the data,
#' using \code{\link[lsmeans]{lsmeans}} directly will generate and average
#' over estimates which are invalid.  \code{nauf_pmmeans} allows the user
#' to generate predicted marginal means for only the subsets of the
#' data where the estimates make sense.  This is done through four subsetting
#' arguments, any combination of which can be used to specify the subset of
#' the reference grid over which to generate estimates:
#'
#' \describe{
#'   \item{keep_level}{A named list where each element is a character vector of
#'     factor levels whose name is the name of the factor.  For example, 
#'     \code{keep_group = list(f1 = c("A", "B"), f2 = "C")} would indicate that 
#'     only the subset of the data where the factor \code{f1} is either \code{A} 
#'     or \code{B} and the factor \code{f2} is \code{C} should be considered.}
#'   \item{drop_level}{A named list with the same structure as \code{keep_level},
#'     but indicating factor levels which should \emph{not} be considered.  For 
#'     example, \code{drop_level = list(f1 = "A")} would indicate that rows in 
#'     the reference grid where the factor \code{f1} is \code{A} should not be 
#'     considered.}
#'   \item{keep_group}{A list whose elements are named lists with the same 
#'     structure as \code{keep_level} and \code{drop_level}.  For example, 
#'     \code{keep_group = list(list(f1 = "A", f2 = c("C", "D")), list(f1 = "B",
#'     f2 = c("D", "E")))}.  For each named list in \code{keep_group}, 
#'     \code{\link[base]{expand.grid}} is called.  Any row in the reference grid
#'     which pertains to one of the resulting groups is kept.  In this example, 
#'     using the format \code{f1:f2}, this would indicate that the groups 
#'     \code{A:C}, \code{A:D}, \code{B:D}, and \code{B:E} should be kept, and 
#'     all other combinations of these factors dropped.  Note that this 
#'     criterion applies to all of the indicated groups simultaneously.  If a 
#'     row of the reference grid belongs to a group specified in the first list 
#'     element of \code{keep_group} but not the second list element, the row 
#'     would still be kept.}
#'   \item{drop_group}{A list with the same structure as \code{keep_group}.  Any
#'     row of the reference grid which belongs to at least one of the groups 
#'     specified in \code{drop_group} is not considered.}
#' }
#'
#' When all four subsetting arguments are empty lists (the default), all cells
#' in the reference grid are used (which makes sense for covariates, and for 
#' factors that are applicable with the same meaning throughout the dataset).  
#' For factors which contain \code{NA} values
#' in the original data frame, \code{NA} is treated as a level and should be
#' specified in the subsetting arguments without quotes.  The subsetting
#' arguments can also be used jointly in any combination.  Only observations
#' which meet the criteria of all subsetting arguments are used.
#'
#' To give an example, consider the \code{\link{plosives}} dataset, and assume
#' we are interested in the model \code{cdur ~ dialect + spont}.
#' The \code{dialect} factor has levels \code{Cuzco}, \code{Lima}, and
#' \code{Valladolid}, and the \code{spont} variable is a logical indicating
#' whether the observation pertains to spontaneous (\code{TRUE}) or read
#' (\code{FALSE}) speech.  However, for the Valladolid speakers, there is
#' only spontaneous speech data, so we code \code{spont} as \code{NA} when
#' \code{dialect = Valladolid}.  The reference grid for this model (fit with
#' \code{\link{nauf_lm}}) will contain the combinations
#' \code{dialect = Cuzco, spont = NA}; \code{dialect = Lima, spont = NA};
#' \code{dialect = Valladolid, spont = TRUE}; and
#' \code{dialect = Valladolid, spont = FALSE}, even though these are not
#' real groups.  If we use no subsetting arguments and call \code{nauf_pmmeans}
#' with \code{specs = pairwise ~ dialect}, then the predicted marginal means
#' will be incorrect.  We want to compare the dialects conditioning on
#' spontaneous speech; i.e. we want to compare
#' \code{dialect = Cuzco, spont = FALSE}; \code{dialect = Lima, spont = FALSE};
#' and \code{dialect = Valladolid, spont = NA}. To do this, we could use either
#' of the follwing subsetting arguments (where \code{rg} is the 
#' \code{\linkS4class{nauf.ref.grid}} for the model):
#'
#' \code{\preformatted{
#' # possibility 1
#' nauf_pmmeans(rg, pairwise ~ dialect, drop_group = list(
#'   list(dialect = "Valladolid", spont = c("TRUE", "FALSE")),
#'   list(dialect = c("Cuzco", "Lima"), spont = c("FALSE", NA))))
#'
#' # possbility 2
#' nauf_pmmeans(rg, pairwise ~ dialect, keep_group = list(
#'   list(dialect = "Valladolid", spont = NA),
#'   list(dialect = c("Cuzco", "Lima"), spont = "TRUE")))
#' }}
#'
#' Both of these possibilities will result in the same three predicted
#' marginal means in the \code{pmmeans} element, with pairwise comparisons
#' in the \code{contrasts} element of the returned \code{lsm.list}.  If we are
#' interested in the effect of \code{spont}, then the situation is substantially 
#' simpler. We just need to drop the Valladolid dialect and the \code{NA} values
#' for \code{spont}, which we could do with
#' \code{drop_level = list(dialect = "Valladolid", spont = NA)}.  See the
#' 'Examples' section.
#'
#' @section Factors:
#' When all variables in \code{specs} are factors, the full interaction term
#' is taken, and combinations of factor levels which are excluded by the
#' the subsetting arguments are dropped.  Predicted marginal means are then
#' computed as they are with \code{\link[lsmeans]{lsmeans}}, averaging over
#' only the relevant subset of the reference grid, and pairwise comparisons
#' for these means are computed as they are normally computed (if \code{specs}
#' has \code{pairwise} on the left hand side).
#'
#' @section Covariates:
#' If \code{specs} indicates a single covariate, then the effect of an increase
#' of \code{1} in the value of the covariate is calculated as the \code{pmmean}.
#' If \code{specs} indicates multiple covariates, then the effect of a
#' simultaneous increase of \code{1} in the value of all of the covariates is
#' calculated as the \code{pmmean}.  In this case, pairwise comparisons are not
#' possible, as there is always one estimate.
#'
#' @section Factor-Covariate Interactions:
#' If \code{specs} indicates a combination of factors and covariates, then
#' the effect of a simultaneous increase of \code{1} in the covariates for
#' each level of the interaction term of the factors is calculated.  If
#' \code{pairwise} is
#' used as the left hand side of \code{specs}, then pairwise comparisons for
#' the effect of the simultaneous increase of \code{1} in the covariates for
#' the different factor interaction levels are computed.
#'
#' @param object A \code{\linkS4class{nauf.ref.grid}} or any \code{nauf} model
#'   (see \code{\link{is.nauf.model}}).  When a regression model is provided,
#'   \code{\link{nauf_ref.grid}} is called first, and the result is used.
#'   Creating the \code{\linkS4class{nauf.ref.grid}} can be time-consuming for
#'   for mixed effects models, and so it is strongly suggested that a
#'   \code{\linkS4class{nauf.ref.grid}} be created first, so that it can be
#'   used for multiple calls to \code{nauf_pmmeans}.
#' @param specs A formula whose right hand side indicates variables for which 
#'   predicted marginal means are to be computed.  Regardless of whether 
#'   \code{*}, \code{:}, or \code{+} is used to separate the variables, the 
#'   full interaction term is considered.  If the formula has a left hand side, 
#'   it must be the key word \code{pairwise}, indicating that pairwise 
#'   comparisons should be computed in addition to the predicted marginal means.
#' @param keep_level A named list of character vectors indicating factor levels 
#'   in the reference grid which should be considered.  See 'Subsetting'.
#' @param drop_level A named list of character vectors indicating factor levels 
#'   in the reference grid which should not be considered.  See 'Subsetting'.
#' @param keep_group A list whose elements are named lists of character vectors 
#'   defining combinations of factor levels in the reference which should be 
#'   considered.  See 'Subsetting'.
#' @param drop_group A list whose elements are named lists of character vectors 
#'   defining combinations of factor levels in the reference which should not be 
#'   considered.  See 'Subsetting'.
#' @param ... Additional arguments. Currently unused and ignored with a warning.
#'
#' @return A \code{lsm.list} containing an element \code{pmmeans} and, if
#'   pairwise comparisons are made, a second element \code{contrasts}, both
#'   of which are \code{lsmobj} objects (see 
#'   \code{\link[lsmeans]{ref.grid-class}}).
#'
#' @examples
#' dat <- droplevels(subset(plosives, voicing == "Voiceless"))
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' sdat <- standardize(cdur ~ dialect + spont, dat)
#' mod <- nauf_lm(sdat$formula, sdat$data)
#' rg <- nauf_ref.grid(mod)
#'
#' ## incorrect estimates
#' # averages over non-existent groups
#' dialect_wrong <- nauf_pmmeans(rg, pairwise ~ dialect)
#' 
#' # treats spont = NA values as a level and generates 3 pmmeans rather than 2
#' spont_wrong <- nauf_pmmeans(rg, pairwise ~ spont)
#'
#' ## correct estimates
#' dialect1 <- nauf_pmmeans(rg, pairwise ~ dialect, drop_group = list(
#'   list(dialect = "Valladolid", spont = c("TRUE", "FALSE")),
#'   list(dialect = c("Cuzco", "Lima"), spont = c("FALSE", NA))))
#'
#' dialect2 <- nauf_pmmeans(rg, pairwise ~ dialect, keep_group = list(
#'   list(dialect = "Valladolid", spont = NA),
#'   list(dialect = c("Cuzco", "Lima"), spont = "TRUE")))
#'
#' spont <- nauf_pmmeans(rg, pairwise ~ spont,
#'   drop_level = list(dialect = "Valladolid", spont = NA))
#'
#' @seealso \code{\link{nauf_ref.grid}}, \code{\link[lsmeans]{lsmeans}}, and
#'   \code{\link[lsmeans]{ref.grid-class}}.
#'
#' @importFrom lsmeans contrast lsmeans pmmeans
#'
#' @export
nauf_pmmeans <- function(object, specs, keep_level = list(),
                         drop_level = list(), keep_group = list(),
                         drop_group = list(), ...) {
  dots <- list(...)
  if (length(dots)) warning("Ignoring arguments: ", add_quote(names(dots)))
  
  if (!is.nauf.ref.grid(object)) object <- nauf_ref.grid(object)
  rg <- object$ref.grid
  
  if (!inherits(specs, "formula")) stop("'specs' must be a formula")
  resp <- has_resp(specs <- stats::terms(specs))
  specs <- varnms(specs)
  if (is.null(specs)) stop("'specs' has no variables")
  if (pw <- (resp && specs[1] == "pairwise")) {
    specs <- specs[-1]
  } else if (resp) {
    stop("If 'specs' has a left hand side, it must be 'pairwise'")
  }
  if (!length(specs)) stop("'specs' has no variables")
  if (length(not_in <- specs[!(specs %in% rg@roles$predictors)])) {
    stop("The following variables are not valid:\n  ", add_quotes(not_in),
      "Valid variables are:\n  ", add_quotes(rg@roles$predictors))
  }
  
  all_facs <- names(xlev <- rg@model.info$xlev)
  vv <- rownames(fmat <- attr((mt <- rg@model.info$terms), "factors"))
  nums <- specs[!(specs %in% all_facs)]
  facs <- specs[specs %in% all_facs]
  
  specs <- vv[vv %in% specs]
  specs_which <- which(colnames(fmat) == (specs_term <- paste(specs,
    collapse = ":")))
  if (!length(specs_which)) {
    warning("The interaction term ", specs_term, " is not in the model.")
  } else {
    specs_order <- attr(mt, "order")[specs_which]
    o_has_specs <- attr(mt, "order")[sapply(fmat[specs, , drop = FALSE], all)]
    higher <- which(o_has_specs > specs_order)
    if (length(higher)) {
      warning("Results may be misleading due to higher order interactions.")
    }
  }
  
  
  check_subset(xlev, keep_level, drop_level, keep_group, drop_group)
  keep <- grid_subset(rg@grid, keep_level, drop_level, keep_group,
    drop_group)
  rg@grid <- rg@grid[keep, , drop = FALSE]
  rg@linfct <- rg@linfct[keep, , drop = FALSE]
  
  est_lvs <- xlev[facs]
  if (length(nums)) {
    g <- rg@grid[, -ncol(rg@grid), drop = FALSE]
    for (v in nums) {
      g[[v]] <- g[[v]] + 1
      est_lvs[[v]] <- "inc_1"
    }
    rg@linfct <- model.matrix(rg@model.info$terms, g) - rg@linfct
  }
  est_grid <- expand.grid(est_lvs)
  if (length(facs)) {
    est_grid <- est_grid[sub_est_grid(rg@grid, est_grid, facs), , drop = FALSE]
  }
  if (!nrow(est_grid)) {
    stop("No possible combinations satisfy subsetting conditions")
  }
  
  pmms <- list()
  est <- estimate_contrasts(rg@grid, est_grid, facs)
  pmms$pmmeans <- marginal_means(rg, est, est_grid)
  if (pw) {
    if (length(est) == 1) {
      pw <- paste("Specified 'pairwise', but there is only one combination",
        "that satisfies subsetting conditions.  Contrasts not computed")
      warning(pw)
    } else {
      pmms$contrasts <- pairwise_contrasts(rg, est, est_grid)
    }
  }
  
  if (length(rg@matlevs) && any(sapply(rg@matlevs, length) > 1)) {
    warning("Some matrix elements have more than one column.  If these are\n",
      "  from a call to poly(), results may be misleading.")
  }

  class(pmms) <- c("lsm.list", "list")
  attr(pmms, "nauf.specs") <- list(specs = specs, keep_level = keep_level,
    drop_level = drop_level, drop_group = drop_group, pairwise = pw)

  return(pmms)
}


check_subset <- function(xlev, keep_level, drop_level, keep_group,
                              drop_group) {
  if (!length(xlev)) {
    if (!isTRUE(all.equal(c(keep_level, drop_level, keep_group, drop_group),
    list()))) {
      stop("Subsetting arguments are not empty lists, but there are no factors",
        " in the model")
    }
  } else {
  
    all_facs <- names(xlev)
    
    validate_level <- function(u) {
      if (!is.list(u)) return("Not a list")
      if (!length(u)) return(TRUE)
      nn <- names(u)
      not_in <- names(u)[!(names(u) %in% all_facs)]
      if (length(not_in)) {
        return(paste0("The following are not factors in the model:\n  ",
          add_quotes(not_in)))
      }
      u <- lapply(nn, function(k) u[[k]][!(u[[k]] %in% xlev[[k]])])
      names(u) <- nn
      u <- u[sapply(u, function(k) length(k) > 0)]
      if (!length(u)) return(TRUE)
      return(paste0("The following are not valid levels for '", names(u)[1],
        "':\n    ", add_quotes(u[[1]])))
    }
    
    validate_group <- function(g) {
      if (!is.list(g)) return("Not a list")
      if (!length(g)) return(TRUE)
      check <- sapply(g, is.list)
      if (any(!check)) return("Contains elements that are not lists")
      check <- sapply(g, length)
      if (any(!check)) return("Contains empty list elments")
      check <- lapply(g, validate_level)
      check <- check[sapply(check, function(k) !isTRUE(k))]
      if (!length(check)) return(TRUE)
      return(check[[1]])
    }
    
    if (!isTRUE(msg <- validate_level(keep_level))) {
      stop("Invalid 'keep_level'.  ", msg)
    }
    if (!isTRUE(msg <- validate_level(drop_level))) {
      stop("Invalid 'drop_level'.  ", msg)
    }
    if (!isTRUE(msg <- validate_group(keep_group))) {
      stop("Invalid 'keep_group'.  ", msg)
    }
    if (!isTRUE(msg <- validate_group(drop_group))) {
      stop("Invalid 'drop_group'.  ", msg)
    }
  }
  
  invisible(NULL)
}


grid_subset <- function(grid, keep_level, drop_level, keep_group,
                            drop_group) {
  keep <- matrix(TRUE, nrow(grid), 4)
  
  if (length(keep_level)) {
    keep[, 1] <- match_row(grid, keep_level)
  }
  
  if (length(drop_level)) {
    keep[, 2] <- !match_row(grid, drop_level, f = any)
  }
  
  if (length(keep_group)) {
    keep[, 3] <- apply(sapply(keep_group, function(g) {
      match_row(grid, g)
    }), 1, any)
  }
  
  if (length(drop_group)) {
    keep[, 4] <- !apply(sapply(drop_group, function(g) {
      match_row(grid, g)
    }), 1, any)
  }
  
  return(apply(keep, 1, all))
}


sub_est_grid <- function(grid, est_grid, facs) {
  u <- unique(grid[, facs, drop = FALSE])
  if (nrow(u) == 1) {
    return(match_row(est_grid, u))
  }
  return(apply(apply(u, 1, function(i) {
    match_row(est_grid, data.frame(t(i)))
  }), 1, any))
}


estimate_contrasts <- function(grid, est_grid, facs) {
  if (!length(facs)) {
    return(list("1" = rep(1 / nrow(grid), nrow(grid))))
  }
  return(list_mat_cols(apply(apply(est_grid, 1,
    function(i) match_row(grid, i, facs)), 2, as_simplex)))
}


marginal_means <- function(rg, est, est_grid) {
  est <- lsmeans::contrast(rg, est)
  est@roles$predictors <- colnames(est_grid)
  for (j in 1:ncol(est_grid)) est_grid[[j]] <- paste(est_grid[[j]])
  est@grid <- est_grid
  est@misc$estName <- "pmmean"
  est@misc$infer <- c(TRUE, FALSE)
  est@misc$famSize <- nrow(est_grid)
  if (length(tran <- est@misc$orig.tran)) est@misc$tran <- tran
  return(est)
}


#' @importFrom utils combn
pairwise_contrasts <- function(rg, est, est_grid) {
  est_grid <- sapply(est_grid, as.character)
  combs <- utils::combn(length(est), 2)
  contr <- list_mat_cols(apply(combs, 2, function(x) est[[x[1]]] - est[[x[2]]]))
  names(contr) <- apply(combs, 2, function(x) paste0(
    paste(est_grid[x[1], ], collapse = ","),
    " - ",
    paste(est_grid[x[2], ], collapse = ",")
  ))
  contr <- lsmeans::contrast(rg, contr)
  contr@misc$estType <- "pairs"
  contr@misc$adjust <- "tukey"
  contr@misc$methDesc <- "pairwise differences"
  contr@misc$famSize <- nrow(est_grid)
  if (length(tran <- contr@misc$orig.tran)) contr@misc$tran <- tran
  return(contr)
}


match_row <- function(d, lvs, nms = names(lvs), f = all) {
  if (nrow(d) == 1) {
    return(f(sapply(nms, function(v) {
      d[[v]] %in% lvs[[v]]
    })))
  }
  return(apply(sapply(nms, function(v) {
    d[[v]] %in% lvs[[v]]
  }), 1, f))
}


as_simplex <- function(x) {
  return(x / sum(x))
}


list_mat_cols <- function(x) {
  return(split(x, c(col(x))))
}

