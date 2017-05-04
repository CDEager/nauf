

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
#' comparisons, the \code{subset} argument can be used in the call to
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
  
  # reorder for consistency
  lvs <- lvs[vars[vars %in% names(lvs)]]
  xlev <- xlev[vars[vars %in% names(xlev)]]
  
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
#' @param subset If \code{TRUE} (the default), the \code{na.info} element of the
#'   \code{\linkS4class{nauf.ref.grid}} is used to determine what subset of the 
#'   reference grid should be used.  If not \code{TRUE}, then either a named 
#'   list describing the subset of the reference grid to be used (see 'Details' 
#'   for structure), or any value besides a list (e.g. \code{FALSE}, \code{NULL}, 
#'   \code{NA}, etc.) to indicate that the entire reference grid should be used.
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
#' @importFrom utils combn
#'
#' @export
nauf_pmmeans <- function(object, specs, pairwise = FALSE, subset = NULL,
                         na_as_level = NULL, ...) {
  if (!is.nauf.ref.grid(object)) stop("must supply a nauf.ref.grid")
  rg <- object$ref.grid

  dots <- list(...)
  if (length(dots)) warning("Ignoring arguments: ", add_quotes(names(dots)))
  
  specs <- check_specs(specs, pairwise, object$ref.grid@model.info$xlev,
    object$ref.grid@model.info$terms, na_as_level)
  subset <- grid_subset(subset, specs)
  
  if (!is.null(subset$keep)) {
    keep <- eval(subset$keep, envir = object$ref.grid@grid)
    object$ref.grid@grid <- object$ref.grid@grid[keep, , drop = FALSE]
    object$ref.grid@linfct <- object$ref.grid@linfct[keep, , drop = FALSE]
  }
  
  if (length(specs$facs)) {
    est_grid <- unique(object$ref.grid@grid[, specs$facs, drop = FALSE])
    if (!nrow(est_grid)) {
      stop("Factors are given in 'specs' but no levels exist in the specified",
        " subset")
    }
    if (!is.null(subset$est_keep)) {
      keep <- eval(subset$est_keep, envir = est_grid)
      est_grid <- est_grid[keep, , drop = FALSE]
    }
    if (!nrow(est_grid)) {
      stop("Factors are given in 'specs' but no levels exist in the specified",
        " subset")
    }
  }
  
  if (length(specs$nums)) {
    if (!length(specs$facs) || !nrow(est_grid)) {
      est_grid <- data.frame(t(rep("inc_1", length(specs$nums))))
      colnames(est_grid) <- specs$nums
    } else {
      est_grid[specs$nums] <- "inc_1"
    }
    for (v in specs$nums) {
      object$ref.grid@grid[[v]] <- object$ref.grid@grid[[v]] + 1
    }
    object$ref.grid@linfct <- model.matrix(object$ref.grid@model.info$terms,
      object$ref.grid@grid) - object$ref.grid@linfct
  }
  if (!nrow(est_grid)) stop("No grid rows satisfy subsetting conditions")
  
  if (length(specs$facs)) {
    est <- estimate_contrasts(est_grid, specs$facs)
    est <- eval(est, envir = object$ref.grid@grid)
    names(est) <- paste(1:length(est))
  } else {
    est <- list("1" = rep(1 / nrow(object$ref.grid@grid),
      nrow(object$ref.grid@grid)))
  }
  
  pmms <- list()
  
  pmms$pmmeans <- lsmeans::contrast(object$ref.grid, est)
  pmms$pmmeans@roles$predictors <- colnames(est_grid)
  pmms$pmmeans@grid <- est_grid
  for (j in 1:ncol(est_grid)) {
    pmms$pmmeans@grid[[j]] <- paste(est_grid[[j]])
  }
  pmms$pmmeans@misc$estName <- "pmmean"
  pmms$pmmeans@misc$infer <- c(TRUE, FALSE)
  pmms$pmmeans@misc$famSize <- nrow(est_grid)
  if (length(tran <- pmms$pmmeans@misc$orig.tran)) {
    pmms$pmmeans@misc$tran <- tran
  }

  if (specs$pw) {
    if (length(est) == 1) {
      specs$pw <- paste("Specified 'pairwise', but there is only one",
        "combination that satisfies subsetting conditions.")
      warning(specs$pw)
      
    } else {
      est_grid <- sapply(est_grid, as.character)
      combs <- utils::combn(length(est), 2)
      est <- list_mat_cols(apply(combs, 2,
        function(x) est[[x[1]]] - est[[x[2]]]))
      names(est) <- apply(combs, 2, function(x) paste0(
        paste(est_grid[x[1], ], collapse = ","),
        " - ",
        paste(est_grid[x[2], ], collapse = ",")
      ))
      
      pmms$contrasts <- lsmeans::contrast(object$ref.grid, est)
      pmms$contrasts@misc$estType <- "pairs"
      pmms$contrasts@misc$adjust <- "tukey"
      pmms$contrasts@misc$methDesc <- "pairwise differences"
      pmms$contrasts@misc$famSize <- nrow(est_grid)
      if (length(tran <- pmms$contrasts@misc$orig.tran)) {
        pmms$contrasts@misc$tran <- tran
      }
    }
  }

  if (length(object$ref.grid@matlevs) && any(sapply(object$ref.grid@matlevs,
  length) > 1)) {
    warning("Some matrix elements have more than one column.  If these are\n",
      "  from a call to poly(), results may be misleading.")
  }

  class(pmms) <- c("nauf.pmm", "lsm.list", "list")
  attr(pmms, "nauf.specs") <- list(
    variables = specs$vars,
    pairwise = specs$pw,
    averaged_over = setdiff(specs$avgfac, subset$cond),
    held_at_mean = specs$avgcov,
    conditioned_on = subset$cond,
    keep_NA = specs$keepNA,
    drop_NA = specs$dropNA,
    subset = subset$subset,
    note = specs$note)
  
  return(pmms)
}


check_named_list <- function(x, nms) {
  if (!is.list(x)) return("Not a list")
  
  if (!length(x)) return("Empty list")
  
  if (is.null(nx <- names(x))) return("Not a named list")
  
  if (any(sapply(x, function(u) !is.character(u) && !all(is.na(u))))) {
    return("Not all list elements are character vectors")
  }
  
  if (is.list(nms)) {
    if (length(not_in <- setdiff(nx, names(nms)))) {
      return(paste(add_quotes(not_in),
        "is/are not unordered factors in the model"))
    }
    
    not_in <- mapply(setdiff, x, nms[nx], SIMPLIFY = FALSE)
    if (length(wrong <- which(lengths(not_in) > 0))) {
      nx <- nx[wrong[1]]
      wrong <- not_in[[wrong[1]]]
      return(paste(add_quotes(wrong), "is/are not valid levels for the factor",
        add_quotes(nx)))
    }
    
  } else {
    if (length(not_in <- setdiff(nx, nms))) {
      return(paste(add_quotes(not_in),
        "is/are not unordered factors in the model"))
    }
  }
  
  if (length(nx) != length(unique(nx))) {
    return("list names are not unique")
  }
  
  return(TRUE)
}


join_and <- function(conds) {
  if (!length(conds)) return(NULL)
  conds <- conds[sapply(conds, function(x) !is.null(x))]
  if (!length(conds)) return(NULL)

  if (length(conds) > 1) {
    cl <- call("&")
    cl[[3]] <- conds[[1]]
    cl[[2]] <- conds[[2]]
    if (length(conds) > 2) {
      for (k in 3:length(conds)) {
        cl <- substitute(a & b, list(a = cl, b = conds[[k]]))
      }
    }
  } else {
    cl <- conds[[1]]
  }
  
  return(cl)
}


join_or <- function(conds) {
  if (!length(conds)) return(NULL)
  conds <- conds[sapply(conds, function(x) !is.null(x))]
  if (!length(conds)) return(NULL)

  if (length(conds) > 1) {
    cl <- call("|")
    cl[[3]] <- conds[[1]]
    cl[[2]] <- conds[[2]]
    if (length(conds) > 2) {
      for (k in 3:length(conds)) {
        cl <- substitute(a | b, list(a = cl, b = conds[[k]]))
      }
    }
  } else {
    cl <- conds[[1]]
  }
  
  return(cl)
}


grid_subset <- function(subset, specs) {
  if (length(specs$dropNA)) {
    est_keep <- join_and(lapply(specs$dropNA, function(fac) 
      substitute(!is.na(F), list(F = as.name(fac)))))
  } else {
    est_keep <- NULL
  }
  
  if (is.null(subset)) {
    return(list(keep = NULL, est_keep = est_keep, cond = character(),
      subset = NULL))
  }
  
  if (!is.list(subset) || !length(subset)) {
    stop("'subset' must be a list specifying groups")
  }
  if (!any(sapply(subset, is.list))) subset <- list(subset)
  
  msg <- mapply(check_named_list, subset, MoreArgs = list(nms = specs$uflev),
    SIMPLIFY = FALSE)
  msg <- msg[sapply(msg, function(u) !isTRUE(u))]
  if (length(msg)) {
    stop("'subset' not properly specified:\n  ", msg[[1]])
  }
  
  cond <- unname(unique(unlist(lapply(subset, names))))
  
  keep <- join_or(lapply(subset, function(group) join_and(mapply(
    function(fac, levs) substitute(F %in% L, list(F = as.name(fac), L = levs)),
    names(group), group, SIMPLIFY = FALSE))))
  
  return(list(keep = keep, est_keep = est_keep, cond = cond, subset = subset))
}


check_specs <- function(specs, pw, xlev, mt, keepNA) {
  vv <- varnms(mt)
  
  if (inherits(specs, "formula")) {
    resp <- has_resp(specs <- stats::terms(specs))
    specs <- varnms(specs)
    if (is.null(specs)) stop("'specs' has no variables")
    
    if (pw <- (resp && specs[1] == "pairwise")) {
      specs <- specs[-1]
    } else if (resp) {
      stop("If 'specs' has a left hand side, it must be 'pairwise'")
    }
  }
  if (!is.character(specs)) {
    stop("'specs' should be a character vector or formula specifying variables")
  }
  if (!length(specs)) stop("'specs' has no variables")
  if (length(not_in <- setdiff(specs, vv))) {
    stop("The following variables are not valid:\n  ", add_quotes(not_in),
      "\nValid variables are:\n  ", add_quotes(vv))
  }
  
  # reorder for consistency; gets rid of possible duplicates
  specs <- vv[vv %in% specs]
  
  all_facs <- names(xlev)
  # intersect since info$uf may contain facs in ranef but not fixef
  all_uf <- names(uflev <- xlev[intersect(names(xlev), names(attr(mt,
    "nauf.info")$uf))])
  all_nauf <- all_facs[sapply(xlev, anyNA)]
  
  facs <- intersect(specs, all_facs)
  avgfac <- setdiff(all_facs, facs)
  uf <- intersect(facs, all_uf)
  nauf <- intersect(uf, all_nauf)
  nums <- setdiff(specs, all_facs)
  avgcov <- setdiff(vv, c(all_facs, nums))
  
  note <- character()
  fmat <- attr(mt, "factors")
  specs_which <- which(colnames(fmat) == (specs_term <- paste(specs,
    collapse = ":")))
  if (!length(specs_which)) {
    note <- paste0("The interaction term '", specs_term,
      "' is not in the model.")
    
  } else {
    ord <- attr(mt, "order")
    specs_order <- ord[specs_which]
    has_specs <- sapply(fmat[specs, , drop = FALSE], all)
    ord[!has_specs] <- 0
    if (length(higher <- which(ord > specs_order))) {
      note <- paste(add_quotes(specs_term), "is included in higher order",
        "interaction(s)", add_quotes(colnames(fmat)[higher]))
    }
  }
  
  if (is.null(keepNA)) keepNA <- character()
  arg_keepNA <- keepNA
  keepNA <- intersect(keepNA, nauf)
  if (length(dropped <- setdiff(arg_keepNA, keepNA))) {
    warning("Dropped from 'keepNA' (not unordered factors with NA values ",
      "included in 'specs':\n  ", add_quotes(dropped))
  }
  dropNA <- setdiff(nauf, keepNA)
  
  
  return(list(vars = specs, facs = facs, nums = nums, uf = uf, uflev = uflev,
    pw = pw, keepNA = keepNA, avgfac = avgfac, avgcov = avgcov, nauf = nauf,
    dropNA = dropNA, note = note))
}


estimate_contrasts <- function(est_grid, facs) {
  est <- call("list")
  est_grid <- as.data.frame(lapply(est_grid[facs], as.character),
    stringsAsFactors = FALSE)
  for (i in 1:nrow(est_grid)) {
    est[[i + 1]] <- call("as_simplex")
    est[[i + 1]][[2]] <- join_and(mapply(function(fac, levs)
      substitute(F %in% L, list(F = as.name(fac), L = levs)),
      facs, est_grid[i, facs, drop = FALSE], SIMPLIFY = FALSE))
  }
  return(est)
}

